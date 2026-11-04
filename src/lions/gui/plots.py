"""Live pyqtgraph plots for LIONS, with persistent zoom/pan and a crosshair
cursor, plus the Plotly 3D view.
"""

import importlib

import numpy as np
import tempfile
import os
import pyqtgraph as pg
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QDoubleSpinBox,
    QCheckBox,
    QPushButton,
    QSpinBox,
    QDialog,
    QFileDialog,
)
from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtGui import QFont, QColor

# OpenGL is optional; it only gates the 3D toolbar action.
try:
    importlib.import_module('pyqtgraph.opengl')
    HAS_GL = True
except ImportError:
    HAS_GL = False

try:
    import plotly.graph_objects as _plotly_go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

try:
    from PyQt5.QtWebEngineWidgets import QWebEngineView

    HAS_WEBENGINE = True
except ImportError:
    HAS_WEBENGINE = False


class IntegerAxis(pg.AxisItem):
    """Axis with integer tick labels."""

    def tickStrings(self, values, scale, spacing):
        return [f'{int(round(v))}' for v in values]


BG_COLOR = '#f0f2f5'
PLOT_BG = '#ffffff'
GRID_COLOR = '#c8d0dc'
TEXT_COLOR = '#2c3e50'
TITLE_COLOR = '#1a1a2e'

FAMILY_COLORS = [
    '#1565C0',
    '#d32f2f',
    '#2e7d32',
    '#ff6f00',
    '#7b1fa2',
    '#00838f',
    '#c62828',
    '#283593',
    '#558b2f',
    '#4e342e',
]

# Antialiasing off: too slow for live updates.
pg.setConfigOptions(
    antialias=False, background=BG_COLOR, foreground=TEXT_COLOR
)


class CrosshairOverlay:
    """Crosshair cursor with a coordinate readout on a PlotItem."""

    def __init__(self, plot_item):
        self.plot = plot_item
        self.vb = plot_item.getViewBox()

        pen = pg.mkPen(color='#888888', width=0.5, style=Qt.DashLine)
        self.v_line = pg.InfiniteLine(angle=90, movable=False, pen=pen)
        self.h_line = pg.InfiniteLine(angle=0, movable=False, pen=pen)

        self.label = pg.TextItem(
            text='',
            color=TEXT_COLOR,
            anchor=(0, 1),
            fill=pg.mkBrush(255, 255, 255, 180),
        )
        self.label.setFont(QFont('sans-serif', 8))

        self.v_line.setVisible(False)
        self.h_line.setVisible(False)
        self.label.setVisible(False)

        self.plot.addItem(self.v_line, ignoreBounds=True)
        self.plot.addItem(self.h_line, ignoreBounds=True)
        self.plot.addItem(self.label, ignoreBounds=True)

        self.proxy = pg.SignalProxy(
            self.plot.scene().sigMouseMoved,
            rateLimit=30,
            slot=self._on_mouse_moved,
        )

    def _on_mouse_moved(self, evt):
        pos = evt[0]
        if self.plot.sceneBoundingRect().contains(pos):
            mouse_point = self.vb.mapSceneToView(pos)
            x = mouse_point.x()
            y = mouse_point.y()
            self.v_line.setPos(x)
            self.h_line.setPos(y)
            self.label.setText(f'x={x:.2f}\ny={y:.2f}')
            self.label.setPos(x, y)
            self.v_line.setVisible(True)
            self.h_line.setVisible(True)
            self.label.setVisible(True)
        else:
            self.v_line.setVisible(False)
            self.h_line.setVisible(False)
            self.label.setVisible(False)


class PlotCanvas(QWidget):
    """pyqtgraph canvas: six plots, each with the reference particle and a
    scatter of all particles.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.graphics_widget = pg.GraphicsLayoutWidget()
        self.graphics_widget.setBackground(BG_COLOR)
        layout.addWidget(self.graphics_widget)

        # Radius limits (mm) for XY, R-theta and R-dphi; 0 means no
        # limit.
        self.rmax_display = 0.0
        self.rmax_rtheta = 0.0
        self.rmax_rdph = 0.0
        # dphi window for R-dphi and dphi vs turn; (0, 0) means no
        # filter.
        self.dphi_min = 0.0
        self.dphi_max = 0.0

        self.plot_xy = self.graphics_widget.addPlot(row=0, col=0)
        self.plot_energy = self.graphics_widget.addPlot(
            row=0,
            col=1,
            axisItems={'bottom': IntegerAxis(orientation='bottom')},
        )
        self.plot_rtheta = self.graphics_widget.addPlot(row=1, col=0)
        self.plot_zpz = self.graphics_widget.addPlot(row=1, col=1)
        self.plot_rdph = self.graphics_widget.addPlot(row=2, col=0)
        self.plot_dphw = self.graphics_widget.addPlot(
            row=2,
            col=1,
            axisItems={'bottom': IntegerAxis(orientation='bottom')},
        )

        self._style_plot(self.plot_xy, 'XY trajectory', 'X (mm)', 'Y (mm)')
        self._style_plot(
            self.plot_energy, 'Energy vs turn', 'Turn', 'Energy (MeV/A)'
        )
        self._style_plot(
            self.plot_rtheta, 'R - \u03b8', '\u03b8 (deg)', 'R (mm)'
        )
        self._style_plot(
            self.plot_zpz, 'Z - Pz phase space', 'Z (mm)', 'Pz (mrad)'
        )
        self._style_plot(
            self.plot_rdph, 'R - d\u03c6', 'R (mm)', 'd\u03c6 (deg)'
        )
        self._style_plot(
            self.plot_dphw, 'd\u03c6 vs turn', 'Turn', 'd\u03c6 (deg)'
        )

        self.plot_xy.setAspectLocked(True, ratio=1.0)

        # Legends are added on the first family when nbrfam > 1.
        self._legends_added = False

        self._deflector_items = []

        self._fam_data = {}  # {k: {'x':[], 'y':[], 'turn':[], ...}}
        # {k: {'xy':..., 'energy':..., 'rtheta':..., 'zpz':...}}
        self._fam_curves = {}

        self.crosshairs = []
        for plot in [
            self.plot_xy,
            self.plot_energy,
            self.plot_rtheta,
            self.plot_zpz,
            self.plot_rdph,
            self.plot_dphw,
        ]:
            self.crosshairs.append(CrosshairOverlay(plot))

        # app.py calls canvas.fig.savefig(), the matplotlib API.
        self.fig = _FigShim(self)

    @staticmethod
    def _style_plot(plot_item, title, xlabel, ylabel):
        """Apply the common style to a PlotItem."""
        plot_item.setTitle(title, color=TITLE_COLOR, size='11pt', bold=True)
        plot_item.setLabel(
            'bottom', xlabel, color=TEXT_COLOR, **{'font-size': '9pt'}
        )
        plot_item.setLabel(
            'left', ylabel, color=TEXT_COLOR, **{'font-size': '9pt'}
        )
        plot_item.getViewBox().setBackgroundColor(PLOT_BG)
        plot_item.showGrid(x=True, y=True, alpha=0.4)
        for axis_name in ['bottom', 'left']:
            axis = plot_item.getAxis(axis_name)
            axis.setPen(pg.mkPen(color='#aab5c8', width=0.8))
            axis.setTextPen(pg.mkPen(color=TEXT_COLOR))
            axis.setStyle(tickLength=-4)

    def _ensure_family(self, k, nbrfam=1):
        """Create the buffers and curves of family k on first use."""
        if k in self._fam_data:
            return

        # Legends must exist before the curves, or pyqtgraph does not
        # register them.
        if nbrfam > 1 and not self._legends_added:
            self._legends_added = True
            for plot in [
                self.plot_xy,
                self.plot_energy,
                self.plot_rtheta,
                self.plot_zpz,
                self.plot_rdph,
                self.plot_dphw,
            ]:
                legend = plot.addLegend(offset=(5, 5))
                legend.setLabelTextSize('8pt')

        color = FAMILY_COLORS[k % len(FAMILY_COLORS)]
        fam_name = f'Fam {k + 1}'

        sc = QColor(color)
        sc.setAlpha(80)
        scatter_brush = pg.mkBrush(sc)

        self._fam_data[k] = {
            'x': [],
            'y': [],
            'turn': [],
            'energy': [],
            'r': [],
            'theta': [],
            'z': [],
            'pz': [],
            'dph': [],
            'scatter_x': [],
            'scatter_y': [],
            'scatter_z': [],
            'scatter_r': [],
            'scatter_theta': [],
            'scatter_pz': [],
            'scatter_dph': [],
            'scatter_wmevpa': [],
            'scatter_turn': [],
        }

        curves = {}
        # Reference particle
        curves['xy'] = self.plot_xy.plot(
            [], [], pen=pg.mkPen(color=color, width=0.5), name=fam_name
        )
        curves['energy'] = self.plot_energy.plot(
            [],
            [],
            pen=pg.mkPen(color=color, width=1.5),
            clipToView=True,
            name=fam_name,
        )
        curves['energy'].setDownsampling(auto=True, method='peak')
        curves['rtheta'] = self.plot_rtheta.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=2,
            symbolPen=None,
            symbolBrush=color,
            clipToView=True,
            name=fam_name,
        )
        curves['zpz'] = self.plot_zpz.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=2,
            symbolPen=None,
            symbolBrush=color,
            clipToView=True,
            name=fam_name,
        )
        curves['rdph'] = self.plot_rdph.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=2,
            symbolPen=None,
            symbolBrush=color,
            clipToView=True,
            name=fam_name,
        )
        curves['dphw'] = self.plot_dphw.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=2,
            symbolPen=None,
            symbolBrush=color,
            clipToView=True,
            name=fam_name,
        )

        # All particles, semi-transparent
        curves['scat_energy'] = self.plot_energy.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=2,
            symbolPen=None,
            symbolBrush=scatter_brush,
        )
        curves['scat_xy'] = self.plot_xy.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=3,
            symbolPen=None,
            symbolBrush=scatter_brush,
        )
        curves['scat_rtheta'] = self.plot_rtheta.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=3,
            symbolPen=None,
            symbolBrush=scatter_brush,
        )
        curves['scat_zpz'] = self.plot_zpz.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=3,
            symbolPen=None,
            symbolBrush=scatter_brush,
        )
        curves['scat_rdph'] = self.plot_rdph.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=3,
            symbolPen=None,
            symbolBrush=scatter_brush,
        )
        curves['scat_dphw'] = self.plot_dphw.plot(
            [],
            [],
            pen=None,
            symbol='o',
            symbolSize=3,
            symbolPen=None,
            symbolBrush=scatter_brush,
        )

        self._fam_curves[k] = curves

    def append_data(self, data):
        """Buffer one data point without redrawing."""
        k = data.get('family', 0)
        nbrfam = data.get('nbrfam', 1)
        self._ensure_family(k, nbrfam)
        fd = self._fam_data[k]
        fd['x'].append(data['x0_mm'])
        fd['y'].append(data['y0_mm'])
        fd['turn'].append(data['numtou'] - 1 + data['tetnew'] / 360.0)
        fd['energy'].append(data['wmevpa'])
        fd['r'].append(data['r0_mm'])
        fd['theta'].append(data['tetnew'])
        fd['z'].append(data['z0_mm'])
        fd['pz'].append(data['tpz'])
        fd['dph'].append(data.get('dph', 0.0))
        for key in [
            'scatter_x',
            'scatter_y',
            'scatter_z',
            'scatter_r',
            'scatter_theta',
            'scatter_pz',
            'scatter_dph',
            'scatter_wmevpa',
            'scatter_turn',
        ]:
            val = data.get(key, None)
            if val is not None:
                fd[key].append(val)

    def update_plots(self, data):
        """Buffer one data point and redraw."""
        self.append_data(data)
        self._refresh_all_curves()

    def _refresh_all_curves(self):
        """Redraw all curves from the buffers."""
        if not self._fam_data:
            return

        for k, fd in self._fam_data.items():
            if not fd['x']:
                continue
            curves = self._fam_curves[k]

            x_arr = np.array(fd['x'], dtype=np.float64)
            y_arr = np.array(fd['y'], dtype=np.float64)

            if self.rmax_display > 0:
                mask = (np.abs(x_arr) > self.rmax_display) | (
                    np.abs(y_arr) > self.rmax_display
                )
                x_filt = np.where(mask, np.nan, x_arr)
                y_filt = np.where(mask, np.nan, y_arr)
                curves['xy'].setData(x_filt, y_filt, connect='finite')
            else:
                curves['xy'].setData(x_arr, y_arr)

            turn_arr = np.array(fd['turn'], dtype=np.float64)
            energy_arr = np.array(fd['energy'], dtype=np.float64)
            curves['energy'].setData(turn_arr, energy_arr)

            theta_arr = np.array(fd['theta'], dtype=np.float64)
            r_arr = np.array(fd['r'], dtype=np.float64)
            if self.rmax_rtheta > 0:
                mask_rt = r_arr <= self.rmax_rtheta
                curves['rtheta'].setData(theta_arr[mask_rt], r_arr[mask_rt])
            else:
                curves['rtheta'].setData(theta_arr, r_arr)

            z_arr = np.array(fd['z'], dtype=np.float64)
            pz_arr = np.array(fd['pz'], dtype=np.float64)
            curves['zpz'].setData(z_arr, pz_arr)

            dph_arr = np.array(fd['dph'], dtype=np.float64)
            mask_rdph = np.ones(len(r_arr), dtype=bool)
            if self.rmax_rdph > 0:
                mask_rdph &= r_arr <= self.rmax_rdph
            if self.dphi_min != 0.0 or self.dphi_max != 0.0:
                mask_rdph &= (dph_arr >= self.dphi_min) & (
                    dph_arr <= self.dphi_max
                )
            curves['rdph'].setData(r_arr[mask_rdph], dph_arr[mask_rdph])

            mask_dphw = np.ones(len(dph_arr), dtype=bool)
            if self.dphi_min != 0.0 or self.dphi_max != 0.0:
                mask_dphw = (dph_arr >= self.dphi_min) & (
                    dph_arr <= self.dphi_max
                )
            curves['dphw'].setData(turn_arr[mask_dphw], dph_arr[mask_dphw])

            sx_list = fd.get('scatter_x', [])
            if sx_list:
                sx = np.concatenate(sx_list)
                sy = np.concatenate(fd['scatter_y'])
                sz = np.concatenate(fd['scatter_z'])
                sr = np.concatenate(fd['scatter_r'])
                stheta = np.concatenate(fd['scatter_theta'])
                spz = np.concatenate(fd['scatter_pz'])
                sdph = np.concatenate(fd['scatter_dph'])
                swmevpa = np.concatenate(fd['scatter_wmevpa'])
                sturn_list = fd.get('scatter_turn', [])
                sturn = (
                    np.concatenate(sturn_list) if sturn_list else np.array([])
                )

                # Decimate to keep redraws fast.
                MAX_SCAT = 200000
                if len(sx) > MAX_SCAT:
                    idx = np.linspace(0, len(sx) - 1, MAX_SCAT, dtype=int)
                    sx, sy, sz = sx[idx], sy[idx], sz[idx]
                    sr, stheta, spz, sdph = (
                        sr[idx],
                        stheta[idx],
                        spz[idx],
                        sdph[idx],
                    )
                    swmevpa = swmevpa[idx]
                    if len(sturn) > 0:
                        sturn = sturn[idx]

                if self.rmax_display > 0:
                    smask = (np.abs(sx) <= self.rmax_display) & (
                        np.abs(sy) <= self.rmax_display
                    )
                    curves['scat_xy'].setData(sx[smask], sy[smask])
                else:
                    curves['scat_xy'].setData(sx, sy)

                if self.rmax_rtheta > 0:
                    smask_rt = sr <= self.rmax_rtheta
                    curves['scat_rtheta'].setData(
                        stheta[smask_rt], sr[smask_rt]
                    )
                else:
                    curves['scat_rtheta'].setData(stheta, sr)

                curves['scat_zpz'].setData(sz, spz)

                smask_rdph = np.ones(len(sr), dtype=bool)
                if self.rmax_rdph > 0:
                    smask_rdph &= sr <= self.rmax_rdph
                if self.dphi_min != 0.0 or self.dphi_max != 0.0:
                    smask_rdph &= (sdph >= self.dphi_min) & (
                        sdph <= self.dphi_max
                    )
                curves['scat_rdph'].setData(sr[smask_rdph], sdph[smask_rdph])

                smask_dphw = np.ones(len(sdph), dtype=bool)
                if self.dphi_min != 0.0 or self.dphi_max != 0.0:
                    smask_dphw = (sdph >= self.dphi_min) & (
                        sdph <= self.dphi_max
                    )
                if len(sturn) > 0:
                    curves['scat_dphw'].setData(
                        sturn[smask_dphw], sdph[smask_dphw]
                    )

                if len(sturn) > 0:
                    curves['scat_energy'].setData(sturn, swmevpa)

        if self.rmax_display > 0:
            self.plot_xy.setXRange(
                -self.rmax_display, self.rmax_display, padding=0
            )
            self.plot_xy.setYRange(
                -self.rmax_display, self.rmax_display, padding=0
            )
        if self.rmax_rtheta > 0:
            self.plot_rtheta.setYRange(0, self.rmax_rtheta, padding=0.02)
        if self.rmax_rdph > 0:
            self.plot_rdph.setXRange(0, self.rmax_rdph, padding=0.02)

    def clear_plots(self):
        """Remove all curves, deflectors and legends."""
        for item in self._deflector_items:
            self.plot_xy.removeItem(item)
        self._deflector_items.clear()

        for k, curves in self._fam_curves.items():
            self.plot_xy.removeItem(curves['xy'])
            self.plot_xy.removeItem(curves['scat_xy'])
            self.plot_energy.removeItem(curves['energy'])
            self.plot_energy.removeItem(curves['scat_energy'])
            self.plot_rtheta.removeItem(curves['rtheta'])
            self.plot_rtheta.removeItem(curves['scat_rtheta'])
            self.plot_zpz.removeItem(curves['zpz'])
            self.plot_zpz.removeItem(curves['scat_zpz'])
            self.plot_rdph.removeItem(curves['rdph'])
            self.plot_rdph.removeItem(curves['scat_rdph'])
            self.plot_dphw.removeItem(curves['dphw'])
            self.plot_dphw.removeItem(curves['scat_dphw'])
        self._fam_data.clear()
        self._fam_curves.clear()

        for plot in [
            self.plot_xy,
            self.plot_energy,
            self.plot_rtheta,
            self.plot_zpz,
            self.plot_rdph,
            self.plot_dphw,
        ]:
            if getattr(plot, 'legend', None) is not None:
                try:
                    plot.legend.clear()
                    scene = plot.legend.scene()
                    if scene is not None:
                        scene.removeItem(plot.legend)
                except Exception:
                    pass
                plot.legend = None
            plot.enableAutoRange()
        self._legends_added = False

    def reset_zoom(self):
        """Reset the zoom on all plots, keeping the radius limits."""
        for plot in [
            self.plot_xy,
            self.plot_energy,
            self.plot_rtheta,
            self.plot_zpz,
            self.plot_rdph,
            self.plot_dphw,
        ]:
            plot.enableAutoRange()
        if self.rmax_display > 0:
            self.plot_xy.setXRange(
                -self.rmax_display, self.rmax_display, padding=0
            )
            self.plot_xy.setYRange(
                -self.rmax_display, self.rmax_display, padding=0
            )
        if self.rmax_rtheta > 0:
            self.plot_rtheta.setYRange(0, self.rmax_rtheta, padding=0.02)
        if self.rmax_rdph > 0:
            self.plot_rdph.setXRange(0, self.rmax_rdph, padding=0.02)

    def draw_deflectors(self, defl_list):
        """Draw the deflectors on the XY plot.

        defl_list: list of dicts with keys
            rdf_m, xdfi_m, ydfi_m, xdfo_m, ydfo_m, tetdfi_deg,
            tetdfo_deg
        """
        for item in self._deflector_items:
            self.plot_xy.removeItem(item)
        self._deflector_items.clear()

        if not defl_list:
            return

        defl_color = '#e53935'
        defl_pen = pg.mkPen(color=defl_color, width=2.5)
        defl_pen_dash = pg.mkPen(
            color=defl_color, width=1.5, style=Qt.DashLine
        )

        for idx, defl in enumerate(defl_list):
            rdf_mm = defl['rdf_m'] * 1000.0
            xdfi_mm = defl['xdfi_m'] * 1000.0
            ydfi_mm = defl['ydfi_m'] * 1000.0
            xdfo_mm = defl['xdfo_m'] * 1000.0
            ydfo_mm = defl['ydfo_m'] * 1000.0
            tetdfi_deg = defl['tetdfi_deg']
            tetdfo_deg = defl['tetdfo_deg']

            # Unwrap so that the exit angle follows the entry angle.
            a1 = tetdfi_deg % 360
            a2 = tetdfo_deg % 360
            if a2 <= a1:
                a2 += 360
            a1_rad = np.radians(a1)
            a2_rad = np.radians(a2)

            ext = rdf_mm * 0.15  # radial depth drawn: 15% of rdf

            # Inner arc at rdf
            angles = np.linspace(a1_rad, a2_rad, 100)
            arc_x = rdf_mm * np.cos(angles)
            arc_y = rdf_mm * np.sin(angles)
            item = self.plot_xy.plot(arc_x, arc_y, pen=defl_pen)
            self._deflector_items.append(item)

            # Entry wall
            cos_a1 = np.cos(a1_rad)
            sin_a1 = np.sin(a1_rad)
            # t_start: where the wall line crosses the rdf circle.
            b_in = xdfi_mm * cos_a1 + ydfi_mm * sin_a1
            c_in = xdfi_mm**2 + ydfi_mm**2 - rdf_mm**2
            disc_in = b_in**2 - c_in
            if disc_in >= 0:
                t_start_in = -b_in + np.sqrt(disc_in)
            else:
                t_start_in = rdf_mm  # fallback
            t_end_in = t_start_in + ext
            x_in = [xdfi_mm + t_start_in * cos_a1, xdfi_mm + t_end_in * cos_a1]
            y_in = [ydfi_mm + t_start_in * sin_a1, ydfi_mm + t_end_in * sin_a1]
            item = self.plot_xy.plot(x_in, y_in, pen=defl_pen)
            self._deflector_items.append(item)

            # Exit wall
            cos_a2 = np.cos(a2_rad)
            sin_a2 = np.sin(a2_rad)
            b_out = xdfo_mm * cos_a2 + ydfo_mm * sin_a2
            c_out = xdfo_mm**2 + ydfo_mm**2 - rdf_mm**2
            disc_out = b_out**2 - c_out
            if disc_out >= 0:
                t_start_out = -b_out + np.sqrt(disc_out)
            else:
                t_start_out = rdf_mm
            t_end_out = t_start_out + ext
            x_out = [
                xdfo_mm + t_start_out * cos_a2,
                xdfo_mm + t_end_out * cos_a2,
            ]
            y_out = [
                ydfo_mm + t_start_out * sin_a2,
                ydfo_mm + t_end_out * sin_a2,
            ]
            item = self.plot_xy.plot(x_out, y_out, pen=defl_pen)
            self._deflector_items.append(item)

            # Outer arc, closes the outline
            outer_x = (rdf_mm + ext) * np.cos(angles)
            outer_y = (rdf_mm + ext) * np.sin(angles)
            item = self.plot_xy.plot(outer_x, outer_y, pen=defl_pen_dash)
            self._deflector_items.append(item)

            mid_angle = (a1_rad + a2_rad) / 2.0
            label_r = rdf_mm + ext + 40
            label = pg.TextItem(
                f'Defl. {idx + 1}', color=defl_color, anchor=(0.5, 0.5)
            )
            label.setFont(QFont('sans-serif', 9, QFont.Bold))
            label.setPos(
                label_r * np.cos(mid_angle), label_r * np.sin(mid_angle)
            )
            self.plot_xy.addItem(label)
            self._deflector_items.append(label)

    def export_image(self, path, dpi=150):
        """Export the plots to SVG, PDF or a raster image."""
        if path.lower().endswith('.svg'):
            from pyqtgraph.exporters import SVGExporter

            exporter = SVGExporter(self.graphics_widget.scene())
            exporter.export(path)
        elif path.lower().endswith('.pdf'):
            from PyQt5.QtPrintSupport import QPrinter
            from PyQt5.QtGui import QPainter

            printer = QPrinter(QPrinter.HighResolution)
            printer.setOutputFormat(QPrinter.PdfFormat)
            printer.setOutputFileName(path)
            painter = QPainter(printer)
            self.graphics_widget.render(painter)
            painter.end()
        else:
            width = int(10 * dpi)
            height = int(8 * dpi)
            pixmap = self.graphics_widget.grab()
            scaled = pixmap.scaled(
                width, height, Qt.KeepAspectRatio, Qt.SmoothTransformation
            )
            scaled.save(path)


class _ToolbarShim:
    """Provides toolbar.home(), the matplotlib API app.py uses."""

    def __init__(self, canvas):
        self._canvas = canvas

    def home(self):
        self._canvas.reset_zoom()


class _FigShim:
    """Provides fig.savefig(), the matplotlib API app.py uses."""

    def __init__(self, canvas):
        self._canvas = canvas

    def savefig(self, path, dpi=150, bbox_inches='tight'):
        self._canvas.export_image(path, dpi=dpi)


class PlotWidget(QWidget):
    """Plot canvas with its control bar."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("background-color: #f0f2f5;")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)

        self.canvas = PlotCanvas(self)
        self.toolbar = _ToolbarShim(self.canvas)

        ctrl_layout = QHBoxLayout()
        ctrl_layout.setContentsMargins(4, 0, 4, 0)

        btn_style = (
            "QPushButton { color: #1a1a2e; background: #f0f4fa; "
            "  border: 1px solid #aab5c8; border-radius: 3px; "
            "  padding: 3px 8px; font-weight: bold; font-size: 10px; }"
            "QPushButton:hover { background: #c0d0f0; }"
            "QPushButton:pressed { background: #a0b8e0; }"
        )

        btn_home = QPushButton("\u2922 Unzoom")
        btn_home.setToolTip("Reset the zoom on all plots")
        btn_home.setStyleSheet(btn_style)
        btn_home.clicked.connect(self.canvas.reset_zoom)
        ctrl_layout.addWidget(btn_home)

        ctrl_layout.addSpacing(12)

        lbl = QLabel("R max (mm):")
        lbl.setStyleSheet(
            "color: #2c3e50; font-weight: bold; font-size: 10px;"
        )
        ctrl_layout.addWidget(lbl)

        self.rmax_spin = QDoubleSpinBox()
        self.rmax_spin.setRange(0, 1e8)
        self.rmax_spin.setDecimals(0)
        self.rmax_spin.setValue(3000)
        self.rmax_spin.setSpecialValueText("Auto")
        self.rmax_spin.setSingleStep(500)
        self.rmax_spin.setToolTip("0 = auto. Increase to see the extraction.")
        self.rmax_spin.setStyleSheet(
            "QDoubleSpinBox { background: #ffffff; color: #1565C0; "
            "  border: 1px solid #aab5c8; border-radius: 3px; padding: 2px 6px; "
            "  font-weight: bold; }"
            "QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { "
            "  background: #e0e8f0; border: 1px solid #aab5c8; }"
        )
        self.rmax_spin.valueChanged.connect(self._on_rmax_changed)
        ctrl_layout.addWidget(self.rmax_spin)

        self.auto_check = QCheckBox("Auto")
        self.auto_check.setChecked(False)
        self.auto_check.setStyleSheet(
            "QCheckBox { color: #2c3e50; font-size: 10px; }"
            "QCheckBox::indicator { width: 14px; height: 14px; }"
        )
        self.auto_check.toggled.connect(self._on_auto_toggled)
        ctrl_layout.addWidget(self.auto_check)

        ctrl_layout.addSpacing(20)

        lbl_rt = QLabel("R max R-\u03b8 (mm):")
        lbl_rt.setStyleSheet(
            "color: #2c3e50; font-weight: bold; font-size: 10px;"
        )
        ctrl_layout.addWidget(lbl_rt)

        self.rmax_rtheta_spin = QDoubleSpinBox()
        self.rmax_rtheta_spin.setRange(0, 1e8)
        self.rmax_rtheta_spin.setDecimals(0)
        self.rmax_rtheta_spin.setValue(3000)
        self.rmax_rtheta_spin.setSpecialValueText("Auto")
        self.rmax_rtheta_spin.setSingleStep(500)
        self.rmax_rtheta_spin.setToolTip(
            "0 = auto (no filter). Hides points with R above this value on the R-\u03b8 plot."
        )
        self.rmax_rtheta_spin.setStyleSheet(
            "QDoubleSpinBox { background: #ffffff; color: #2e7d32; "
            "  border: 1px solid #aab5c8; border-radius: 3px; padding: 2px 6px; "
            "  font-weight: bold; }"
            "QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { "
            "  background: #e0e8f0; border: 1px solid #aab5c8; }"
        )
        self.rmax_rtheta_spin.valueChanged.connect(
            self._on_rmax_rtheta_changed
        )
        ctrl_layout.addWidget(self.rmax_rtheta_spin)

        self.auto_rtheta_check = QCheckBox("Auto")
        self.auto_rtheta_check.setChecked(False)
        self.auto_rtheta_check.setStyleSheet(
            "QCheckBox { color: #2c3e50; font-size: 10px; }"
            "QCheckBox::indicator { width: 14px; height: 14px; }"
        )
        self.auto_rtheta_check.toggled.connect(self._on_auto_rtheta_toggled)
        ctrl_layout.addWidget(self.auto_rtheta_check)

        ctrl_layout.addSpacing(20)

        lbl_rdph = QLabel("R max R-d\u03c6 (mm):")
        lbl_rdph.setStyleSheet(
            "color: #2c3e50; font-weight: bold; font-size: 10px;"
        )
        ctrl_layout.addWidget(lbl_rdph)

        self.rmax_rdph_spin = QDoubleSpinBox()
        self.rmax_rdph_spin.setRange(0, 1e8)
        self.rmax_rdph_spin.setDecimals(0)
        self.rmax_rdph_spin.setValue(3000)
        self.rmax_rdph_spin.setSpecialValueText("Auto")
        self.rmax_rdph_spin.setSingleStep(500)
        self.rmax_rdph_spin.setToolTip(
            "0 = auto (no filter). Hides points with R above this value on the R-d\u03c6 plot."
        )
        self.rmax_rdph_spin.setStyleSheet(
            "QDoubleSpinBox { background: #ffffff; color: #2e7d32; "
            "  border: 1px solid #aab5c8; border-radius: 3px; padding: 2px 6px; "
            "  font-weight: bold; }"
            "QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { "
            "  background: #e0e8f0; border: 1px solid #aab5c8; }"
        )
        self.rmax_rdph_spin.valueChanged.connect(self._on_rmax_rdph_changed)
        ctrl_layout.addWidget(self.rmax_rdph_spin)

        self.auto_rdph_check = QCheckBox("Auto")
        self.auto_rdph_check.setChecked(False)
        self.auto_rdph_check.setStyleSheet(
            "QCheckBox { color: #2c3e50; font-size: 10px; }"
            "QCheckBox::indicator { width: 14px; height: 14px; }"
        )
        self.auto_rdph_check.toggled.connect(self._on_auto_rdph_toggled)
        ctrl_layout.addWidget(self.auto_rdph_check)

        ctrl_layout.addSpacing(20)

        dph_spin_style = (
            "QDoubleSpinBox { background: #ffffff; color: #7b1fa2; "
            "  border: 1px solid #aab5c8; border-radius: 3px; padding: 2px 6px; "
            "  font-weight: bold; }"
            "QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { "
            "  background: #e0e8f0; border: 1px solid #aab5c8; }"
        )

        lbl_dph = QLabel("d\u03c6:")
        lbl_dph.setStyleSheet(
            "color: #2c3e50; font-weight: bold; font-size: 10px;"
        )
        ctrl_layout.addWidget(lbl_dph)

        lbl_dph_min = QLabel("min")
        lbl_dph_min.setStyleSheet("color: #7b1fa2; font-size: 10px;")
        ctrl_layout.addWidget(lbl_dph_min)

        self.dphi_min_spin = QDoubleSpinBox()
        self.dphi_min_spin.setRange(-1e6, 1e6)
        self.dphi_min_spin.setDecimals(1)
        self.dphi_min_spin.setValue(0)
        self.dphi_min_spin.setSingleStep(10)
        self.dphi_min_spin.setSuffix("\u00b0")
        self.dphi_min_spin.setToolTip(
            "d\u03c6 minimum (deg). No filter when min = max = 0."
        )
        self.dphi_min_spin.setStyleSheet(dph_spin_style)
        self.dphi_min_spin.valueChanged.connect(self._on_dphi_filter_changed)
        ctrl_layout.addWidget(self.dphi_min_spin)

        lbl_dph_max = QLabel("max")
        lbl_dph_max.setStyleSheet("color: #7b1fa2; font-size: 10px;")
        ctrl_layout.addWidget(lbl_dph_max)

        self.dphi_max_spin = QDoubleSpinBox()
        self.dphi_max_spin.setRange(-1e6, 1e6)
        self.dphi_max_spin.setDecimals(1)
        self.dphi_max_spin.setValue(0)
        self.dphi_max_spin.setSingleStep(10)
        self.dphi_max_spin.setSuffix("\u00b0")
        self.dphi_max_spin.setToolTip(
            "d\u03c6 maximum (deg). No filter when min = max = 0."
        )
        self.dphi_max_spin.setStyleSheet(dph_spin_style)
        self.dphi_max_spin.valueChanged.connect(self._on_dphi_filter_changed)
        ctrl_layout.addWidget(self.dphi_max_spin)

        ctrl_layout.addStretch()

        layout.addLayout(ctrl_layout)
        layout.addWidget(self.canvas)

        self.canvas.rmax_display = 3000.0
        self.canvas.rmax_rtheta = 3000.0
        self.canvas.rmax_rdph = 3000.0

    def _on_rmax_changed(self, val):
        self.canvas.rmax_display = val
        self.canvas._refresh_all_curves()

    def _on_auto_toggled(self, checked):
        if checked:
            self.rmax_spin.setValue(0)
        else:
            self.rmax_spin.setValue(3000)

    def _on_rmax_rtheta_changed(self, val):
        self.canvas.rmax_rtheta = val
        self.canvas._refresh_all_curves()

    def _on_auto_rtheta_toggled(self, checked):
        if checked:
            self.rmax_rtheta_spin.setValue(0)
        else:
            self.rmax_rtheta_spin.setValue(3000)

    def _on_rmax_rdph_changed(self, val):
        self.canvas.rmax_rdph = val
        self.canvas._refresh_all_curves()

    def _on_auto_rdph_toggled(self, checked):
        if checked:
            self.rmax_rdph_spin.setValue(0)
        else:
            self.rmax_rdph_spin.setValue(3000)

    def _on_dphi_filter_changed(self, _=None):
        self.canvas.dphi_min = self.dphi_min_spin.value()
        self.canvas.dphi_max = self.dphi_max_spin.value()
        self.canvas._refresh_all_curves()


class Plot3DWindow(QDialog):
    """Separate window with the 3D trajectory view (Plotly WebGL)."""

    SCAT_COLORS = [
        '#ff0000',
        '#00cc00',
        '#ff8000',
        '#cc00cc',
        '#00cccc',
        '#ffcc00',
        '#0066ff',
        '#ff3399',
        '#66ff33',
        '#9933ff',
    ]
    TRAJ_COLORS = [
        '#1565C0',
        '#d32f2f',
        '#2e7d32',
        '#ff6f00',
        '#7b1fa2',
        '#00838f',
        '#c62828',
        '#283593',
        '#558b2f',
        '#4e342e',
    ]

    def __init__(self, canvas, parent=None):
        super().__init__(parent)
        self.setWindowTitle("LIONS - 3D View")
        self.resize(1000, 750)
        self.canvas = canvas
        self.z_scale = 200.0
        self._defl_list = []
        self._available = HAS_PLOTLY and HAS_WEBENGINE

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        if not self._available:
            msg = []
            if not HAS_PLOTLY:
                msg.append("pip install plotly")
            if not HAS_WEBENGINE:
                msg.append("pip install PyQtWebEngine")
            lbl = QLabel("The 3D view requires:\n" + "\n".join(msg))
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet(
                "color: #c62828; font-size: 14px; font-weight: bold;"
            )
            layout.addWidget(lbl)
            return

        ctrl = QHBoxLayout()
        ctrl.setContentsMargins(4, 0, 4, 0)

        _spin_style = (
            "QSpinBox, QDoubleSpinBox { background: #ffffff; color: #1565C0; "
            "  border: 1px solid #aab5c8; border-radius: 3px; "
            "  padding: 2px 6px; font-weight: bold; font-size: 11px; }"
            "QSpinBox::up-button, QSpinBox::down-button, "
            "QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { "
            "  background: #e0e8f0; border: 1px solid #aab5c8; }"
        )
        _btn_style = (
            "QPushButton { color: #1a1a2e; background: #f0f4fa; "
            "  border: 1px solid #aab5c8; border-radius: 3px; "
            "  padding: 3px 10px; font-weight: bold; font-size: 10px; }"
            "QPushButton:hover { background: #c0d0f0; }"
        )

        lbl_z = QLabel("Z scale:")
        lbl_z.setStyleSheet(
            "color: #2c3e50; font-weight: bold; font-size: 11px;"
        )
        ctrl.addWidget(lbl_z)

        self.z_scale_spin = QSpinBox()
        self.z_scale_spin.setRange(1, 10000)
        self.z_scale_spin.setValue(200)
        self.z_scale_spin.setSingleStep(50)
        self.z_scale_spin.setSuffix("x")
        self.z_scale_spin.setToolTip("Z magnification factor")
        self.z_scale_spin.setStyleSheet(_spin_style)
        self.z_scale_spin.valueChanged.connect(self._on_z_scale_changed)
        ctrl.addWidget(self.z_scale_spin)

        ctrl.addSpacing(12)

        lbl_rmin = QLabel("Rmin:")
        lbl_rmin.setStyleSheet(
            "color: #2c3e50; font-weight: bold; font-size: 10px;"
        )
        ctrl.addWidget(lbl_rmin)

        self.rmin_3d_spin = QDoubleSpinBox()
        self.rmin_3d_spin.setRange(0, 9999)
        self.rmin_3d_spin.setValue(0)
        self.rmin_3d_spin.setSingleStep(100)
        self.rmin_3d_spin.setSuffix(" mm")
        self.rmin_3d_spin.setToolTip("Minimum radius (0 = no limit)")
        self.rmin_3d_spin.setStyleSheet(_spin_style)
        self.rmin_3d_spin.valueChanged.connect(self._on_range_changed)
        ctrl.addWidget(self.rmin_3d_spin)

        lbl_rmax = QLabel("Rmax:")
        lbl_rmax.setStyleSheet(
            "color: #2c3e50; font-weight: bold; font-size: 10px;"
        )
        ctrl.addWidget(lbl_rmax)

        self.rmax_3d_spin = QDoubleSpinBox()
        self.rmax_3d_spin.setRange(0, 9999)
        self.rmax_3d_spin.setValue(0)
        self.rmax_3d_spin.setSingleStep(100)
        self.rmax_3d_spin.setSuffix(" mm")
        self.rmax_3d_spin.setToolTip("Maximum radius (0 = no limit)")
        self.rmax_3d_spin.setStyleSheet(_spin_style)
        self.rmax_3d_spin.valueChanged.connect(self._on_range_changed)
        ctrl.addWidget(self.rmax_3d_spin)

        ctrl.addSpacing(12)

        btn_refresh = QPushButton("\u21bb Refresh")
        btn_refresh.setToolTip("Rebuild the 3D view")
        btn_refresh.setStyleSheet(_btn_style)
        btn_refresh.clicked.connect(self.refresh_3d)
        ctrl.addWidget(btn_refresh)

        btn_export = QPushButton("\u2913 Export HTML")
        btn_export.setToolTip("Save the 3D view as an interactive HTML file")
        btn_export.setStyleSheet(_btn_style)
        btn_export.clicked.connect(self._on_export)
        ctrl.addWidget(btn_export)

        ctrl.addStretch()
        layout.addLayout(ctrl)

        self._web = QWebEngineView()
        layout.addWidget(self._web)

    def _on_z_scale_changed(self, val):
        self.z_scale = float(val)
        self.refresh_3d()

    def _on_range_changed(self, _=None):
        self.refresh_3d()

    def _on_export(self):
        if not self._available:
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export 3D view",
            "lions_3d.html",
            "Interactive HTML (*.html);;PNG (*.png)",
        )
        if not path:
            return
        fig = self._build_figure()
        if path.lower().endswith('.png'):
            try:
                fig.write_image(path, width=1600, height=1200, scale=2)
            except Exception:
                # PNG export needs kaleido; fall back to HTML.
                path = path.rsplit('.', 1)[0] + '.html'
                fig.write_html(path, include_plotlyjs='cdn')
        else:
            fig.write_html(path, include_plotlyjs='cdn')
        self.setWindowTitle(f"LIONS - 3D View  [Exported: {path}]")

    def _filter_by_radius(self, x, y, *arrays):
        """Keep the points with Rmin <= r <= Rmax; return x, y and the
        arrays, filtered.
        """
        rmin = self.rmin_3d_spin.value()
        rmax = self.rmax_3d_spin.value()
        if rmin <= 0 and rmax <= 0:
            return (x, y) + arrays
        r2 = x**2 + y**2
        mask = np.ones(len(x), dtype=bool)
        if rmin > 0:
            mask &= r2 >= rmin**2
        if rmax > 0:
            mask &= r2 <= rmax**2
        return tuple(a[mask] for a in (x, y) + arrays)

    def _build_figure(self):
        """Build the Plotly 3D figure."""
        go = _plotly_go
        fig = go.Figure()

        fam_data = self.canvas._fam_data

        for k, fd in fam_data.items():
            if not fd['x'] or len(fd['x']) < 2:
                continue

            traj_color = self.TRAJ_COLORS[k % len(self.TRAJ_COLORS)]
            scat_color = self.SCAT_COLORS[k % len(self.SCAT_COLORS)]
            fam_name = f'Fam {k + 1}'

            # Reference particle
            x_arr = np.array(fd['x'], dtype=np.float64)
            y_arr = np.array(fd['y'], dtype=np.float64)
            z_arr = np.array(fd['z'], dtype=np.float64) * self.z_scale

            x_arr, y_arr, z_arr = self._filter_by_radius(x_arr, y_arr, z_arr)

            fig.add_trace(
                go.Scatter3d(
                    x=x_arr,
                    y=y_arr,
                    z=z_arr,
                    mode='lines',
                    line=dict(color=traj_color, width=3),
                    name=f'{fam_name} trajectory',
                    legendgroup=fam_name,
                )
            )

            # All particles
            sx_list = fd.get('scatter_x', [])
            if sx_list:
                sx = np.concatenate(sx_list)
                sy = np.concatenate(fd['scatter_y'])
                sz = np.concatenate(fd['scatter_z']) * self.z_scale

                # WebGL stays responsive up to about 100k points.
                MAX_SCAT = 80000
                if len(sx) > MAX_SCAT:
                    idx = np.random.choice(len(sx), MAX_SCAT, replace=False)
                    idx.sort()
                    sx, sy, sz = sx[idx], sy[idx], sz[idx]

                sx, sy, sz = self._filter_by_radius(sx, sy, sz)

                fig.add_trace(
                    go.Scatter3d(
                        x=sx,
                        y=sy,
                        z=sz,
                        mode='markers',
                        marker=dict(size=2, color=scat_color, opacity=0.6),
                        name=f'{fam_name} bunch',
                        legendgroup=fam_name,
                    )
                )

        # Deflectors
        for i, defl in enumerate(self._defl_list):
            rdf_mm = defl['rdf_m'] * 1000.0
            tetdfi_deg = defl['tetdfi_deg']
            tetdfo_deg = defl['tetdfo_deg']

            a1 = tetdfi_deg % 360
            a2 = tetdfo_deg % 360
            if a2 <= a1:
                a2 += 360

            ext = rdf_mm * 0.15
            angles = np.linspace(np.radians(a1), np.radians(a2), 80)

            arc_x = rdf_mm * np.cos(angles)
            arc_y = rdf_mm * np.sin(angles)
            # Outer arc, reversed to close the outline
            out_x = (rdf_mm + ext) * np.cos(angles[::-1])
            out_y = (rdf_mm + ext) * np.sin(angles[::-1])
            cx = np.concatenate([arc_x, out_x, arc_x[:1]])
            cy = np.concatenate([arc_y, out_y, arc_y[:1]])
            cz = np.zeros_like(cx)

            fig.add_trace(
                go.Scatter3d(
                    x=cx,
                    y=cy,
                    z=cz,
                    mode='lines',
                    line=dict(color='#e53935', width=4),
                    name=f'Defl. {i + 1}',
                    showlegend=True,
                )
            )

        # Z is drawn scaled; label the ticks in real mm.
        all_z_scaled = []
        for tr in fig.data:
            if tr.z is not None and len(tr.z) > 0:
                zv = np.asarray(tr.z, dtype=np.float64)
                zv = zv[np.isfinite(zv)]
                if len(zv) > 0:
                    all_z_scaled.append(zv)
        if all_z_scaled:
            z_all = np.concatenate(all_z_scaled)
            z_real_min = z_all.min() / self.z_scale
            z_real_max = z_all.max() / self.z_scale
            z_range = z_real_max - z_real_min
            if z_range > 0:
                step = z_range / 8.0
                mag = 10 ** np.floor(np.log10(step))
                nice = [1, 2, 2.5, 5, 10]
                step = min(nice, key=lambda x: abs(x * mag - step)) * mag
            else:
                step = 1.0
            t_start = np.floor(z_real_min / step) * step
            t_end = np.ceil(z_real_max / step) * step
            tick_real = np.arange(t_start, t_end + step * 0.5, step)
            tick_scaled = tick_real * self.z_scale
            tick_text = [f'{v:.3g}' for v in tick_real]
            zaxis_dict = dict(
                title='Z (mm)',
                backgroundcolor='#f0f0f5',
                gridcolor='#dee2e6',
                showbackground=True,
                tickvals=tick_scaled.tolist(),
                ticktext=tick_text,
            )
        else:
            zaxis_dict = dict(
                title='Z (mm)',
                backgroundcolor='#f0f0f5',
                gridcolor='#dee2e6',
                showbackground=True,
            )

        fig.update_layout(
            title=dict(
                text='LIONS - 3D trajectory  (Z = f(X, Y))',
                font=dict(size=16, color='#1a1a2e'),
            ),
            scene=dict(
                xaxis=dict(
                    title='X (mm)',
                    backgroundcolor='#f8f9fa',
                    gridcolor='#dee2e6',
                    showbackground=True,
                ),
                yaxis=dict(
                    title='Y (mm)',
                    backgroundcolor='#f8f9fa',
                    gridcolor='#dee2e6',
                    showbackground=True,
                ),
                zaxis=zaxis_dict,
                aspectmode='data',
                camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
            ),
            legend=dict(
                bgcolor='rgba(255,255,255,0.85)',
                bordercolor='#aab5c8',
                borderwidth=1,
                font=dict(size=11),
            ),
            margin=dict(l=0, r=0, t=40, b=0),
            paper_bgcolor='#ffffff',
        )

        return fig

    def refresh_3d(self):
        """Rebuild and show the Plotly 3D figure."""
        if not self._available:
            return
        fig = self._build_figure()
        # Through a file: setHtml is limited to 2 MB, smaller than
        # inlined plotly.js.
        tmp_dir = tempfile.gettempdir()
        tmp_path = os.path.join(tmp_dir, '_lions_3d_view.html')
        fig.write_html(
            tmp_path,
            include_plotlyjs=True,
            config={
                'displayModeBar': True,
                'scrollZoom': True,
            },
        )
        self._web.load(QUrl.fromLocalFile(tmp_path))

    def clear_3d(self):
        """Clear the 3D view."""
        if not self._available:
            return
        self._defl_list = []
        self._web.setHtml(
            '<html><body style="background:#f0f2f5;"></body></html>'
        )

    def draw_deflectors_3d(self, defl_list):
        """Store the deflectors for the next refresh."""
        self._defl_list = defl_list or []
