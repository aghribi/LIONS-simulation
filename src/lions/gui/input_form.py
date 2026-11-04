"""Input form for the LIONS simulation: every data file parameter plus
the runtime answers.
"""

import json
from dataclasses import astuple
import math
import os
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QLabel,
    QLineEdit,
    QDoubleSpinBox,
    QSpinBox,
    QComboBox,
    QPushButton,
    QFileDialog,
    QScrollArea,
    QToolButton,
    QFrame,
    QGroupBox,
    QSizePolicy,
    QPlainTextEdit,
    QMessageBox,
)
from PyQt5.QtCore import Qt, pyqtSignal

from ..constants import AMU, C_LIGHT, ELECTRON_MASS_U, MEV
from ..deck import (
    OPTIONS,
    Deck,
    Deflector,
    Dipole,
    Family,
    FieldHarmonic,
    Grid,
    read_deck,
    write_deck,
)
from ..fortran_text import parse_floats


class CollapsibleSection(QWidget):
    """Collapsible section with a clickable title."""

    def __init__(self, title, parent=None, collapsed=False):
        super().__init__(parent)
        self._is_collapsed = collapsed

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.toggle_btn = QToolButton()
        self.toggle_btn.setStyleSheet(
            "QToolButton { font-weight: bold; border: none; "
            "background: #e0e0e0; padding: 6px; }"
            "QToolButton:hover { background: #d0d0d0; }"
        )
        self.toggle_btn.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle_btn.setArrowType(
            Qt.RightArrow if collapsed else Qt.DownArrow
        )
        self.toggle_btn.setText(title)
        self.toggle_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.toggle_btn.clicked.connect(self._toggle)

        self.content = QWidget()
        self.content.setVisible(not collapsed)
        self.content_layout = QFormLayout(self.content)
        self.content_layout.setContentsMargins(10, 5, 5, 5)

        layout.addWidget(self.toggle_btn)
        layout.addWidget(self.content)

    def _toggle(self):
        self._is_collapsed = not self._is_collapsed
        self.content.setVisible(not self._is_collapsed)
        self.toggle_btn.setArrowType(
            Qt.RightArrow if self._is_collapsed else Qt.DownArrow
        )

    def add_row(self, label, widget):
        """Add a labelled widget to the section, return the widget."""
        self.content_layout.addRow(label, widget)
        return widget


def _make_double_spin(value=0.0, minimum=-1e12, maximum=1e12, decimals=6):
    spin = QDoubleSpinBox()
    spin.setRange(minimum, maximum)
    spin.setDecimals(decimals)
    spin.setValue(value)
    spin.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    return spin


def _make_int_spin(value=0, minimum=-999999, maximum=999999):
    spin = QSpinBox()
    spin.setRange(minimum, maximum)
    spin.setValue(value)
    spin.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    return spin


def _make_file_row(default='', browse_cb=None):
    """QLineEdit with a browse button."""
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    edit = QLineEdit(default)
    btn = QPushButton('...')
    btn.setFixedWidth(30)
    if browse_cb:
        btn.clicked.connect(browse_cb)
    lay.addWidget(edit)
    lay.addWidget(btn)
    return w, edit, btn


def _numbers(values):
    """Values as text that reads back to the same floats."""
    return ' '.join(repr(float(v)) for v in values)


def _rows(text, width, count):
    """count rows of width floats from a text box, padded with zeros."""
    rows = [
        (parse_floats(line) + [0.0] * width)[:width]
        for line in text.splitlines()
        if line.strip()
    ]
    rows += [[0.0] * width for _ in range(count - len(rows))]
    return rows[:count]


def _fortran_fmt(val):
    """Format a float in Fortran d-exponent notation."""
    if val == 0.0:
        return '0.d0'
    s = f'{val:.10e}'.replace('e', 'd')
    return s


class _FamilyWidget(QGroupBox):
    """Spinboxes for one extra ion family."""

    def __init__(self, index, parent=None):
        super().__init__(f"Family {index}", parent)
        lay = QFormLayout(self)
        self.nbrion_spin = _make_int_spin(1, 1, 999999)
        self.amat_spin = _make_double_spin(100.0, 0.0, decimals=8)
        self.a_mass_spin = _make_double_spin(100.0, 0.0, decimals=8)
        self.qn_spin = _make_double_spin(22.0, 0.0, decimals=2)
        self.r0_spin = _make_double_spin(0.033, 0.0, decimals=12)
        self.teta0_spin = _make_double_spin(266.9, decimals=8)
        self.z0_spin = _make_double_spin(0.0, decimals=6)
        self.pr0_spin = _make_double_spin(426.28, decimals=8)
        self.pz0_spin = _make_double_spin(0.0, decimals=6)
        self.w0mev_spin = _make_double_spin(5.665e-3, 0.0, decimals=8)
        self.phi0_spin = _make_double_spin(-250.0, decimals=4)
        lay.addRow("Ions :", self.nbrion_spin)
        lay.addRow("Atomic mass :", self.amat_spin)
        lay.addRow("Mass number A :", self.a_mass_spin)
        lay.addRow("Charge Q :", self.qn_spin)
        lay.addRow("r0 (m) :", self.r0_spin)
        lay.addRow("teta0 (deg) :", self.teta0_spin)
        lay.addRow("z0 (m) :", self.z0_spin)
        lay.addRow("pr0 (mrad) :", self.pr0_spin)
        lay.addRow("pz0 (mrad) :", self.pz0_spin)
        lay.addRow("W0 (MeV/A) :", self.w0mev_spin)
        lay.addRow("phi0 (deg) :", self.phi0_spin)

    def get_values(self):
        return {
            'nbrion': self.nbrion_spin.value(),
            'amat': self.amat_spin.value(),
            'a_mass': self.a_mass_spin.value(),
            'qn': self.qn_spin.value(),
            'r0': self.r0_spin.value(),
            'teta0': self.teta0_spin.value(),
            'z0': self.z0_spin.value(),
            'pr0': self.pr0_spin.value(),
            'pz0': self.pz0_spin.value(),
            'w0mev': self.w0mev_spin.value(),
            'phi0': self.phi0_spin.value(),
        }

    def set_values(self, d):
        self.nbrion_spin.setValue(d.get('nbrion', 1))
        self.amat_spin.setValue(d.get('amat', 100.0))
        self.a_mass_spin.setValue(d.get('a_mass', 100.0))
        self.qn_spin.setValue(d.get('qn', 22.0))
        self.r0_spin.setValue(d.get('r0', 0.033))
        self.teta0_spin.setValue(d.get('teta0', 266.9))
        self.z0_spin.setValue(d.get('z0', 0.0))
        self.pr0_spin.setValue(d.get('pr0', 426.28))
        self.pz0_spin.setValue(d.get('pz0', 0.0))
        self.w0mev_spin.setValue(d.get('w0mev', 5.665e-3))
        self.phi0_spin.setValue(d.get('phi0', -250.0))


class InputFormPanel(QWidget):
    """Form panel with every LIONS simulation input."""

    params_changed = pyqtSignal()
    file_loaded = pyqtSignal(str)  # path of the data file just loaded

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(280)
        self.setMaximumWidth(400)
        self._loading = False  # guards against re-entrant loads
        self._loaded_deck = None

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        container = QWidget()
        main_layout = QVBoxLayout(container)
        main_layout.setSpacing(5)

        # ============================================================
        # Section 1: General
        # ============================================================
        self.sec_general = CollapsibleSection("General")
        main_layout.addWidget(self.sec_general)

        self.repswap_combo = QComboBox()
        self.repswap_combo.addItems(['y', 'n'])
        self.repswap_combo.setCurrentText('y')
        self.sec_general.add_row("Byte swap :", self.repswap_combo)

        file_widget = QWidget()
        file_layout = QHBoxLayout(file_widget)
        file_layout.setContentsMargins(0, 0, 0, 0)
        self.ficin_edit = QLineEdit('don.120.ref_phi10_mars_2017')
        btn_browse = QPushButton('...')
        btn_browse.setFixedWidth(30)
        btn_browse.clicked.connect(self._browse_ficin)
        file_layout.addWidget(self.ficin_edit)
        file_layout.addWidget(btn_browse)
        self.sec_general.add_row("Data file :", file_widget)

        btn_load_data = QPushButton("Load data file")
        btn_load_data.clicked.connect(
            lambda: self.load_data_file(self.ficin_edit.text())
        )
        self.sec_general.add_row("", btn_load_data)

        # ============================================================
        # Section 2: Description
        # ============================================================
        self.sec_desc = CollapsibleSection("Description", collapsed=True)
        main_layout.addWidget(self.sec_desc)

        self.title_edit = QLineEdit()
        self.sec_desc.add_row("Title :", self.title_edit)

        # ============================================================
        # Section 3: RF / Cyclotron
        # ============================================================
        self.sec_rf = CollapsibleSection("RF / Cyclotron", collapsed=True)
        main_layout.addWidget(self.sec_rf)

        self.npaspt_spin = _make_int_spin(1080, 1, 99999)
        self.fhf_spin = _make_double_spin(14400000.0, 0.0, decimals=1)
        self.pothf_spin = _make_double_spin(-51890.0, decimals=1)
        self.teta1_spin = _make_double_spin(20.0, 0.0, 180.0, decimals=4)
        self.harm_spin = _make_double_spin(3.0, 0.0, 20.0, decimals=1)
        self.dbsb_spin = _make_double_spin(1.0e-4, decimals=8)
        self.sec_rf.add_row("Steps/turn (npaspt) :", self.npaspt_spin)
        self.sec_rf.add_row("RF freq (Hz) :", self.fhf_spin)
        self.sec_rf.add_row("RF voltage (V) :", self.pothf_spin)
        self.sec_rf.add_row("Dee angle (deg) :", self.teta1_spin)
        self.sec_rf.add_row("Harmonic :", self.harm_spin)
        self.sec_rf.add_row("dB/B :", self.dbsb_spin)

        # ============================================================
        # Section 4: Radius / Turns
        # ============================================================
        self.sec_radius = CollapsibleSection("Radius / Turns", collapsed=True)
        main_layout.addWidget(self.sec_radius)

        self.r1th_spin = _make_double_spin(1.50, 0.0, decimals=6)
        self.ntours_spin = _make_int_spin(250, 1, 99999)
        self.sec_radius.add_row("Final radius (m) :", self.r1th_spin)
        self.sec_radius.add_row("Turns :", self.ntours_spin)

        # ============================================================
        # Section 5: Dees
        # ============================================================
        self.sec_dees = CollapsibleSection("Dees", collapsed=True)
        main_layout.addWidget(self.sec_dees)

        self.dphid_spin = _make_double_spin(0.0, decimals=4)
        self.dphig_spin = _make_double_spin(0.0, decimals=4)
        self.pothfd_spin = _make_double_spin(51890.0, 0.0, decimals=1)
        self.pothfg_spin = _make_double_spin(51890.0, 0.0, decimals=1)
        self.sec_dees.add_row("dphid (deg) :", self.dphid_spin)
        self.sec_dees.add_row("dphig (deg) :", self.dphig_spin)
        self.sec_dees.add_row("Right dee voltage (V) :", self.pothfd_spin)
        self.sec_dees.add_row("Left dee voltage (V) :", self.pothfg_spin)

        # ============================================================
        # Section 6: Harmonics
        # ============================================================
        self.sec_harmonics = CollapsibleSection("Harmonics", collapsed=True)
        main_layout.addWidget(self.sec_harmonics)

        self.ncoh0_spin = _make_int_spin(0, 0, 10)
        self.sec_harmonics.add_row("ncoh channel 0 :", self.ncoh0_spin)
        self.coh0_detail = QWidget()
        coh0_lay = QFormLayout(self.coh0_detail)
        coh0_lay.setContentsMargins(0, 0, 0, 0)
        self.coh0_edit = QLineEdit()
        self.coh0_edit.setPlaceholderText("space-separated coefficients")
        coh0_lay.addRow("Coeff coh[0] :", self.coh0_edit)
        self.phihar0_spin = _make_double_spin(0.0, decimals=4)
        coh0_lay.addRow("Phase (deg) :", self.phihar0_spin)
        self.sec_harmonics.content_layout.addRow(self.coh0_detail)
        self.coh0_detail.setVisible(False)
        self.ncoh0_spin.valueChanged.connect(
            lambda v: self.coh0_detail.setVisible(v != 0)
        )

        sep_h = QFrame()
        sep_h.setFrameShape(QFrame.HLine)
        self.sec_harmonics.content_layout.addRow(sep_h)

        self.ncoh1_spin = _make_int_spin(0, 0, 10)
        self.sec_harmonics.add_row("ncoh channel 1 :", self.ncoh1_spin)
        self.coh1_detail = QWidget()
        coh1_lay = QFormLayout(self.coh1_detail)
        coh1_lay.setContentsMargins(0, 0, 0, 0)
        self.coh1_edit = QLineEdit()
        self.coh1_edit.setPlaceholderText("space-separated coefficients")
        coh1_lay.addRow("Coeff coh[1] :", self.coh1_edit)
        self.phihar1_spin = _make_double_spin(0.0, decimals=4)
        coh1_lay.addRow("Phase (deg) :", self.phihar1_spin)
        self.sec_harmonics.content_layout.addRow(self.coh1_detail)
        self.coh1_detail.setVisible(False)
        self.ncoh1_spin.valueChanged.connect(
            lambda v: self.coh1_detail.setVisible(v != 0)
        )

        # ============================================================
        # Section 7: Z geometry
        # ============================================================
        self.sec_geomz = CollapsibleSection("Z geometry", collapsed=True)
        main_layout.addWidget(self.sec_geomz)

        self.zmin_spin = _make_double_spin(-0.5, decimals=4)
        self.zmax_spin = _make_double_spin(0.5, decimals=4)
        self.rmax1_spin = _make_double_spin(0.164, 0.0, decimals=6)
        self.sec_geomz.add_row("zmin (m) :", self.zmin_spin)
        self.sec_geomz.add_row("zmax (m) :", self.zmax_spin)
        self.sec_geomz.add_row("rmax1 (m) :", self.rmax1_spin)

        # ============================================================
        # Section 8: Extraction
        # ============================================================
        self.sec_extract = CollapsibleSection("Extraction", collapsed=True)
        main_layout.addWidget(self.sec_extract)

        self.rex_spin = _make_double_spin(-1.0, decimals=4)
        self.tet1ex_spin = _make_double_spin(-1.0, decimals=4)
        self.tet2ex_spin = _make_double_spin(-1.0, decimals=4)
        self.eex_spin = _make_double_spin(-1.0, decimals=4)
        self.sec_extract.add_row("rex (m) :", self.rex_spin)
        self.sec_extract.add_row("tet1ex (deg) :", self.tet1ex_spin)
        self.sec_extract.add_row("tet2ex (deg) :", self.tet2ex_spin)
        self.sec_extract.add_row("eex :", self.eex_spin)

        # Shown only when rex < 0
        self.extract_detail = QWidget()
        ext_lay = QFormLayout(self.extract_detail)
        ext_lay.setContentsMargins(0, 0, 0, 0)

        self.ndf_spin = _make_int_spin(0, 0, 10)
        ext_lay.addRow("Deflectors (count) :", self.ndf_spin)
        self.edf_spin = _make_double_spin(0.0, decimals=0)
        self.edf_spin.setRange(-1e9, 1e9)
        self.edf_spin.setSingleStep(10000)
        self.edf_spin.setSuffix(" V")
        ext_lay.addRow("Extraction voltage :", self.edf_spin)
        self.deflectors_text = QPlainTextEdit()
        self.deflectors_text.setPlaceholderText(
            "1 line per deflector: rdf edf xdfi ydfi xdfo ydfo "
            "tetdfi_deg tetdfo_deg"
        )
        self.deflectors_text.setMaximumHeight(80)
        ext_lay.addRow("Deflectors :", self.deflectors_text)

        sep_e = QFrame()
        sep_e.setFrameShape(QFrame.HLine)
        ext_lay.addRow(sep_e)

        self.ndm_spin = _make_int_spin(0, 0, 10)
        ext_lay.addRow("Mag. dipoles (count) :", self.ndm_spin)
        self.magdip_text = QPlainTextEdit()
        self.magdip_text.setPlaceholderText(
            "1 line per dipole: rdm bdm xdm ydm tetdmi_deg tetdmo_deg"
        )
        self.magdip_text.setMaximumHeight(60)
        ext_lay.addRow("Mag. dipoles :", self.magdip_text)

        # Shown only when tet2ex > 0
        self.grids_detail = QWidget()
        grid_lay = QFormLayout(self.grids_detail)
        grid_lay.setContentsMargins(0, 0, 0, 0)
        self.ngex_spin = _make_int_spin(0, 0, 10)
        grid_lay.addRow("Grids (count) :", self.ngex_spin)
        self.grids_text = QPlainTextEdit()
        self.grids_text.setPlaceholderText(
            "1 line per grid: tetagex xgex_mm ygex_mm xgexc_mm ygexc_mm"
        )
        self.grids_text.setMaximumHeight(60)
        grid_lay.addRow("Grids :", self.grids_text)
        self.ficdebgex_edit = QLineEdit()
        grid_lay.addRow("Grid base file :", self.ficdebgex_edit)
        ext_lay.addRow(self.grids_detail)
        self.grids_detail.setVisible(False)

        self.sec_extract.content_layout.addRow(self.extract_detail)
        self.extract_detail.setVisible(False)

        self.rex_spin.valueChanged.connect(self._on_rex_changed)
        self.tet2ex_spin.valueChanged.connect(self._on_tet2ex_changed)

        # ============================================================
        # Section 9: Ion(s)
        # ============================================================
        self.sec_ions = CollapsibleSection("Ion(s)", collapsed=True)
        main_layout.addWidget(self.sec_ions)

        self.nbrfam_spin = _make_int_spin(1, 1, 10)
        self.sec_ions.add_row("Families :", self.nbrfam_spin)

        fam1_group = QGroupBox("Family 1")
        fam1_lay = QFormLayout(fam1_group)
        self.nbrion_spin = _make_int_spin(1, 1, 999999)
        self.amat_spin = _make_double_spin(100.0, 0.0, decimals=8)
        self.a_mass_spin = _make_double_spin(100.0, 0.0, decimals=8)
        self.qn_spin = _make_double_spin(22.0, 0.0, decimals=2)
        self.r0_spin = _make_double_spin(0.033, 0.0, decimals=12)
        self.teta0_spin = _make_double_spin(266.9, decimals=8)
        self.z0_spin = _make_double_spin(0.0, decimals=6)
        self.pr0_spin = _make_double_spin(426.28, decimals=8)
        self.pz0_spin = _make_double_spin(0.0, decimals=6)
        self.w0mev_spin = _make_double_spin(5.665e-3, 0.0, decimals=8)
        self.phi0_spin = _make_double_spin(-250.0, decimals=4)
        fam1_lay.addRow("Ions :", self.nbrion_spin)
        fam1_lay.addRow("Atomic mass :", self.amat_spin)
        fam1_lay.addRow("Mass number A :", self.a_mass_spin)
        fam1_lay.addRow("Charge Q :", self.qn_spin)
        fam1_lay.addRow("r0 (m) :", self.r0_spin)
        fam1_lay.addRow("teta0 (deg) :", self.teta0_spin)
        fam1_lay.addRow("z0 (m) :", self.z0_spin)
        fam1_lay.addRow("pr0 (mrad) :", self.pr0_spin)
        fam1_lay.addRow("pz0 (mrad) :", self.pz0_spin)
        fam1_lay.addRow("W0 (MeV/A) :", self.w0mev_spin)
        fam1_lay.addRow("phi0 (deg) :", self.phi0_spin)
        self.sec_ions.content_layout.addRow(fam1_group)

        # Shown only when nbrfam > 1
        self.n_total_widget = QWidget()
        ntot_lay = QFormLayout(self.n_total_widget)
        ntot_lay.setContentsMargins(0, 0, 0, 0)
        self.n_total_spin = _make_int_spin(1, 1, 999999)
        ntot_lay.addRow("Total particles :", self.n_total_spin)
        self.sec_ions.content_layout.addRow(self.n_total_widget)
        self.n_total_widget.setVisible(False)

        self.extra_fam_container = QWidget()
        self.extra_fam_layout = QVBoxLayout(self.extra_fam_container)
        self.extra_fam_layout.setContentsMargins(0, 0, 0, 0)
        self.extra_fam_layout.setSpacing(4)
        self.sec_ions.content_layout.addRow(self.extra_fam_container)
        self.extra_families = []  # list of _FamilyWidget

        self.nbrfam_spin.valueChanged.connect(self._on_nbrfam_changed)

        # ============================================================
        # Section 10: Correlations
        # ============================================================
        self.sec_correlations = CollapsibleSection(
            "Correlations", collapsed=True
        )
        main_layout.addWidget(self.sec_correlations)

        self.iflcor_spin = _make_int_spin(0, 0, 10)
        self.sec_correlations.add_row("iflcor :", self.iflcor_spin)

        self.correl_detail = QWidget()
        cor_lay = QFormLayout(self.correl_detail)
        cor_lay.setContentsMargins(0, 0, 0, 0)
        w_ficsig, self.ficsig_edit, _btn = _make_file_row(
            '', lambda: self._browse_open(self.ficsig_edit)
        )
        cor_lay.addRow("Sigma file :", w_ficsig)
        self.phase_space_edit = QLineEdit('6d')
        cor_lay.addRow("Phase space :", self.phase_space_edit)

        sep_sig = QFrame()
        sep_sig.setFrameShape(QFrame.HLine)
        cor_lay.addRow(sep_sig)
        lbl_sig = QLabel("Diagonal \u03c3 (no corr.) :")
        lbl_sig.setStyleSheet("font-weight: bold;")
        cor_lay.addRow(lbl_sig)
        self.sig_dr = _make_double_spin(0.1, 0.0, 1000.0, decimals=6)
        self.sig_ds = _make_double_spin(1.0, 0.0, 1e6, decimals=6)
        self.sig_dz = _make_double_spin(0.1, 0.0, 1000.0, decimals=6)
        self.sig_dpz = _make_double_spin(1.0, 0.0, 1e6, decimals=6)
        self.sig_dphi = _make_double_spin(0.1, 0.0, 1000.0, decimals=6)
        self.sig_dp = _make_double_spin(1.0, 0.0, 1e6, decimals=6)
        cor_lay.addRow("\u03c3 dr (mm) :", self.sig_dr)
        cor_lay.addRow("\u03c3 ds (mrad) :", self.sig_ds)
        cor_lay.addRow("\u03c3 dz (mm) :", self.sig_dz)
        cor_lay.addRow("\u03c3 dpz (mrad) :", self.sig_dpz)
        cor_lay.addRow("\u03c3 d\u03c6 (mm) :", self.sig_dphi)
        cor_lay.addRow("\u03c3 dp (mrad) :", self.sig_dp)
        btn_apply_sig = QPushButton("Apply \u2192 sigma file")
        btn_apply_sig.setToolTip(
            "Write the diagonal values to the sigma file "
            "(overwrites the existing file)"
        )
        btn_apply_sig.clicked.connect(self._write_sigma_from_gui)
        cor_lay.addRow("", btn_apply_sig)

        sep_sig2 = QFrame()
        sep_sig2.setFrameShape(QFrame.HLine)
        cor_lay.addRow(sep_sig2)
        self.ficcor_edit = QLineEdit('no')
        cor_lay.addRow("Correlation file :", self.ficcor_edit)
        self.ficmat_edit = QLineEdit('no')
        cor_lay.addRow("Matrix file :", self.ficmat_edit)
        self.ficdc_edit = QLineEdit('no')
        cor_lay.addRow("dc file :", self.ficdc_edit)
        self.sec_correlations.content_layout.addRow(self.correl_detail)
        self.correl_detail.setVisible(False)
        self.iflcor_spin.valueChanged.connect(
            lambda v: self.correl_detail.setVisible(v != 0)
        )

        # ============================================================
        # Section 11: Simulation
        # ============================================================
        self.sec_simulation = CollapsibleSection("Simulation", collapsed=True)
        main_layout.addWidget(self.sec_simulation)

        self.nbdt_spin = _make_int_spin(350000, 1, 99999999)
        self.ifreq_spin = _make_int_spin(30, 1, 999999)
        self.sec_simulation.add_row("Steps (nbdt) :", self.nbdt_spin)
        self.sec_simulation.add_row("Display freq :", self.ifreq_spin)

        # ============================================================
        # Section 12: Output files
        # ============================================================
        self.sec_output_files = CollapsibleSection(
            "Output files", collapsed=True
        )
        main_layout.addWidget(self.sec_output_files)

        self.ficz_edit = QLineEdit('no')
        self.sec_output_files.add_row("z file :", self.ficz_edit)
        self.ficxy_edit = QLineEdit('no')
        self.sec_output_files.add_row("xy file :", self.ficxy_edit)

        # ============================================================
        # Section 13: Field files
        # ============================================================
        self.sec_field_files = CollapsibleSection(
            "Field files", collapsed=True
        )
        main_layout.addWidget(self.sec_field_files)

        w_ficb, self.ficb_edit, _btn = _make_file_row(
            '', lambda: self._browse_ficb()
        )
        self.sec_field_files.add_row("Magnetic map (ficb) :", w_ficb)
        w_ficv1, self.ficv1_edit, _btn = _make_file_row(
            '', lambda: self._browse_open(self.ficv1_edit)
        )
        self.sec_field_files.add_row("Central potential (ficv1) :", w_ficv1)
        w_ficv2, self.ficv2_edit, _btn = _make_file_row(
            '', lambda: self._browse_open(self.ficv2_edit)
        )
        self.sec_field_files.add_row("Outer potential (ficv2) :", w_ficv2)

        self.b_factor_spin = _make_double_spin(1.0, 0.0, 100.0, decimals=6)
        self.b_factor_spin.setSuffix(" T")
        self.b_factor_spin.setToolTip(
            "B factor (tesla). The magnetic map is normalised to 1,\n"
            "then multiplied by this factor. Example: 1.56 T for CIME."
        )
        self.sec_field_files.add_row("B factor (T) :", self.b_factor_spin)

        self.field_interp_combo = QComboBox()
        self.field_interp_combo.addItems(OPTIONS['interpolation'][1])
        self.field_interp_combo.setToolTip(
            "Magnetic map interpolation. bilinear: the Fortran model.\n"
            "hermite: bicubic Hermite surface, C1 across cells, "
            "about 1.8 times the field cost\n"
            "(benchmarks/DESIGN.md, case H)."
        )
        self.sec_field_files.add_row(
            "Field interpolation :", self.field_interp_combo
        )

        # ============================================================
        # Section 14: Vertical sorter
        # ============================================================
        self.sec_vertical = CollapsibleSection("Vertical sorter")
        main_layout.addWidget(self.sec_vertical)

        self.rmin_ver = _make_double_spin(1000.0, decimals=1)
        self.rmax_ver = _make_double_spin(2000.0, decimals=1)
        self.teta1_ver = _make_double_spin(0.0, decimals=2)
        self.teta2_ver = _make_double_spin(0.0, decimals=2)
        self.sec_vertical.add_row("rmin_ver :", self.rmin_ver)
        self.sec_vertical.add_row("rmax_ver :", self.rmax_ver)
        self.sec_vertical.add_row("teta1_ver (deg) :", self.teta1_ver)
        self.sec_vertical.add_row("teta2_ver (deg) :", self.teta2_ver)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        self.sec_vertical.content_layout.addRow(sep)

        self.pot_ver = _make_double_spin(0.0)
        self.h_ver = _make_double_spin(0.0)
        self.phi_ver = _make_double_spin(0.0, decimals=2)
        self.anuz = _make_double_spin(0.0)
        self.coupnuz = _make_double_spin(0.0)
        self.dpotspot = _make_double_spin(0.0)
        self.sec_vertical.add_row("pot_ver :", self.pot_ver)
        self.sec_vertical.add_row("h_ver :", self.h_ver)
        self.sec_vertical.add_row("phi_ver (deg) :", self.phi_ver)
        self.sec_vertical.add_row("anuz :", self.anuz)
        self.sec_vertical.add_row("coupnuz :", self.coupnuz)
        self.sec_vertical.add_row("dpotspot :", self.dpotspot)

        # ============================================================
        # Section 15: Bidim / Diagnostics
        # ============================================================
        self.sec_bidim = CollapsibleSection("Bidim / Diagnostics")
        main_layout.addWidget(self.sec_bidim)

        self.numbid_spin = _make_int_spin(0, 0, 999)
        self.sec_bidim.add_row("First bidim number :", self.numbid_spin)
        self.numbid_spin.valueChanged.connect(self._on_numbid_changed)

        self.bidim_detail = QWidget()
        bidim_layout = QFormLayout(self.bidim_detail)
        bidim_layout.setContentsMargins(0, 0, 0, 0)

        self.ifreb_spin = _make_int_spin(1, 1, 999999)
        bidim_layout.addRow("Plot freq :", self.ifreb_spin)

        self.baryc_combo = QComboBox()
        self.baryc_combo.addItems(['c', 'b'])
        bidim_layout.addRow("Centres/Barycentres (c/b) :", self.baryc_combo)

        self.bid_xmin = _make_double_spin(-0.5, decimals=4)
        self.bid_xmax = _make_double_spin(0.5, decimals=4)
        self.bid_ymin = _make_double_spin(-0.5, decimals=4)
        self.bid_ymax = _make_double_spin(0.5, decimals=4)
        bidim_layout.addRow("xmin (m) :", self.bid_xmin)
        bidim_layout.addRow("xmax (m) :", self.bid_xmax)
        bidim_layout.addRow("ymin (m) :", self.bid_ymin)
        bidim_layout.addRow("ymax (m) :", self.bid_ymax)

        self.multi_group = QGroupBox("Correlations (n > 1)")
        multi_layout = QFormLayout(self.multi_group)
        self.dphi_max = _make_double_spin(10.0, 0.0, decimals=2)
        self.dwsw_max = _make_double_spin(0.01, 0.0, decimals=6)
        self.dr_max = _make_double_spin(5.0, 0.0, decimals=2)
        self.dpr_max = _make_double_spin(0.01, 0.0, decimals=6)
        self.dz_max = _make_double_spin(5.0, 0.0, decimals=2)
        self.dpz_max = _make_double_spin(0.01, 0.0, decimals=6)
        self.itebid_spin = _make_int_spin(1, 1, 999999)
        multi_layout.addRow("dphi_max :", self.dphi_max)
        multi_layout.addRow("dw/w_max :", self.dwsw_max)
        multi_layout.addRow("dr_max :", self.dr_max)
        multi_layout.addRow("dpr/p0_max :", self.dpr_max)
        multi_layout.addRow("dz_max :", self.dz_max)
        multi_layout.addRow("dpz/p0_max :", self.dpz_max)
        multi_layout.addRow("Bidim start iteration :", self.itebid_spin)
        bidim_layout.addRow(self.multi_group)

        self.bary_xmin = _make_double_spin(-0.2, decimals=4)
        self.bary_xmax = _make_double_spin(0.2, decimals=4)
        self.bary_ymin = _make_double_spin(-0.2, decimals=4)
        self.bary_ymax = _make_double_spin(0.2, decimals=4)
        bidim_layout.addRow("Bary xmin (m) :", self.bary_xmin)
        bidim_layout.addRow("Bary xmax (m) :", self.bary_xmax)
        bidim_layout.addRow("Bary ymin (m) :", self.bary_ymin)
        bidim_layout.addRow("Bary ymax (m) :", self.bary_ymax)

        ficrep_widget = QWidget()
        ficrep_layout = QHBoxLayout(ficrep_widget)
        ficrep_layout.setContentsMargins(0, 0, 0, 0)
        self.ficrep_edit = QLineEdit('no')
        btn_ficrep = QPushButton('...')
        btn_ficrep.setFixedWidth(30)
        btn_ficrep.clicked.connect(self._browse_ficrep)
        ficrep_layout.addWidget(self.ficrep_edit)
        ficrep_layout.addWidget(btn_ficrep)
        bidim_layout.addRow("Marker file :", ficrep_widget)

        self.sec_bidim.content_layout.addRow(self.bidim_detail)
        self.bidim_detail.setVisible(False)

        # ============================================================
        # Section 16: Post-simulation
        # ============================================================
        self.sec_post = CollapsibleSection("Post-Simulation")
        main_layout.addWidget(self.sec_post)

        self.save_trapped_combo = QComboBox()
        self.save_trapped_combo.addItems(['n', 'o'])
        self.sec_post.add_row("Save trapped (o/n) :", self.save_trapped_combo)

        ficpie_widget = QWidget()
        ficpie_layout = QHBoxLayout(ficpie_widget)
        ficpie_layout.setContentsMargins(0, 0, 0, 0)
        self.ficpie_edit = QLineEdit('trapped.dat')
        btn_ficpie = QPushButton('...')
        btn_ficpie.setFixedWidth(30)
        btn_ficpie.clicked.connect(self._browse_ficpie)
        ficpie_layout.addWidget(self.ficpie_edit)
        ficpie_layout.addWidget(btn_ficpie)
        self.sec_post.add_row("Trapped file :", ficpie_widget)

        # ============================================================
        # Save / load buttons
        # ============================================================
        btn_layout = QHBoxLayout()
        btn_save = QPushButton("Save params")
        btn_load = QPushButton("Load params")
        btn_save_data = QPushButton("Save data file")
        btn_save.clicked.connect(self.save_params)
        btn_load.clicked.connect(self.load_params)
        btn_save_data.clicked.connect(self._save_data_file_dialog)
        btn_layout.addWidget(btn_save)
        btn_layout.addWidget(btn_load)
        btn_layout.addWidget(btn_save_data)
        main_layout.addLayout(btn_layout)

        main_layout.addStretch()
        scroll.setWidget(container)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(scroll)

    # ================================================================
    # Conditional visibility
    # ================================================================

    def _on_numbid_changed(self, val):
        self.bidim_detail.setVisible(val != 0)

    def _on_nbrfam_changed(self, val):
        self.n_total_widget.setVisible(val > 1)
        needed = val - 1  # families beyond the first
        while len(self.extra_families) < needed:
            idx = len(self.extra_families) + 2
            fw = _FamilyWidget(idx)
            self.extra_fam_layout.addWidget(fw)
            self.extra_families.append(fw)
        for i, fw in enumerate(self.extra_families):
            fw.setVisible(i < needed)
        # Keep n_total >= sum of nbrion
        if val > 1:
            total = self.nbrion_spin.value()
            for i in range(needed):
                total += self.extra_families[i].nbrion_spin.value()
            if self.n_total_spin.value() < total:
                self.n_total_spin.setValue(total)

    def _on_rex_changed(self, val):
        self.extract_detail.setVisible(val < 0)

    def _on_tet2ex_changed(self, val):
        self.grids_detail.setVisible(val > 0)

    # ================================================================
    # Browse callbacks
    # ================================================================

    def _browse_ficin(self):
        path, _ = QFileDialog.getOpenFileName(self, "Data file")
        if path:
            self.ficin_edit.setText(path)
            self.load_data_file(path)

    def _browse_ficrep(self):
        path, _ = QFileDialog.getOpenFileName(self, "Marker file")
        if path:
            self.ficrep_edit.setText(path)

    def _browse_ficpie(self):
        path, _ = QFileDialog.getSaveFileName(self, "Trapped particles file")
        if path:
            self.ficpie_edit.setText(path)

    def _compute_betagamma(self):
        """Total beta*gamma of family 1 from the form, as lions.io.load
        computes it; 0.0 if undefined.
        """
        amat = self.amat_spin.value()
        mass_number = self.a_mass_spin.value()
        w0mev = self.w0mev_spin.value()
        if amat <= 0 or w0mev <= 0:
            return 0.0
        rest_energy = (
            AMU * (amat - ELECTRON_MASS_U * self.qn_spin.value()) * C_LIGHT**2
        )
        gamma = 1.0 + w0mev * MEV * mass_number / rest_energy
        return math.sqrt(gamma**2 - 1.0)

    def _load_sigma_diag(self, ficsig, data_path):
        """Load the 6 diagonal values of a sigma file into the
        spinboxes.

        The file holds sigma^2 (variance) in m^2 / (bg)^2. The GUI shows
        sigma (std dev) in mm (positions) / mrad (momenta).
        """
        try:
            # Relative paths are resolved against the data file
            # directory
            if data_path and os.path.dirname(data_path):
                sig_path = os.path.join(os.path.dirname(data_path), ficsig)
            else:
                sig_path = ficsig
            if not os.path.exists(sig_path):
                sig_path = ficsig
            if not os.path.exists(sig_path):
                return
            bg_total = self._compute_betagamma()
            with open(sig_path, 'r') as sf:
                spins = [
                    self.sig_dr,
                    self.sig_ds,
                    self.sig_dz,
                    self.sig_dpz,
                    self.sig_dphi,
                    self.sig_dp,
                ]
                # True = position (m -> mm), False = momentum (bg ->
                # mrad)
                is_pos = [True, False, True, False, True, False]
                for i, spin in enumerate(spins):
                    line = sf.readline()
                    if not line:
                        break
                    vals = line.replace('d', 'e').replace('D', 'E').split()
                    if len(vals) > i:
                        sig2 = float(vals[i])
                        sig = math.sqrt(max(sig2, 0.0))
                        if is_pos[i]:
                            sig *= 1000.0  # m -> mm
                        else:
                            # Small-angle approximation of arcsin(p_comp
                            # / p_total) * 1000
                            if bg_total > 0:
                                sig = sig / bg_total * 1000.0
                        spin.setValue(sig)
        except Exception:
            pass  # unreadable sigma file: keep the current values

    def _write_sigma_from_gui(self):
        """Write a diagonal 6x6 sigma file from the GUI values.

        The GUI holds sigma (std dev) in mm (positions) / mrad
        (momenta). The file holds sigma^2 (variance) in m^2 / (bg)^2.
        """
        ficsig = self.ficsig_edit.text().strip()
        if not ficsig:
            ficsig = 'sigma6d_no_corr.dat'
            self.ficsig_edit.setText(ficsig)
        ficin = self.ficin_edit.text().strip()
        if ficin and os.path.dirname(ficin):
            sig_path = os.path.join(os.path.dirname(ficin), ficsig)
        else:
            sig_path = ficsig
        try:
            bg_total = self._compute_betagamma()
            sig_gui = [
                self.sig_dr.value(),
                self.sig_ds.value(),
                self.sig_dz.value(),
                self.sig_dpz.value(),
                self.sig_dphi.value(),
                self.sig_dp.value(),
            ]
            # True = position (mm -> m), False = momentum (mrad -> bg)
            is_pos = [True, False, True, False, True, False]
            # Floor sigma^2 above 0: the Cholesky decomposition needs a
            # positive definite matrix
            eps_pos = 1e-30  # m^2 (~ 1e-12 mm)
            eps_mom = 1e-30  # (bg)^2
            sig2 = []
            for i, s in enumerate(sig_gui):
                if is_pos[i]:
                    s_m = s * 1e-3  # mm -> m
                    sig2.append(max(s_m * s_m, eps_pos))
                else:
                    # mrad -> bg : sig_bg = sig_mrad * 1e-3 * bg_total
                    s_bg = s * 1e-3 * bg_total if bg_total > 0 else s
                    sig2.append(max(s_bg * s_bg, eps_mom))
            with open(sig_path, 'w') as f:
                for i in range(6):
                    row = []
                    for j in range(6):
                        if i == j:
                            row.append(_fortran_fmt(sig2[i]))
                        else:
                            row.append('0.d0')
                    f.write('   '.join(row) + '\n')
            print(f'Sigma file written: {sig_path}')
        except Exception as e:
            QMessageBox.warning(self, "Sigma error", str(e))

    def _browse_open(self, edit_widget):
        path, _ = QFileDialog.getOpenFileName(self, "Open file")
        if path:
            edit_widget.setText(path)

    def _browse_ficb(self):
        """Browse for the magnetic field map, then detect Bmax."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open magnetic map",
            "",
            "Field maps (*.dat *.lio.dat);;All (*)",
        )
        if path:
            self.ficb_edit.setText(path)
            self._auto_detect_bmax(path)

    def _auto_detect_bmax(self, filepath):
        """Read the magnetic field map and set the B factor to its
        Bmax."""
        try:
            with open(filepath, 'r') as f:
                f.readline()  # comment line
                vals = f.readline().split()
                nr = int(vals[2])
                ntet = int(vals[3])
                # First block is Bz
                bz_vals = []
                while len(bz_vals) < nr * ntet:
                    line = f.readline()
                    if not line:
                        break
                    for tok in line.split():
                        bz_vals.append(
                            float(tok.replace('d', 'e').replace('D', 'e'))
                        )
                if bz_vals:
                    bmax = max(abs(v) for v in bz_vals[: nr * ntet])
                    self.b_factor_spin.setValue(bmax)
                    print(f'Bmax detected: {bmax:.6e} T, B factor updated')
        except Exception as e:
            print(f'Bmax detection failed: {e}')

    def _save_data_file_dialog(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save data file", "", "All (*)"
        )
        if path:
            self.write_data_file(path)

    # ================================================================
    # Data file loading
    # ================================================================

    def load_data_file(self, path):
        """Parse the data file and fill the GUI fields."""
        if self._loading:
            return
        if not path or not os.path.exists(path):
            return

        self._loading = True
        try:
            self._do_load_data_file(path)
            self.file_loaded.emit(path)
        except Exception as e:
            QMessageBox.warning(
                self, "Load error", f"Error while loading {path}:\n{e}"
            )
        finally:
            self._loading = False

    def _do_load_data_file(self, path):
        """Fill the form from a deck read by lions.deck.read_deck."""
        deck = read_deck(path)
        # Kept for the entries the form does not edit (gap offset, bunch
        # files)
        self._loaded_deck = deck

        self.title_edit.setText(deck.title)
        self.npaspt_spin.setValue(deck.steps_per_turn)
        self.r1th_spin.setValue(deck.design_r_final)
        self.ntours_spin.setValue(deck.design_turns)
        self.fhf_spin.setValue(deck.rf_frequency)
        self.pothf_spin.setValue(deck.rf_voltage)
        self.teta1_spin.setValue(deck.dee_half_angle_deg)
        self.harm_spin.setValue(deck.harmonic)
        self.dbsb_spin.setValue(deck.db_over_b)
        self.dphid_spin.setValue(deck.phase_right_deg)
        self.dphig_spin.setValue(deck.phase_left_deg)
        self.pothfd_spin.setValue(deck.v_dee_right)
        self.pothfg_spin.setValue(deck.v_dee_left)
        for h, n_spin, coef_edit, phase_spin in (
            (
                deck.harmonic1,
                self.ncoh0_spin,
                self.coh0_edit,
                self.phihar0_spin,
            ),
            (
                deck.harmonic2,
                self.ncoh1_spin,
                self.coh1_edit,
                self.phihar1_spin,
            ),
        ):
            n_spin.setValue(0 if h is None else len(h.coefs))
            coef_edit.setText('' if h is None else _numbers(h.coefs))
            phase_spin.setValue(0.0 if h is None else h.phase_deg)

        self.zmin_spin.setValue(deck.zmin)
        self.zmax_spin.setValue(deck.zmax)
        self.rmax1_spin.setValue(deck.central_rmax)
        self.rex_spin.setValue(deck.extraction)
        self.tet1ex_spin.setValue(deck.cut_bunch)
        self.tet2ex_spin.setValue(deck.grids_on)
        self.eex_spin.setValue(deck.grid_start)
        self.ndf_spin.setValue(len(deck.deflectors))
        self.edf_spin.setValue(
            deck.deflectors[0].e if deck.deflectors else 0.0
        )
        self.deflectors_text.setPlainText(
            '\n'.join(_numbers(astuple(e)) for e in deck.deflectors)
        )
        self.ndm_spin.setValue(len(deck.dipoles))
        self.magdip_text.setPlainText(
            '\n'.join(_numbers(astuple(m)) for m in deck.dipoles)
        )
        if deck.grids:
            self.ngex_spin.setValue(len(deck.grids))
            self.grids_text.setPlainText(
                '\n'.join(_numbers(astuple(g)) for g in deck.grids)
            )
            self.ficdebgex_edit.setText(deck.grid_file_prefix)

        self.nbrfam_spin.setValue(len(deck.families))
        if deck.n_particles is not None:
            self.n_total_spin.setValue(deck.n_particles)
        if len(self.extra_families) < len(deck.families) - 1:
            self._on_nbrfam_changed(len(deck.families))
        self._set_family_values(self, deck.families[0])
        for widget, fam in zip(self.extra_families, deck.families[1:]):
            self._set_family_values(widget, fam)

        self.iflcor_spin.setValue(deck.correlations)
        if deck.correlations != 0:
            self.ficsig_edit.setText(deck.sigma_file)
            self.phase_space_edit.setText(deck.phase_space)
            self.ficcor_edit.setText(deck.correlation_file)
            self.ficmat_edit.setText(deck.matrix_file)
            self.ficdc_edit.setText(deck.dc_matrix_file)
            self._load_sigma_diag(deck.sigma_file, path)

        self.nbdt_spin.setValue(deck.n_steps)
        self.ifreq_spin.setValue(deck.output_every)
        self.ficz_edit.setText(deck.z_output_file)
        self.ficxy_edit.setText(deck.xy_output_file)
        self.ficb_edit.setText(deck.bfield_file)
        self.ficv1_edit.setText(deck.central_potential_file)
        self.ficv2_edit.setText(deck.outer_potential_file)
        self.field_interp_combo.setCurrentText(deck.field_interpolation)
        if deck.b_factor is not None:
            self.b_factor_spin.setValue(deck.b_factor)
        else:
            self.b_factor_spin.setValue(1.0)
            ficb_path = self.ficb_edit.text().strip()
            if ficb_path and ficb_path != 'no' and os.path.isfile(ficb_path):
                self._auto_detect_bmax(ficb_path)

    @staticmethod
    def _set_family_values(widget, fam):
        widget.nbrion_spin.setValue(fam.size)
        widget.amat_spin.setValue(fam.atomic_mass)
        widget.a_mass_spin.setValue(fam.mass_number)
        widget.qn_spin.setValue(fam.charge_state)
        widget.r0_spin.setValue(fam.r0)
        widget.teta0_spin.setValue(fam.theta0_deg)
        widget.z0_spin.setValue(fam.z0)
        widget.pr0_spin.setValue(fam.pr0_mrad)
        widget.pz0_spin.setValue(fam.pz0_mrad)
        widget.w0mev_spin.setValue(fam.w0_mev_per_u)
        widget.phi0_spin.setValue(fam.phi0_deg)

    @staticmethod
    def _family(widget):
        return Family(
            widget.nbrion_spin.value(),
            widget.amat_spin.value(),
            widget.a_mass_spin.value(),
            widget.qn_spin.value(),
            widget.r0_spin.value(),
            widget.teta0_spin.value(),
            widget.z0_spin.value(),
            widget.pr0_spin.value(),
            widget.pz0_spin.value(),
            widget.w0mev_spin.value(),
            widget.phi0_spin.value(),
        )

    # ================================================================
    # Data file writing
    # ================================================================

    def to_deck(self):
        """Deck built from the form values."""
        harmonics = []
        for n_spin, coef_edit, phase_spin in (
            (self.ncoh0_spin, self.coh0_edit, self.phihar0_spin),
            (self.ncoh1_spin, self.coh1_edit, self.phihar1_spin),
        ):
            n = n_spin.value()
            coefs = (parse_floats(coef_edit.text()) + [0.0] * n)[:n]
            harmonics.append(
                FieldHarmonic(coefs, phase_spin.value()) if n else None
            )

        deflectors, dipoles, grids = [], [], []
        grid_prefix = ''
        if self.rex_spin.value() < 0:
            rows = _rows(
                self.deflectors_text.toPlainText(), 8, self.ndf_spin.value()
            )
            if rows:
                rows[0][1] = (
                    self.edf_spin.value()
                )  # the field box overrides deflector 1
            deflectors = [Deflector(*row) for row in rows]
            dipoles = [
                Dipole(*row)
                for row in _rows(
                    self.magdip_text.toPlainText(), 6, self.ndm_spin.value()
                )
            ]
            if self.tet2ex_spin.value() > 0:
                grids = [
                    Grid(*row)
                    for row in _rows(
                        self.grids_text.toPlainText(),
                        5,
                        self.ngex_spin.value(),
                    )
                ]
                grid_prefix = self.ficdebgex_edit.text().strip() or 'grex'

        n_families = self.nbrfam_spin.value()
        families = [self._family(self)]
        families += [
            self._family(w) for w in self.extra_families[: n_families - 1]
        ]
        families += [self._family(self)] * (n_families - len(families))

        kept = self._loaded_deck
        return Deck(
            title=self.title_edit.text(),
            steps_per_turn=self.npaspt_spin.value(),
            design_r_final=self.r1th_spin.value(),
            design_turns=self.ntours_spin.value(),
            rf_frequency=self.fhf_spin.value(),
            rf_voltage=self.pothf_spin.value(),
            dee_half_angle_deg=self.teta1_spin.value(),
            harmonic=self.harm_spin.value(),
            db_over_b=self.dbsb_spin.value(),
            phase_right_deg=self.dphid_spin.value(),
            phase_left_deg=self.dphig_spin.value(),
            v_dee_right=self.pothfd_spin.value(),
            v_dee_left=self.pothfg_spin.value(),
            zmin=self.zmin_spin.value(),
            zmax=self.zmax_spin.value(),
            central_rmax=self.rmax1_spin.value(),
            families=families,
            n_particles=self.n_total_spin.value() if n_families > 1 else None,
            n_steps=self.nbdt_spin.value(),
            output_every=self.ifreq_spin.value(),
            z_output_file=self.ficz_edit.text().strip(),
            xy_output_file=self.ficxy_edit.text().strip(),
            bfield_file=self.ficb_edit.text().strip() or 'no',
            central_potential_file=self.ficv1_edit.text().strip() or 'no',
            outer_potential_file=self.ficv2_edit.text().strip() or 'no',
            harmonic1=harmonics[0],
            harmonic2=harmonics[1],
            gap_offset=kept.gap_offset if kept else None,
            extraction=self.rex_spin.value(),
            cut_bunch=self.tet1ex_spin.value(),
            grids_on=self.tet2ex_spin.value(),
            grid_start=self.eex_spin.value(),
            deflectors=deflectors,
            dipoles=dipoles,
            grids=grids,
            grid_file_prefix=grid_prefix,
            correlations=self.iflcor_spin.value(),
            sigma_file=self.ficsig_edit.text().strip() or 'no',
            phase_space=self.phase_space_edit.text().strip() or '6d',
            correlation_file=self.ficcor_edit.text().strip() or 'no',
            matrix_file=self.ficmat_edit.text().strip() or 'no',
            dc_matrix_file=self.ficdc_edit.text().strip() or 'no',
            b_factor=self.b_factor_spin.value(),
            bunch_load_file=kept.bunch_load_file if kept else 'no',
            bunch_save_file=kept.bunch_save_file if kept else 'no',
            field_interpolation=self.field_interp_combo.currentText(),
            field_azimuth=kept.field_azimuth if kept else 'arccos',
        )

    def write_data_file(self, path):
        """Write the form values as a deck (lions.deck.write_deck)."""
        write_deck(self.to_deck(), path)
        self._log_family_summary(self.nbrfam_spin.value())
        return path

    def _log_family_summary(self, nbrfam):
        """Print the main values of each ion family."""
        summary = [f"[write_data_file] {nbrfam} family(ies):"]
        summary.append(
            f"  Family 1: nbrion={self.nbrion_spin.value()}, "
            f"amat={self.amat_spin.value()}, qn={self.qn_spin.value()}, "
            f"w0={self.w0mev_spin.value()}"
        )
        for k in range(nbrfam - 1):
            if k < len(self.extra_families):
                fw = self.extra_families[k]
                summary.append(
                    f"  Family {k + 2}: nbrion={fw.nbrion_spin.value()}, "
                    f"amat={fw.amat_spin.value()}, qn={fw.qn_spin.value()}, "
                    f"w0={fw.w0mev_spin.value()}"
                )
        print('\n'.join(summary))

    # ================================================================
    # Deflector geometry for the XY plot
    # ================================================================

    def get_deflector_geometry(self):
        """Deflector geometry for the XY plot.

        Returns:
            list of dicts with keys: rdf_m, xdfi_m, ydfi_m, xdfo_m,
            ydfo_m, tetdfi_deg, tetdfo_deg
        """
        if self.rex_spin.value() >= 0:
            return []
        ndf = self.ndf_spin.value()
        if ndf <= 0:
            return []
        text = self.deflectors_text.toPlainText().strip()
        if not text:
            return []

        from ..fortran_text import fortran_float

        result = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                result.append(
                    {
                        'rdf_m': fortran_float(parts[0]),
                        'xdfi_m': fortran_float(parts[2]),
                        'ydfi_m': fortran_float(parts[3]),
                        'xdfo_m': fortran_float(parts[4]),
                        'ydfo_m': fortran_float(parts[5]),
                        'tetdfi_deg': fortran_float(parts[6]),
                        'tetdfo_deg': fortran_float(parts[7]),
                    }
                )
            except (ValueError, IndexError):
                continue
        return result

    # ================================================================
    # Answer sequences for the worker (input queue)
    # ================================================================

    def get_input_sequence(self):
        """Ordered answers fed to the engine in place of input().

        The GUI values go to a temporary data file next to ficin.
        """
        ficin_orig = self.ficin_edit.text()
        if ficin_orig:
            work_dir = os.path.dirname(os.path.abspath(ficin_orig))
        else:
            work_dir = '.'
        temp_path = os.path.join(work_dir, '_lions_gui_temp.dat')
        if self.iflcor_spin.value() != 0:
            self._write_sigma_from_gui()
        self.write_data_file(temp_path)

        answers = []
        # 1. repswap
        answers.append(self.repswap_combo.currentText())
        # 2. ficin (temporary file)
        answers.append(temp_path)
        # 3. rmin_ver rmax_ver teta1_ver teta2_ver
        answers.append(
            f"{self.rmin_ver.value()} {self.rmax_ver.value()} "
            f"{self.teta1_ver.value()} {self.teta2_ver.value()}"
        )
        # 4. pot_ver h_ver phi_ver anuz coupnuz dpotspot
        answers.append(
            f"{self.pot_ver.value()} {self.h_ver.value()} "
            f"{self.phi_ver.value()} {self.anuz.value()} "
            f"{self.coupnuz.value()} {self.dpotspot.value()}"
        )
        # 5. numbid
        answers.append(str(self.numbid_spin.value()))

        if self.numbid_spin.value() != 0:
            answers.append(str(self.ifreb_spin.value()))
            answers.append(self.baryc_combo.currentText())
            answers.append(
                f"{self.bid_xmin.value()} {self.bid_xmax.value()} "
                f"{self.bid_ymin.value()} {self.bid_ymax.value()}"
            )
            # Multi-particle correlations (initgra reads them when d.n >
            # 1)
            total_n = self.nbrion_spin.value()
            for fw in self.extra_families:
                total_n += fw.nbrion_spin.value()
            if total_n > 1:
                answers.append(
                    f"{self.dphi_max.value()} {self.dwsw_max.value()}"
                )
                answers.append(f"{self.dr_max.value()} {self.dpr_max.value()}")
                answers.append(f"{self.dz_max.value()} {self.dpz_max.value()}")
                answers.append(str(self.itebid_spin.value()))
            answers.append(
                f"{self.bary_xmin.value()} {self.bary_xmax.value()} "
                f"{self.bary_ymin.value()} {self.bary_ymax.value()}"
            )
            answers.append(self.ficrep_edit.text() or "no")

        return answers

    def get_post_answers(self):
        """Answers for the post-simulation phase (trapped particles)."""
        answers = []
        answers.append(self.save_trapped_combo.currentText())
        if self.save_trapped_combo.currentText() == 'o':
            answers.append(self.ficpie_edit.text())
        return answers

    # ================================================================
    # Save / load JSON parameters
    # ================================================================

    def save_params(self):
        """Ask for a JSON file and save the form parameters to it."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save parameters", "", "JSON (*.json)"
        )
        if not path:
            return
        data = self._collect_all_params()
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)

    def load_params(self):
        """Ask for a JSON file and apply its parameters to the form."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Load parameters", "", "JSON (*.json)"
        )
        if not path:
            return
        with open(path, 'r') as f:
            data = json.load(f)
        self._apply_all_params(data)

    def _collect_all_params(self):
        """All parameters as a dict (keys are the saved JSON format)."""
        return {
            # General
            'repswap': self.repswap_combo.currentText(),
            'ficin': self.ficin_edit.text(),
            # Description
            'corps': self.title_edit.text(),
            # RF
            'npaspt': self.npaspt_spin.value(),
            'fhf': self.fhf_spin.value(),
            'pothf': self.pothf_spin.value(),
            'teta1': self.teta1_spin.value(),
            'harm': self.harm_spin.value(),
            'dbsb': self.dbsb_spin.value(),
            # Radius / turns
            'r1th': self.r1th_spin.value(),
            'ntours': self.ntours_spin.value(),
            # Dees
            'dphid': self.dphid_spin.value(),
            'dphig': self.dphig_spin.value(),
            'pothfd': self.pothfd_spin.value(),
            'pothfg': self.pothfg_spin.value(),
            # Harmonics
            'ncoh0': self.ncoh0_spin.value(),
            'coh0': self.coh0_edit.text(),
            'phihar0': self.phihar0_spin.value(),
            'ncoh1': self.ncoh1_spin.value(),
            'coh1': self.coh1_edit.text(),
            'phihar1': self.phihar1_spin.value(),
            # Z geometry
            'zmin': self.zmin_spin.value(),
            'zmax': self.zmax_spin.value(),
            'rmax1': self.rmax1_spin.value(),
            # Extraction
            'rex': self.rex_spin.value(),
            'tet1ex': self.tet1ex_spin.value(),
            'tet2ex': self.tet2ex_spin.value(),
            'eex': self.eex_spin.value(),
            'ndf': self.ndf_spin.value(),
            'edf': self.edf_spin.value(),
            'deflectors': self.deflectors_text.toPlainText(),
            'ndm': self.ndm_spin.value(),
            'magdip': self.magdip_text.toPlainText(),
            'ngex': self.ngex_spin.value(),
            'grids': self.grids_text.toPlainText(),
            'ficdebgex': self.ficdebgex_edit.text(),
            # Ions
            'nbrfam': self.nbrfam_spin.value(),
            'n_total': self.n_total_spin.value(),
            'nbrion': self.nbrion_spin.value(),
            'amat': self.amat_spin.value(),
            'a_mass': self.a_mass_spin.value(),
            'qn': self.qn_spin.value(),
            'r0': self.r0_spin.value(),
            'teta0': self.teta0_spin.value(),
            'z0': self.z0_spin.value(),
            'pr0': self.pr0_spin.value(),
            'pz0': self.pz0_spin.value(),
            'w0mev': self.w0mev_spin.value(),
            'phi0': self.phi0_spin.value(),
            'extra_families': [
                fw.get_values() for fw in self.extra_families if fw.isVisible()
            ],
            # Correlations
            'iflcor': self.iflcor_spin.value(),
            'ficsig': self.ficsig_edit.text(),
            'espace': self.phase_space_edit.text(),
            'ficcor': self.ficcor_edit.text(),
            'ficmat': self.ficmat_edit.text(),
            'ficdc': self.ficdc_edit.text(),
            'sig_dr': self.sig_dr.value(),
            'sig_ds': self.sig_ds.value(),
            'sig_dz': self.sig_dz.value(),
            'sig_dpz': self.sig_dpz.value(),
            'sig_dphi': self.sig_dphi.value(),
            'sig_dp': self.sig_dp.value(),
            # Simulation
            'nbdt': self.nbdt_spin.value(),
            'ifreq': self.ifreq_spin.value(),
            # Output files
            'ficz': self.ficz_edit.text(),
            'ficxy': self.ficxy_edit.text(),
            # Field files
            'ficb': self.ficb_edit.text(),
            'ficv1': self.ficv1_edit.text(),
            'ficv2': self.ficv2_edit.text(),
            'b_factor': self.b_factor_spin.value(),
            'field_interpolation': self.field_interp_combo.currentText(),
            # Vertical sorter
            'rmin_ver': self.rmin_ver.value(),
            'rmax_ver': self.rmax_ver.value(),
            'teta1_ver': self.teta1_ver.value(),
            'teta2_ver': self.teta2_ver.value(),
            'pot_ver': self.pot_ver.value(),
            'h_ver': self.h_ver.value(),
            'phi_ver': self.phi_ver.value(),
            'anuz': self.anuz.value(),
            'coupnuz': self.coupnuz.value(),
            'dpotspot': self.dpotspot.value(),
            # Bidim
            'numbid': self.numbid_spin.value(),
            'ifreb': self.ifreb_spin.value(),
            'baryc': self.baryc_combo.currentText(),
            'bid_xmin': self.bid_xmin.value(),
            'bid_xmax': self.bid_xmax.value(),
            'bid_ymin': self.bid_ymin.value(),
            'bid_ymax': self.bid_ymax.value(),
            'dphi_max': self.dphi_max.value(),
            'dwsw_max': self.dwsw_max.value(),
            'dr_max': self.dr_max.value(),
            'dpr_max': self.dpr_max.value(),
            'dz_max': self.dz_max.value(),
            'dpz_max': self.dpz_max.value(),
            'itebid': self.itebid_spin.value(),
            'bary_xmin': self.bary_xmin.value(),
            'bary_xmax': self.bary_xmax.value(),
            'bary_ymin': self.bary_ymin.value(),
            'bary_ymax': self.bary_ymax.value(),
            'ficrep': self.ficrep_edit.text(),
            # Post-simulation
            'save_trapped': self.save_trapped_combo.currentText(),
            'ficpie': self.ficpie_edit.text(),
        }

    def _apply_all_params(self, data):
        """Apply a parameter dict to the GUI fields."""
        # General
        self.repswap_combo.setCurrentText(data.get('repswap', 'y'))
        self.ficin_edit.setText(data.get('ficin', ''))
        # Description
        self.title_edit.setText(data.get('corps', ''))
        # RF
        self.npaspt_spin.setValue(data.get('npaspt', 1080))
        self.fhf_spin.setValue(data.get('fhf', 14400000.0))
        self.pothf_spin.setValue(data.get('pothf', -51890.0))
        self.teta1_spin.setValue(data.get('teta1', 20.0))
        self.harm_spin.setValue(data.get('harm', 3.0))
        self.dbsb_spin.setValue(data.get('dbsb', 1e-4))
        # Radius / turns
        self.r1th_spin.setValue(data.get('r1th', 1.50))
        self.ntours_spin.setValue(data.get('ntours', 250))
        # Dees
        self.dphid_spin.setValue(data.get('dphid', 0.0))
        self.dphig_spin.setValue(data.get('dphig', 0.0))
        self.pothfd_spin.setValue(data.get('pothfd', 51890.0))
        self.pothfg_spin.setValue(data.get('pothfg', 51890.0))
        # Harmonics
        self.ncoh0_spin.setValue(data.get('ncoh0', 0))
        self.coh0_edit.setText(data.get('coh0', ''))
        self.phihar0_spin.setValue(data.get('phihar0', 0.0))
        self.ncoh1_spin.setValue(data.get('ncoh1', 0))
        self.coh1_edit.setText(data.get('coh1', ''))
        self.phihar1_spin.setValue(data.get('phihar1', 0.0))
        # Z geometry
        self.zmin_spin.setValue(data.get('zmin', -0.5))
        self.zmax_spin.setValue(data.get('zmax', 0.5))
        self.rmax1_spin.setValue(data.get('rmax1', 0.164))
        # Extraction
        self.rex_spin.setValue(data.get('rex', -1.0))
        self.tet1ex_spin.setValue(data.get('tet1ex', -1.0))
        self.tet2ex_spin.setValue(data.get('tet2ex', -1.0))
        self.eex_spin.setValue(data.get('eex', -1.0))
        self.ndf_spin.setValue(data.get('ndf', 0))
        self.edf_spin.setValue(data.get('edf', 0.0))
        self.deflectors_text.setPlainText(data.get('deflectors', ''))
        self.ndm_spin.setValue(data.get('ndm', 0))
        self.magdip_text.setPlainText(data.get('magdip', ''))
        self.ngex_spin.setValue(data.get('ngex', 0))
        self.grids_text.setPlainText(data.get('grids', ''))
        self.ficdebgex_edit.setText(data.get('ficdebgex', ''))
        # Ions
        self.nbrfam_spin.setValue(data.get('nbrfam', 1))
        self.n_total_spin.setValue(data.get('n_total', 1))
        self.nbrion_spin.setValue(data.get('nbrion', 1))
        self.amat_spin.setValue(data.get('amat', 100.0))
        self.a_mass_spin.setValue(data.get('a_mass', 100.0))
        self.qn_spin.setValue(data.get('qn', 22.0))
        self.r0_spin.setValue(data.get('r0', 0.033))
        self.teta0_spin.setValue(data.get('teta0', 266.9))
        self.z0_spin.setValue(data.get('z0', 0.0))
        self.pr0_spin.setValue(data.get('pr0', 426.28))
        self.pz0_spin.setValue(data.get('pz0', 0.0))
        self.w0mev_spin.setValue(data.get('w0mev', 5.665e-3))
        self.phi0_spin.setValue(data.get('phi0', -250.0))
        # Extra families (nbrfam_spin.setValue created the widgets)
        extra_fams = data.get('extra_families', [])
        for i, fd in enumerate(extra_fams):
            if i < len(self.extra_families):
                self.extra_families[i].set_values(fd)
        # Correlations
        self.iflcor_spin.setValue(data.get('iflcor', 0))
        self.ficsig_edit.setText(data.get('ficsig', ''))
        self.phase_space_edit.setText(data.get('espace', '6d'))
        self.ficcor_edit.setText(data.get('ficcor', 'no'))
        self.ficmat_edit.setText(data.get('ficmat', 'no'))
        self.ficdc_edit.setText(data.get('ficdc', 'no'))
        self.sig_dr.setValue(data.get('sig_dr', 0.1))
        self.sig_ds.setValue(data.get('sig_ds', 1.0))
        self.sig_dz.setValue(data.get('sig_dz', 0.1))
        self.sig_dpz.setValue(data.get('sig_dpz', 1.0))
        self.sig_dphi.setValue(data.get('sig_dphi', 0.1))
        self.sig_dp.setValue(data.get('sig_dp', 1.0))
        # Simulation
        self.nbdt_spin.setValue(data.get('nbdt', 350000))
        self.ifreq_spin.setValue(data.get('ifreq', 30))
        # Output files
        self.ficz_edit.setText(data.get('ficz', 'no'))
        self.ficxy_edit.setText(data.get('ficxy', 'no'))
        # Field files
        self.ficb_edit.setText(data.get('ficb', ''))
        self.ficv1_edit.setText(data.get('ficv1', ''))
        self.ficv2_edit.setText(data.get('ficv2', ''))
        self.b_factor_spin.setValue(data.get('b_factor', 1.0))
        self.field_interp_combo.setCurrentText(
            data.get('field_interpolation', 'bilinear')
        )
        # Vertical sorter
        self.rmin_ver.setValue(data.get('rmin_ver', 1000.0))
        self.rmax_ver.setValue(data.get('rmax_ver', 2000.0))
        self.teta1_ver.setValue(data.get('teta1_ver', 0.0))
        self.teta2_ver.setValue(data.get('teta2_ver', 0.0))
        self.pot_ver.setValue(data.get('pot_ver', 0.0))
        self.h_ver.setValue(data.get('h_ver', 0.0))
        self.phi_ver.setValue(data.get('phi_ver', 0.0))
        self.anuz.setValue(data.get('anuz', 0.0))
        self.coupnuz.setValue(data.get('coupnuz', 0.0))
        self.dpotspot.setValue(data.get('dpotspot', 0.0))
        # Bidim
        self.numbid_spin.setValue(data.get('numbid', 0))
        self.ifreb_spin.setValue(data.get('ifreb', 1))
        self.baryc_combo.setCurrentText(data.get('baryc', 'c'))
        self.bid_xmin.setValue(data.get('bid_xmin', -0.5))
        self.bid_xmax.setValue(data.get('bid_xmax', 0.5))
        self.bid_ymin.setValue(data.get('bid_ymin', -0.5))
        self.bid_ymax.setValue(data.get('bid_ymax', 0.5))
        self.dphi_max.setValue(data.get('dphi_max', 10.0))
        self.dwsw_max.setValue(data.get('dwsw_max', 0.01))
        self.dr_max.setValue(data.get('dr_max', 5.0))
        self.dpr_max.setValue(data.get('dpr_max', 0.01))
        self.dz_max.setValue(data.get('dz_max', 5.0))
        self.dpz_max.setValue(data.get('dpz_max', 0.01))
        self.itebid_spin.setValue(data.get('itebid', 1))
        self.bary_xmin.setValue(data.get('bary_xmin', -0.2))
        self.bary_xmax.setValue(data.get('bary_xmax', 0.2))
        self.bary_ymin.setValue(data.get('bary_ymin', -0.2))
        self.bary_ymax.setValue(data.get('bary_ymax', 0.2))
        self.ficrep_edit.setText(data.get('ficrep', 'no'))
        # Post-simulation
        self.save_trapped_combo.setCurrentText(data.get('save_trapped', 'n'))
        self.ficpie_edit.setText(data.get('ficpie', 'trapped.dat'))
