"""
lions_bidim_viewer.py - Visualiseur de fichiers bidim GANIL (.bid)
Peut etre utilise seul ou integre dans la GUI LIONS.
"""

import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QFileDialog, QComboBox, QLabel, QDialog, QListWidget,
    QSplitter
)
from PyQt5.QtCore import Qt


def read_bidim(filepath):
    """Lit un fichier bidim GANIL et retourne le header + les points (ix, iy).

    Returns:
        dict avec les cles:
            num, titre, labelx, labely, maxx, maxy,
            xmin, xmax, ymin, ymax, ix_array, iy_array
    """
    with open(filepath, 'rb') as f:
        # Header
        num = struct.unpack('i', f.read(4))[0]
        f.read(10 + 12 + 12 + 10 + 10)  # type, date, time, place, author
        titre = f.read(25).decode('ascii', errors='replace').strip()
        labelx = f.read(25).decode('ascii', errors='replace').strip()
        labely = f.read(25).decode('ascii', errors='replace').strip()
        maxx, maxy = struct.unpack('2i', f.read(8))
        f.read(10)  # format
        cxmi = f.read(20).decode('ascii', errors='replace').strip()
        cxma = f.read(20).decode('ascii', errors='replace').strip()
        cymi = f.read(20).decode('ascii', errors='replace').strip()
        cyma = f.read(20).decode('ascii', errors='replace').strip()
        f.read(120)  # padding

        try:
            xmin = float(cxmi)
            xmax = float(cxma)
            ymin = float(cymi)
            ymax = float(cyma)
        except ValueError:
            xmin, xmax, ymin, ymax = 0, maxx, 0, maxy

        # Lire les blocs de donnees (4096 int32 par bloc)
        end_marker = 9999 * 65536 + 9999
        ix_list = []
        iy_list = []

        while True:
            block_data = f.read(4096 * 4)
            if not block_data:
                break
            block = np.frombuffer(block_data, dtype=np.int32)
            for val in block:
                if val == end_marker:
                    break
                ix = val // 65536
                iy = val % 65536
                ix_list.append(ix)
                iy_list.append(iy)
            else:
                continue
            break

    return {
        'num': num,
        'titre': titre,
        'labelx': labelx,
        'labely': labely,
        'maxx': maxx,
        'maxy': maxy,
        'xmin': xmin,
        'xmax': xmax,
        'ymin': ymin,
        'ymax': ymax,
        'ix': np.array(ix_list, dtype=np.int32),
        'iy': np.array(iy_list, dtype=np.int32),
    }


def plot_bidim(ax, bidim_data, nbins=256, cmap='hot'):
    """Trace un bidim sur un axe matplotlib.

    Args:
        ax: axe matplotlib
        bidim_data: dict retourne par read_bidim()
        nbins: nombre de bins pour l'histogramme 2D
        cmap: colormap matplotlib
    """
    ax.clear()
    d = bidim_data

    if len(d['ix']) == 0:
        ax.set_title(f"{d['titre']} (vide)")
        return

    # Convertir indices en coordonnees physiques
    if d['maxx'] > 0:
        x = d['xmin'] + (d['xmax'] - d['xmin']) * d['ix'] / d['maxx']
        y = d['ymin'] + (d['ymax'] - d['ymin']) * d['iy'] / d['maxy']
    else:
        x = d['ix'].astype(float)
        y = d['iy'].astype(float)

    # Histogramme 2D
    h, xedges, yedges = np.histogram2d(x, y, bins=nbins)
    h = h.T  # pour que l'axe y soit correct

    # Afficher avec log scale pour mieux voir
    h_display = np.where(h > 0, np.log10(h + 1), 0)

    ax.imshow(h_display, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
              origin='lower', aspect='auto', cmap=cmap, interpolation='nearest')

    ax.set_xlabel(d['labelx'] if d['labelx'] else 'x')
    ax.set_ylabel(d['labely'] if d['labely'] else 'y')
    ax.set_title(f"Bidim {d['num']} : {d['titre']} ({len(d['ix'])} pts)")


class BidimViewerWidget(QWidget):
    """Widget Qt pour visualiser les fichiers bidim."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.bidim_files = []
        self.bidim_data = {}

        layout = QVBoxLayout(self)

        # Toolbar
        toolbar = QHBoxLayout()
        btn_open = QPushButton("Ouvrir bidim(s)...")
        btn_open.clicked.connect(self.open_files)
        toolbar.addWidget(btn_open)

        toolbar.addWidget(QLabel("Colormap:"))
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(['hot', 'viridis', 'plasma', 'inferno', 'jet', 'gray'])
        self.cmap_combo.currentTextChanged.connect(self.replot)
        toolbar.addWidget(self.cmap_combo)

        toolbar.addWidget(QLabel("Bins:"))
        self.bins_combo = QComboBox()
        self.bins_combo.addItems(['128', '256', '512', '1024'])
        self.bins_combo.setCurrentText('256')
        self.bins_combo.currentTextChanged.connect(self.replot)
        toolbar.addWidget(self.bins_combo)

        btn_refresh = QPushButton("Rafraichir")
        btn_refresh.clicked.connect(self.reload_and_replot)
        toolbar.addWidget(btn_refresh)

        toolbar.addStretch()
        layout.addLayout(toolbar)

        # Splitter : liste + canvas
        splitter = QSplitter(Qt.Horizontal)

        self.file_list = QListWidget()
        self.file_list.setMaximumWidth(200)
        self.file_list.currentRowChanged.connect(self.on_selection_changed)
        splitter.addWidget(self.file_list)

        # Canvas matplotlib
        canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(canvas_widget)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(8, 6), tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.nav_toolbar = NavigationToolbar2QT(self.canvas, self)
        canvas_layout.addWidget(self.nav_toolbar)
        canvas_layout.addWidget(self.canvas)
        splitter.addWidget(canvas_widget)

        splitter.setSizes([150, 600])
        layout.addWidget(splitter)

    def open_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Ouvrir fichier(s) bidim", "", "Bidim (*.bid);;Tous (*)")
        if paths:
            self.bidim_files = paths
            self.file_list.clear()
            self.bidim_data.clear()
            for p in paths:
                name = p.split('/')[-1].split('\\')[-1]
                self.file_list.addItem(name)
                self.bidim_data[name] = read_bidim(p)
            if self.file_list.count() > 0:
                self.file_list.setCurrentRow(0)

    def load_from_paths(self, paths):
        """Charge des fichiers bidim depuis une liste de chemins."""
        self.bidim_files = paths
        self.file_list.clear()
        self.bidim_data.clear()
        for p in paths:
            name = p.split('/')[-1].split('\\')[-1]
            self.file_list.addItem(name)
            try:
                self.bidim_data[name] = read_bidim(p)
            except Exception as e:
                print(f"Erreur lecture {p}: {e}")
        if self.file_list.count() > 0:
            self.file_list.setCurrentRow(0)

    def on_selection_changed(self, row):
        if row >= 0:
            self.replot()

    def replot(self):
        item = self.file_list.currentItem()
        if not item:
            return
        name = item.text()
        if name not in self.bidim_data:
            return
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        nbins = int(self.bins_combo.currentText())
        cmap = self.cmap_combo.currentText()
        plot_bidim(ax, self.bidim_data[name], nbins=nbins, cmap=cmap)
        self.canvas.draw()

    def reload_and_replot(self):
        """Recharge les fichiers depuis le disque et re-affiche."""
        for i, p in enumerate(self.bidim_files):
            name = p.split('/')[-1].split('\\')[-1]
            try:
                self.bidim_data[name] = read_bidim(p)
            except Exception:
                pass
        self.replot()


class BidimViewerDialog(QDialog):
    """Dialogue standalone pour visualiser les bidims."""

    def __init__(self, parent=None, paths=None):
        super().__init__(parent)
        self.setWindowTitle("LIONS - Visualiseur Bidim")
        self.setMinimumSize(900, 600)

        layout = QVBoxLayout(self)
        self.viewer = BidimViewerWidget(self)
        layout.addWidget(self.viewer)

        if paths:
            self.viewer.load_from_paths(paths)


# ---------------------------------------------------------------
# Utilisation standalone
# ---------------------------------------------------------------
if __name__ == '__main__':
    import sys
    from PyQt5.QtWidgets import QApplication

    # Mode simple : afficher avec matplotlib directement
    if len(sys.argv) > 1 and '--simple' in sys.argv:
        files = [f for f in sys.argv[1:] if f != '--simple']
        for fpath in files:
            data = read_bidim(fpath)
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            plot_bidim(ax, data)
        plt.tight_layout()
        plt.show()
    else:
        # Mode Qt
        app = QApplication(sys.argv)
        dlg = BidimViewerDialog()
        if len(sys.argv) > 1:
            dlg.viewer.load_from_paths(sys.argv[1:])
        dlg.show()
        sys.exit(app.exec_())
