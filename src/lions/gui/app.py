"""
lions_gui.py - Interface graphique PyQt5 pour la simulation LIONS
Simulateur de dynamique d'ions dans un cyclotron.
"""

import sys
import os
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QDockWidget, QTextEdit,
    QToolBar, QAction, QProgressBar, QMessageBox,
    QInputDialog, QStatusBar, QFileDialog
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QTextCursor

from .bidim_viewer import BidimViewerDialog
from .field_viewer import FieldViewerDialog
from .input_form import InputFormPanel
from .plots import HAS_GL, Plot3DWindow, PlotWidget
from .worker import SimulationWorker

# Style clair-bleu
APP_STYLE = """
QMainWindow {
    background-color: #e8edf2;
}
QToolBar {
    background: #d0d8e8;
    border-bottom: 1px solid #aab5c8;
    spacing: 6px;
    padding: 4px 8px;
}
QToolBar QToolButton {
    color: #1a1a2e;
    background: #f0f4fa;
    border: 1px solid #aab5c8;
    border-radius: 4px;
    padding: 6px 14px;
    font-weight: bold;
    font-size: 11px;
    min-width: 60px;
}
QToolBar QToolButton:hover {
    background: #c0d0f0;
    border: 1px solid #6688bb;
}
QToolBar QToolButton:pressed {
    background: #a0b8e0;
}
QToolBar QToolButton:disabled {
    color: #999999;
    background: #dde2ea;
    border: 1px solid #c0c8d4;
}
QToolBar::separator {
    width: 2px;
    background: #aab5c8;
    margin: 4px 8px;
}
QProgressBar {
    background: #f0f4fa;
    border: 1px solid #aab5c8;
    border-radius: 4px;
    text-align: center;
    color: #1a1a2e;
    font-weight: bold;
    font-size: 10px;
    min-height: 20px;
}
QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #2196F3, stop:1 #1565C0);
    border-radius: 3px;
}
QDockWidget {
    color: #1a1a2e;
    titlebar-close-icon: none;
    font-weight: bold;
}
QDockWidget::title {
    background: #c0cce0;
    padding: 6px;
    border: 1px solid #aab5c8;
    color: #1a1a2e;
}
QStatusBar {
    background: #d0d8e8;
    color: #1a1a2e;
    font-weight: bold;
    font-size: 10px;
    border-top: 1px solid #aab5c8;
}
"""


class LionsMainWindow(QMainWindow):
    """Fenetre principale de l'interface LIONS."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("LIONS - Simulateur Cyclotron")
        self.setMinimumSize(1200, 800)
        self.setStyleSheet(APP_STYLE)
        self.worker = None
        self._pending_data = []
        self._3d_window = None

        self._setup_toolbar()
        self._setup_input_panel()
        self._setup_plot_area()
        self._setup_console()
        self._setup_statusbar()
        self._setup_plot_timer()

    # ----------------------------------------------------------------
    # Setup de l'interface
    # ----------------------------------------------------------------

    def _setup_toolbar(self):
        toolbar = QToolBar("Simulation")
        toolbar.setMovable(False)
        toolbar.setIconSize(toolbar.iconSize())
        self.addToolBar(toolbar)

        self.action_run = QAction("\u25b6 Run", self)
        self.action_run.setToolTip("Lancer la simulation")
        self.action_run.triggered.connect(self.on_run)
        toolbar.addAction(self.action_run)

        self.action_stop = QAction("\u25a0 Stop", self)
        self.action_stop.setToolTip("Arreter la simulation")
        self.action_stop.setEnabled(False)
        self.action_stop.triggered.connect(self.on_stop)
        toolbar.addAction(self.action_stop)

        self.action_reset = QAction("\u21ba Reset", self)
        self.action_reset.setToolTip("Reinitialiser les plots et la console")
        self.action_reset.triggered.connect(self.on_reset)
        toolbar.addAction(self.action_reset)

        toolbar.addSeparator()

        self.action_unzoom = QAction("\u2922 Unzoom", self)
        self.action_unzoom.setToolTip("Reinitialiser le zoom sur tous les plots")
        self.action_unzoom.triggered.connect(self.on_unzoom)
        toolbar.addAction(self.action_unzoom)

        self.action_export = QAction("\u2913 Exporter", self)
        self.action_export.setToolTip("Sauvegarder les plots en PNG")
        self.action_export.triggered.connect(self.on_export_plots)
        toolbar.addAction(self.action_export)

        self.action_bidim = QAction("\u25a6 Bidims", self)
        self.action_bidim.setToolTip("Ouvrir le visualiseur de fichiers bidim (.bid)")
        self.action_bidim.triggered.connect(self.on_view_bidims)
        toolbar.addAction(self.action_bidim)

        self.action_fields = QAction("\u2261 Champs", self)
        self.action_fields.setToolTip("Visualiser les champs magnetiques et electriques")
        self.action_fields.triggered.connect(self.on_view_fields)
        toolbar.addAction(self.action_fields)

        if HAS_GL:
            self.action_3d = QAction("\u2b1a 3D", self)
            self.action_3d.setToolTip("Ouvrir la vue 3D des trajectoires (X, Y, Z)")
            self.action_3d.triggered.connect(self.on_view_3d)
            toolbar.addAction(self.action_3d)

        toolbar.addSeparator()

        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximumWidth(400)
        self.progress_bar.setMinimumWidth(200)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("%p%")
        toolbar.addWidget(self.progress_bar)

    def _setup_input_panel(self):
        self.input_panel = InputFormPanel()
        self.input_panel.file_loaded.connect(self._on_data_file_loaded)
        dock = QDockWidget("Parametres", self)
        dock.setWidget(self.input_panel)
        dock.setFeatures(QDockWidget.DockWidgetMovable)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

    def _setup_plot_area(self):
        self.plot_widget = PlotWidget()
        self.setCentralWidget(self.plot_widget)

    def _setup_console(self):
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setFont(QFont("Consolas", 9))
        self.console.setMaximumHeight(200)
        self.console.setStyleSheet(
            "QTextEdit { background-color: #1e1e2e; color: #e0e0e0; "
            "border: 1px solid #aab5c8; }")

        dock = QDockWidget("Console", self)
        dock.setWidget(self.console)
        dock.setFeatures(QDockWidget.DockWidgetMovable)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)

    def _setup_statusbar(self):
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
        self.statusbar.showMessage("Pret")

    def _setup_plot_timer(self):
        """Timer continu pour rafraichir plots + console (10 FPS)."""
        self._plot_timer = QTimer()
        self._plot_timer.setInterval(100)
        self._plot_timer.timeout.connect(self._flush_all)
        self._pending_console = []

    # ----------------------------------------------------------------
    # Actions
    # ----------------------------------------------------------------

    def on_run(self):
        """Lance la simulation."""
        if self.worker and self.worker.isRunning():
            QMessageBox.warning(self, "Attention",
                                "Une simulation est deja en cours.")
            return

        # Verifier que les fichiers obligatoires sont renseignes
        ficin = self.input_panel.ficin_edit.text().strip()
        if not ficin:
            QMessageBox.warning(self, "Fichier manquant",
                                "Veuillez selectionner un fichier de donnees (ficin).")
            return
        ficb = self.input_panel.ficb_edit.text().strip()
        ficv2 = self.input_panel.ficv2_edit.text().strip()
        if not ficb or not ficv2:
            missing = []
            if not ficb:
                missing.append("ficb (carte magnetique)")
            if not ficv2:
                missing.append("ficv2 (potentiel peripherie)")
            QMessageBox.warning(
                self, "Fichiers champs manquants",
                "Les fichiers suivants sont obligatoires :\n- "
                + "\n- ".join(missing)
                + "\n\nChargez un fichier de donnees ou renseignez-les manuellement.")
            return

        # Construire la liste des reponses
        answers = self.input_panel.get_input_sequence()

        # Ajouter les reponses post-simulation
        answers.extend(self.input_panel.get_post_answers())

        # Determiner le repertoire de travail
        ficin = self.input_panel.ficin_edit.text()
        work_dir = os.path.dirname(os.path.abspath(ficin)) if ficin else None

        # Creer le worker
        self.worker = SimulationWorker(answers, work_dir)
        self.worker.progress.connect(self.on_progress)
        self.worker.iteration_data.connect(self.on_iteration_data)
        self.worker.console_output.connect(self.on_console_output)
        self.worker.simulation_finished.connect(self.on_finished)
        self.worker.simulation_error.connect(self.on_error)
        self.worker.input_needed.connect(self.on_input_needed)

        # UI state
        self.action_run.setEnabled(False)
        self.action_stop.setEnabled(True)
        self.progress_bar.setValue(0)
        self.statusbar.showMessage("Simulation en cours...")

        # Dessiner les deflecteurs sur le plot XY
        self._draw_deflectors()

        # Demarrer le timer de rafraichissement continu
        self._plot_timer.start()

        self.worker.start()

    def on_stop(self):
        """Arrete la simulation."""
        if self.worker:
            self.worker.stop_flag = True
            self.statusbar.showMessage("Arret en cours...")

    def on_reset(self):
        """Reinitialise les plots et la console."""
        self.plot_widget.canvas.clear_plots()
        if self._3d_window is not None:
            self._3d_window.clear_3d()
        self.console.clear()
        self.progress_bar.setValue(0)
        self.statusbar.showMessage("Pret")

    def on_unzoom(self):
        """Reinitialise le zoom sur tous les plots."""
        self.plot_widget.toolbar.home()
        self.statusbar.showMessage("Zoom reinitialise")

    def _on_data_file_loaded(self, path):
        """Affiche un message quand un fichier de donnees est charge."""
        name = os.path.basename(path)
        self.statusbar.showMessage(f"Fichier de donnees charge : {name}")
        # Dessiner les deflecteurs sur le plot XY
        self._draw_deflectors()

    def _draw_deflectors(self):
        """Dessine les deflecteurs sur le plot XY et la vue 3D."""
        defl = self.input_panel.get_deflector_geometry()
        self.plot_widget.canvas.draw_deflectors(defl)
        if self._3d_window is not None:
            self._3d_window.draw_deflectors_3d(defl)

    def on_view_3d(self):
        """Ouvre ou met au premier plan la fenetre 3D."""
        if self._3d_window is None:
            self._3d_window = Plot3DWindow(self.plot_widget.canvas, parent=self)
        self._3d_window.show()
        self._3d_window.raise_()
        self._3d_window.activateWindow()
        # Rafraichir avec les donnees existantes
        self._3d_window.refresh_3d()
        # Dessiner les deflecteurs
        defl = self.input_panel.get_deflector_geometry()
        self._3d_window.draw_deflectors_3d(defl)

    def on_export_plots(self):
        """Exporte les plots en image."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Exporter les plots", "lions_plots.png",
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if path:
            self.plot_widget.canvas.fig.savefig(path, dpi=150, bbox_inches='tight')
            self.statusbar.showMessage(f"Plots exportes : {path}")

    def on_view_bidims(self):
        """Ouvre le visualiseur de fichiers bidim."""
        import glob
        # Chercher les .bid dans le repertoire du fichier de donnees
        ficin = self.input_panel.ficin_edit.text()
        search_dir = os.path.dirname(os.path.abspath(ficin)) if ficin else '.'
        bid_files = sorted(glob.glob(os.path.join(search_dir, '*.bid')))
        dlg = BidimViewerDialog(self, paths=bid_files if bid_files else None)
        dlg.show()

    def on_view_fields(self):
        """Ouvre le visualiseur de champs magnetiques et electriques."""
        ficin = self.input_panel.ficin_edit.text()
        work_dir = os.path.dirname(os.path.abspath(ficin)) if ficin else '.'
        swap = self.input_panel.repswap_combo.currentText()

        # Extraire ficb, ficv1, ficv2 depuis le fichier de donnees
        ficb, ficv1, ficv2 = '', '', ''
        if ficin and os.path.exists(ficin):
            try:
                from ..io import FortranFile
                with open(ficin, 'r') as f:
                    ff = FortranFile(f)
                    ff.readline()                    # titre
                    ff.readline()                    # n
                    ff.readline()                    # r0 nr
                    ff.readline()                    # parametres physiques
                    ff.readline()                    # ncoh
                    ff.readline()                    # coh
                    ff.readline()                    # nbrpie
                    ff.readline()                    # ndef
                    ff.readline()                    # positions initiales
                    ff.readline()                    # rex etc
                    nbrfam_line = ff.readline()       # nbrfam
                    nbrfam = int(nbrfam_line.split()[0])
                    for _ in range(nbrfam):
                        ff.readline()                # famille
                    ff.readline()                    # ncoh2 (harmoniques)
                    ff.readline()                    # nbdt ifreq
                    ff.read_tokens(2)                # ficz ficxy
                    vals = ff.read_tokens(3)         # ficb ficv1 ficv2
                    ficb_name = vals[0]
                    ficv1_name = vals[1]
                    ficv2_name = vals[2]
                    # Construire les chemins absolus
                    if ficb_name and ficb_name != 'no':
                        ficb = os.path.join(work_dir, ficb_name)
                    if ficv1_name and ficv1_name != 'no':
                        ficv1 = os.path.join(work_dir, ficv1_name)
                    if ficv2_name and ficv2_name != 'no':
                        ficv2 = os.path.join(work_dir, ficv2_name)
            except Exception:
                pass

        dlg = FieldViewerDialog(self, ficb=ficb, ficv1=ficv1, ficv2=ficv2, swap=swap)
        dlg.show()

    # ----------------------------------------------------------------
    # Slots pour les signaux du worker
    # ----------------------------------------------------------------

    def on_progress(self, current, total):
        if total > 0:
            self.progress_bar.setValue(int(100 * current / total))

    def on_iteration_data(self, data):
        self._pending_data.append(data)

    def on_console_output(self, text):
        self._pending_console.append(text)

    def _flush_all(self):
        """Rafraichit plots + console en un seul coup (appele par timer)."""
        # Plots
        if self._pending_data:
            canvas = self.plot_widget.canvas
            for d in self._pending_data:
                canvas.append_data(d)
            self._pending_data.clear()
            canvas._refresh_all_curves()
            # Mettre a jour la vue 3D si ouverte
            if self._3d_window is not None and self._3d_window.isVisible():
                self._3d_window.refresh_3d()
        # Console
        if self._pending_console:
            text = ''.join(self._pending_console)
            self._pending_console.clear()
            self.console.moveCursor(QTextCursor.End)
            self.console.insertPlainText(text)
            self.console.moveCursor(QTextCursor.End)

    def on_finished(self):
        self._plot_timer.stop()
        self._flush_all()
        self.action_run.setEnabled(True)
        self.action_stop.setEnabled(False)
        self.progress_bar.setValue(100)
        self.statusbar.showMessage("Simulation terminee")

    def on_error(self, msg):
        self._plot_timer.stop()
        self._flush_all()
        self.action_run.setEnabled(True)
        self.action_stop.setEnabled(False)
        self.statusbar.showMessage("Erreur !")
        QMessageBox.critical(self, "Erreur de simulation", msg)

    def on_input_needed(self, prompt):
        """Gere un input() inattendu pendant la simulation."""
        # Proposer des valeurs par defaut selon le prompt
        defaults = {
            'dphi': f"{self.input_panel.dphi_max.value()} {self.input_panel.dwsw_max.value()}",
            'dr-dpr': f"{self.input_panel.dr_max.value()} {self.input_panel.dpr_max.value()}",
            'dz-dpz': f"{self.input_panel.dz_max.value()} {self.input_panel.dpz_max.value()}",
            'debut_iteration': str(self.input_panel.itebid_spin.value()),
            'iteration max': str(1000),
            'flag de choix': '0',
            'pieges': 'n',
            'nom du fichier': self.input_panel.ficpie_edit.text(),
        }
        default_val = ''
        for key, val in defaults.items():
            if key in prompt.lower():
                default_val = val
                break
        text, ok = QInputDialog.getText(
            self, "Entree requise", prompt, text=default_val)
        if ok and self.worker:
            self.worker.provide_input(text)
        elif self.worker:
            self.worker.provide_input(default_val or "0")

    # ----------------------------------------------------------------
    # Fermeture
    # ----------------------------------------------------------------

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.worker.stop_flag = True
            self.worker.wait(5000)
        event.accept()


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("LIONS Cyclotron Simulator")
    app.setStyle("Fusion")

    window = LionsMainWindow()
    window.show()

    sys.exit(app.exec_())

