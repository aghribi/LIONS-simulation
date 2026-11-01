"""Main window of the LIONS GUI: ion dynamics in a cyclotron."""

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
    """Main window of the LIONS GUI."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("LIONS - Cyclotron Simulator")
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

    def _setup_toolbar(self):
        toolbar = QToolBar("Simulation")
        toolbar.setMovable(False)
        toolbar.setIconSize(toolbar.iconSize())
        self.addToolBar(toolbar)

        self.action_run = QAction("\u25b6 Run", self)
        self.action_run.setToolTip("Start the simulation")
        self.action_run.triggered.connect(self.on_run)
        toolbar.addAction(self.action_run)

        self.action_stop = QAction("\u25a0 Stop", self)
        self.action_stop.setToolTip("Stop the simulation")
        self.action_stop.setEnabled(False)
        self.action_stop.triggered.connect(self.on_stop)
        toolbar.addAction(self.action_stop)

        self.action_reset = QAction("\u21ba Reset", self)
        self.action_reset.setToolTip("Clear the plots and the console")
        self.action_reset.triggered.connect(self.on_reset)
        toolbar.addAction(self.action_reset)

        toolbar.addSeparator()

        self.action_unzoom = QAction("\u2922 Unzoom", self)
        self.action_unzoom.setToolTip("Reset the zoom on all plots")
        self.action_unzoom.triggered.connect(self.on_unzoom)
        toolbar.addAction(self.action_unzoom)

        self.action_export = QAction("\u2913 Export", self)
        self.action_export.setToolTip("Save the plots as PNG, PDF or SVG")
        self.action_export.triggered.connect(self.on_export_plots)
        toolbar.addAction(self.action_export)

        self.action_bidim = QAction("\u25a6 Bidims", self)
        self.action_bidim.setToolTip("Open the bidim file viewer (.bid)")
        self.action_bidim.triggered.connect(self.on_view_bidims)
        toolbar.addAction(self.action_bidim)

        self.action_fields = QAction("\u2261 Fields", self)
        self.action_fields.setToolTip("View the magnetic and electric fields")
        self.action_fields.triggered.connect(self.on_view_fields)
        toolbar.addAction(self.action_fields)

        self.action_arcade = QAction("Arcade", self)
        self.action_arcade.setToolTip("ION RUSH: the cyclotron trainer (needs the game extra)")
        self.action_arcade.triggered.connect(self.on_open_arcade)
        toolbar.addAction(self.action_arcade)

        if HAS_GL:
            self.action_3d = QAction("\u2b1a 3D", self)
            self.action_3d.setToolTip("Open the 3D trajectory view (X, Y, Z)")
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
        dock = QDockWidget("Parameters", self)
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
        self.statusbar.showMessage("Ready")

    def _setup_plot_timer(self):
        """Refresh plots and console at 10 Hz while a run is active."""
        self._plot_timer = QTimer()
        self._plot_timer.setInterval(100)
        self._plot_timer.timeout.connect(self._flush_all)
        self._pending_console = []

    def on_run(self):
        """Start the simulation."""
        if self.worker and self.worker.isRunning():
            QMessageBox.warning(self, "Warning",
                                "A simulation is already running.")
            return

        ficin = self.input_panel.ficin_edit.text().strip()
        if not ficin:
            QMessageBox.warning(self, "Missing file",
                                "Select a data file (ficin).")
            return
        ficb = self.input_panel.ficb_edit.text().strip()
        ficv2 = self.input_panel.ficv2_edit.text().strip()
        if not ficb or not ficv2:
            missing = []
            if not ficb:
                missing.append("ficb (magnetic field map)")
            if not ficv2:
                missing.append("ficv2 (outer potential)")
            QMessageBox.warning(
                self, "Missing field files",
                "The following files are required:\n- "
                + "\n- ".join(missing)
                + "\n\nLoad a data file or fill them in by hand.")
            return

        answers = self.input_panel.get_input_sequence()
        answers.extend(self.input_panel.get_post_answers())

        ficin = self.input_panel.ficin_edit.text()
        work_dir = os.path.dirname(os.path.abspath(ficin)) if ficin else None

        self.worker = SimulationWorker(answers, work_dir)
        self.worker.progress.connect(self.on_progress)
        self.worker.iteration_data.connect(self.on_iteration_data)
        self.worker.console_output.connect(self.on_console_output)
        self.worker.simulation_finished.connect(self.on_finished)
        self.worker.simulation_error.connect(self.on_error)
        self.worker.input_needed.connect(self.on_input_needed)

        self.action_run.setEnabled(False)
        self.action_stop.setEnabled(True)
        self.progress_bar.setValue(0)
        self.statusbar.showMessage("Simulation running...")

        self._draw_deflectors()
        self._plot_timer.start()

        self.worker.start()

    def on_stop(self):
        """Ask the worker to stop."""
        if self.worker:
            self.worker.stop_flag = True
            self.statusbar.showMessage("Stopping...")

    def on_reset(self):
        """Clear the plots and the console."""
        self.plot_widget.canvas.clear_plots()
        if self._3d_window is not None:
            self._3d_window.clear_3d()
        self.console.clear()
        self.progress_bar.setValue(0)
        self.statusbar.showMessage("Ready")

    def on_unzoom(self):
        """Reset the zoom on all plots."""
        self.plot_widget.toolbar.home()
        self.statusbar.showMessage("Zoom reset")

    def _on_data_file_loaded(self, path):
        """Show the loaded data file and redraw the deflectors."""
        name = os.path.basename(path)
        self.statusbar.showMessage(f"Data file loaded: {name}")
        self._draw_deflectors()

    def _draw_deflectors(self):
        """Draw the deflectors on the XY plot and the 3D view."""
        defl = self.input_panel.get_deflector_geometry()
        self.plot_widget.canvas.draw_deflectors(defl)
        if self._3d_window is not None:
            self._3d_window.draw_deflectors_3d(defl)

    def on_view_3d(self):
        """Open the 3D window or bring it to the front."""
        if self._3d_window is None:
            self._3d_window = Plot3DWindow(self.plot_widget.canvas, parent=self)
        self._3d_window.show()
        self._3d_window.raise_()
        self._3d_window.activateWindow()
        self._3d_window.refresh_3d()
        defl = self.input_panel.get_deflector_geometry()
        self._3d_window.draw_deflectors_3d(defl)

    def on_export_plots(self):
        """Export the plots to an image file."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Export plots", "lions_plots.png",
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if path:
            self.plot_widget.canvas.fig.savefig(path, dpi=150, bbox_inches='tight')
            self.statusbar.showMessage(f"Plots exported: {path}")

    def on_view_bidims(self):
        """Open the bidim viewer on the .bid files next to the data file."""
        import glob
        ficin = self.input_panel.ficin_edit.text()
        search_dir = os.path.dirname(os.path.abspath(ficin)) if ficin else '.'
        bid_files = sorted(glob.glob(os.path.join(search_dir, '*.bid')))
        dlg = BidimViewerDialog(self, paths=bid_files if bid_files else None)
        dlg.show()

    def on_open_arcade(self):
        try:
            from lions.game.arcade import Arcade
        except ImportError as exc:
            QMessageBox.warning(self, "Arcade", f"Needs the game extra: uv run --extra game lions-gui ({exc})")
            return
        self.arcade = Arcade()
        self.arcade.show()

    def on_view_fields(self):
        """Open the field viewer on the map files named in the form."""
        ficin = self.input_panel.ficin_edit.text()
        work_dir = os.path.dirname(os.path.abspath(ficin)) if ficin else '.'
        swap = self.input_panel.repswap_combo.currentText()

        def path(edit):
            name = edit.text().strip()
            return os.path.join(work_dir, name) if name and name != 'no' else ''

        panel = self.input_panel
        ficb, ficv1, ficv2 = path(panel.ficb_edit), path(panel.ficv1_edit), path(panel.ficv2_edit)
        dlg = FieldViewerDialog(self, ficb=ficb, ficv1=ficv1, ficv2=ficv2, swap=swap)
        dlg.show()

    def on_progress(self, current, total):
        if total > 0:
            self.progress_bar.setValue(int(100 * current / total))

    def on_iteration_data(self, data):
        self._pending_data.append(data)

    def on_console_output(self, text):
        self._pending_console.append(text)

    def _flush_all(self):
        """Flush buffered plot data and console text in one pass (timer slot)."""
        if self._pending_data:
            canvas = self.plot_widget.canvas
            for d in self._pending_data:
                canvas.append_data(d)
            self._pending_data.clear()
            canvas._refresh_all_curves()
            if self._3d_window is not None and self._3d_window.isVisible():
                self._3d_window.refresh_3d()
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
        self.statusbar.showMessage("Simulation finished")

    def on_error(self, msg):
        self._plot_timer.stop()
        self._flush_all()
        self.action_run.setEnabled(True)
        self.action_stop.setEnabled(False)
        self.statusbar.showMessage("Error")
        QMessageBox.critical(self, "Simulation error", msg)

    def on_input_needed(self, prompt):
        """Ask the user for an answer the input sequence did not provide."""
        text, ok = QInputDialog.getText(self, "Input required", prompt)
        if ok and self.worker:
            self.worker.provide_input(text)
        elif self.worker:
            self.worker.provide_input("0")

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

