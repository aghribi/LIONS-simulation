"""QThread that runs the LIONS engine and routes input() and print() to the
GUI.
"""

import sys
import io
import os
import queue
import builtins
import traceback
import numpy as np
from PyQt5.QtCore import QThread, pyqtSignal

from .. import cli
from ..bidim import BidimWriter, setup_bidims
from ..constants import MEV
from ..io import load
from ..state import State


class InputQueue:
    """Stands in for builtins.input and returns the preloaded answers in
    order.
    """

    def __init__(self, answers):
        self.q = queue.Queue()
        for a in answers:
            self.q.put(a)
        self.fallback_signal = (
            None  # pyqtSignal(str) emitted when the answers run out
        )
        self.response_queue = queue.Queue()

    def __call__(self, prompt=''):
        try:
            answer = self.q.get_nowait()
            return answer
        except queue.Empty:
            if self.fallback_signal:
                self.fallback_signal.emit(prompt or "Input required:")
                try:
                    return self.response_queue.get(timeout=300)
                except queue.Empty:
                    raise RuntimeError("Timed out waiting for user input")
            raise RuntimeError(f"Unexpected input: {prompt}")


class OutputInterceptor(io.TextIOBase):
    """Replaces sys.stdout: emits each write as a Qt signal and echoes it
    to the terminal.
    """

    def __init__(self, signal):
        super().__init__()
        self.signal = signal
        self._original = sys.__stdout__

    def write(self, text):
        if text:
            self.signal.emit(text)
            if self._original:
                self._original.write(text)
        return len(text) if text else 0

    def flush(self):
        if self._original:
            self._original.flush()


class SimulationWorker(QThread):
    """Runs the LIONS simulation with GUI hooks."""

    progress = pyqtSignal(int, int)  # (iteration, total)
    iteration_data = pyqtSignal(dict)  # plot data
    console_output = pyqtSignal(str)  # console text
    simulation_finished = pyqtSignal()  # normal end
    simulation_error = pyqtSignal(str)  # error message
    input_needed = pyqtSignal(str)  # unexpected input() prompt

    def __init__(self, input_answers, work_dir=None):
        super().__init__()
        self.input_answers = input_answers
        self.work_dir = work_dir
        self.stop_flag = False
        self._input_queue = None

    def provide_input(self, answer):
        """Answer an unexpected input() call from the GUI."""
        if self._input_queue:
            self._input_queue.response_queue.put(answer)

    def run(self):
        original_input = builtins.input
        original_stdout = sys.stdout
        original_dir = os.getcwd()

        try:
            # Paths in the data file are relative to its directory.
            if self.work_dir:
                os.chdir(self.work_dir)

            self._input_queue = InputQueue(self.input_answers)
            self._input_queue.fallback_signal = self.input_needed
            interceptor = OutputInterceptor(self.console_output)
            builtins.input = self._input_queue
            sys.stdout = interceptor

            # Same input() sequence as lions.cli.main.
            cli.banner()
            d = State()
            d.run.swap_bytes = input().strip()
            load(d, input().strip())
            cli.read_vertical_sorter(d)
            bidim = BidimWriter()
            setup_bidims(d, bidim)
            cli.track(d, bidim, on_step=self._on_step)
            cli.write_trapped(d.bunch, d.traps)

            if not self.stop_flag:
                print('\n=== Simulation finished ===')

            self.simulation_finished.emit()

        except (Exception, SystemExit) as e:
            err_msg = f"{str(e)}\n{traceback.format_exc()}"
            self.simulation_error.emit(err_msg)

        finally:
            builtins.input = original_input
            sys.stdout = original_stdout
            os.chdir(original_dir)

    def _on_step(self, d):
        """lions.cli.track hook: progress, plot data, stop request."""
        ite = d.run.step
        # About 200 progress signals per run, to keep the Qt event loop
        # responsive.
        if (
            ite % max(1, d.run.n_steps // 200) == 0
            or ite == 1
            or ite == d.run.n_steps
        ):
            self.progress.emit(ite, d.run.n_steps)
        if ite == 0 or (
            d.run.output_every > 0
            and (
                ite % d.run.output_every == 0
                or ite == 1
                or ite == d.run.n_steps
            )
        ):
            self._emit_data(d)
        if self.stop_flag:
            print(f'\n--- Simulation stopped at step {ite + 1} ---')
            return False
        return True

    def _emit_data(self, d):
        """Emit one plot-data dict per family: reference particle plus all
        live particles.
        """
        try:
            nbrfam = len(d.bunch.families)
            for k in range(nbrfam):
                i1 = d.bunch.families[k].start
                rr = np.sqrt(d.bunch.x[i1] ** 2 + d.bunch.y[i1] ** 2)
                wmevpa = (
                    d.bunch.rest_energy[i1]
                    * (d.bunch.gamma[i1] - 1.0)
                    / MEV
                    / d.bunch.families[k].mass_number
                )
                ppp = np.sqrt(
                    d.bunch.px[i1] ** 2
                    + d.bunch.py[i1] ** 2
                    + d.bunch.pz[i1] ** 2
                )
                tpz = 0.0
                if ppp > 0:
                    tpz = (
                        np.arcsin(np.clip(d.bunch.pz[i1] / ppp, -1.0, 1.0))
                        * 1000.0
                    )

                # Azimuth of the reference particle, in degrees.
                if rr > 0:
                    tetnew_k = np.arccos(
                        np.clip(d.bunch.x[i1] / rr, -1.0, 1.0)
                    )
                    if d.bunch.y[i1] < 0.0:
                        tetnew_k = 2.0 * np.pi - tetnew_k
                    tetnew_k = 180.0 / np.pi * tetnew_k
                else:
                    tetnew_k = 0.0

                # RF phase slip (deg).
                phihf = 180.0 / np.pi * d.rf.omega * d.run.t
                tetall_k = tetnew_k + 360.0 * (d.run.turn - 1)
                dph = phihf - d.rf.harmonic * tetall_k - 90.0

                data = {
                    'family': k,
                    'nbrfam': nbrfam,
                    'ite': d.run.step,
                    'numtou': d.run.turn,
                    'x0_mm': float(d.bunch.x[i1] * 1000.0),
                    'y0_mm': float(d.bunch.y[i1] * 1000.0),
                    'r0_mm': float(rr * 1000.0),
                    'tetnew': float(tetnew_k),
                    'z0_mm': float(d.bunch.z[i1] * 1000.0),
                    'tpz': float(tpz),
                    'bz0': float(d.bunch.bz[i1]),
                    'wmevpa': float(wmevpa),
                    'dph': float(dph),
                }

                # Scatter: all live particles of the family.
                n_k = d.bunch.families[k].size
                if n_k > 1:
                    idx_s = d.bunch.families[k].start
                    idx_e = idx_s + n_k
                    ok_k = d.bunch.alive[idx_s:idx_e]

                    x_a = d.bunch.x[idx_s:idx_e][ok_k] * 1000.0
                    y_a = d.bunch.y[idx_s:idx_e][ok_k] * 1000.0
                    z_a = d.bunch.z[idx_s:idx_e][ok_k] * 1000.0
                    r_a = np.sqrt(x_a**2 + y_a**2)

                    # Theta (deg)
                    safe_r = np.where(r_a > 0, r_a, 1.0)
                    theta_a = np.arccos(np.clip(x_a / safe_r, -1.0, 1.0))
                    theta_a[y_a < 0.0] = 2.0 * np.pi - theta_a[y_a < 0.0]
                    theta_a = 180.0 / np.pi * theta_a
                    theta_a[r_a == 0] = 0.0

                    # Pz angle (mrad)
                    px_a = d.bunch.px[idx_s:idx_e][ok_k]
                    py_a = d.bunch.py[idx_s:idx_e][ok_k]
                    pz_a = d.bunch.pz[idx_s:idx_e][ok_k]
                    ppp_a = np.sqrt(px_a**2 + py_a**2 + pz_a**2)
                    safe_p = np.where(ppp_a > 0, ppp_a, 1.0)
                    tpz_a = (
                        np.arcsin(np.clip(pz_a / safe_p, -1.0, 1.0)) * 1000.0
                    )
                    tpz_a[ppp_a == 0] = 0.0

                    # Energy (MeV/A)
                    ga_a = d.bunch.gamma[idx_s:idx_e][ok_k]
                    am0c2_a = d.bunch.rest_energy[idx_s:idx_e][ok_k]
                    wmevpa_a = (
                        am0c2_a
                        * (ga_a - 1.0)
                        / MEV
                        / d.bunch.families[k].mass_number
                    )

                    # dphi (deg)
                    tetall_a = theta_a + 360.0 * (d.run.turn - 1)
                    dph_a = phihf - d.rf.harmonic * tetall_a - 90.0

                    data['scatter_x'] = x_a
                    data['scatter_y'] = y_a
                    data['scatter_z'] = z_a
                    data['scatter_r'] = r_a
                    data['scatter_theta'] = theta_a
                    data['scatter_pz'] = tpz_a
                    data['scatter_dph'] = dph_a
                    data['scatter_wmevpa'] = wmevpa_a
                    data['scatter_turn'] = d.run.turn - 1.0 + theta_a / 360.0

                self.iteration_data.emit(data)
        except Exception:
            pass  # a plotting error must not stop the run
