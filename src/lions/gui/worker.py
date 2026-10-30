"""
lions_worker.py - Worker QThread pour la simulation LIONS
Intercepte input() et print() pour les router vers la GUI.
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
    """Remplace builtins.input : fournit les reponses pre-chargees."""

    def __init__(self, answers):
        self.q = queue.Queue()
        for a in answers:
            self.q.put(a)
        self.fallback_signal = None  # pyqtSignal(str) pour inputs inattendus
        self.response_queue = queue.Queue()

    def __call__(self, prompt=''):
        try:
            answer = self.q.get_nowait()
            return answer
        except queue.Empty:
            if self.fallback_signal:
                self.fallback_signal.emit(prompt or "Entree requise :")
                try:
                    return self.response_queue.get(timeout=300)
                except queue.Empty:
                    raise RuntimeError("Timeout en attente de reponse utilisateur")
            raise RuntimeError(f"Input inattendu: {prompt}")


class OutputInterceptor(io.TextIOBase):
    """Intercepte sys.stdout pour emettre un signal Qt."""

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
    """Thread qui execute la simulation LIONS avec hooks pour la GUI."""

    progress = pyqtSignal(int, int)          # (iteration, total)
    iteration_data = pyqtSignal(dict)        # donnees pour les plots
    console_output = pyqtSignal(str)         # texte de sortie
    simulation_finished = pyqtSignal()       # fin normale
    simulation_error = pyqtSignal(str)       # erreur
    input_needed = pyqtSignal(str)           # input inattendu

    def __init__(self, input_answers, work_dir=None):
        super().__init__()
        self.input_answers = input_answers
        self.work_dir = work_dir
        self.stop_flag = False
        self._input_queue = None

    def provide_input(self, answer):
        """Appele par la GUI quand un input inattendu est necessaire."""
        if self._input_queue:
            self._input_queue.response_queue.put(answer)

    def run(self):
        # Importer lions ici (dans le thread) pour eviter les problemes
        # de modules charges au mauvais moment
        original_input = builtins.input
        original_stdout = sys.stdout
        original_dir = os.getcwd()

        try:
            # Changer de repertoire vers le fichier de donnees
            if self.work_dir:
                os.chdir(self.work_dir)

            # Installer les intercepteurs I/O
            self._input_queue = InputQueue(self.input_answers)
            self._input_queue.fallback_signal = self.input_needed
            interceptor = OutputInterceptor(self.console_output)
            builtins.input = self._input_queue
            sys.stdout = interceptor

            # Import du moteur de simulation

            # Reponses fournies par la file d'entree (meme sequence que le CLI)
            cli.banner()
            d = State()
            d.swap_bytes = input().strip()
            load(d, input().strip())
            cli.read_vertical_sorter(d)
            bidim = BidimWriter()
            setup_bidims(d, bidim)
            cli.track(d, bidim, on_step=self._on_step)
            cli.write_trapped(d)

            if not self.stop_flag:
                print('\n=== Simulation terminee ===')

            self.simulation_finished.emit()

        except (Exception, SystemExit) as e:
            err_msg = f"{str(e)}\n{traceback.format_exc()}"
            self.simulation_error.emit(err_msg)

        finally:
            builtins.input = original_input
            sys.stdout = original_stdout
            os.chdir(original_dir)

    def _on_step(self, d):
        """Hook de lions.cli.track : progression, donnees de plot, arret."""
        ite = d.step
        # Progression limitee a ~200 emissions pour ne pas noyer la boucle Qt
        if ite % max(1, d.n_steps // 200) == 0 or ite == 1 or ite == d.n_steps:
            self.progress.emit(ite, d.n_steps)
        if ite == 0 or (d.output_every > 0 and
                        (ite % d.output_every == 0 or ite == 1 or ite == d.n_steps)):
            self._emit_data(d)
        if self.stop_flag:
            print(f'\n--- Simulation arretee a l iteration {ite + 1} ---')
            return False
        return True

    def _emit_data(self, d):
        """Extrait les donnees de la simulation pour les plots (par famille)."""
        try:
            nbrfam = d.n_families
            for k in range(nbrfam):
                i1 = d.family_start[k]
                rr = np.sqrt(d.x[i1]**2 + d.y[i1]**2)
                wmevpa = d.rest_energy[i1] * (d.gamma[i1] - 1.0) / MEV / d.mass_number[k]
                ppp = np.sqrt(d.px[i1]**2 + d.py[i1]**2 + d.pz[i1]**2)
                tpz = 0.0
                if ppp > 0:
                    tpz = np.arcsin(np.clip(d.pz[i1] / ppp, -1.0, 1.0)) * 1000.0

                # Calculer theta pour cette famille (comme info() le fait pour fam 0)
                if rr > 0:
                    tetnew_k = np.arccos(np.clip(d.x[i1] / rr, -1.0, 1.0))
                    if d.y[i1] < 0.0:
                        tetnew_k = 2.0 * np.pi - tetnew_k
                    tetnew_k = 180.0 / np.pi * tetnew_k
                else:
                    tetnew_k = 0.0

                # Calculer dphi (dephasage HF)
                phihf = 180.0 / np.pi * d.omega_rf * d.t
                tetall_k = tetnew_k + 360.0 * (d.turn - 1)
                dph = phihf - d.harmonic * tetall_k - 90.0

                data = {
                    'family': k,
                    'nbrfam': nbrfam,
                    'ite': d.step,
                    'numtou': d.turn,
                    'x0_mm': float(d.x[i1] * 1000.0),
                    'y0_mm': float(d.y[i1] * 1000.0),
                    'r0_mm': float(rr * 1000.0),
                    'tetnew': float(tetnew_k),
                    'z0_mm': float(d.z[i1] * 1000.0),
                    'tpz': float(tpz),
                    'bz0': float(d.bz[i1]),
                    'wmevpa': float(wmevpa),
                    'dph': float(dph),
                }

                # === Scatter: toutes les particules actives de la famille ===
                n_k = d.family_size[k]
                if n_k > 1:
                    idx_s = d.family_start[k]
                    idx_e = idx_s + n_k
                    ok_k = d.alive[idx_s:idx_e]

                    x_a = d.x[idx_s:idx_e][ok_k] * 1000.0
                    y_a = d.y[idx_s:idx_e][ok_k] * 1000.0
                    z_a = d.z[idx_s:idx_e][ok_k] * 1000.0
                    r_a = np.sqrt(x_a**2 + y_a**2)

                    # Theta (deg)
                    safe_r = np.where(r_a > 0, r_a, 1.0)
                    theta_a = np.arccos(np.clip(x_a / safe_r, -1.0, 1.0))
                    theta_a[y_a < 0.0] = 2.0 * np.pi - theta_a[y_a < 0.0]
                    theta_a = 180.0 / np.pi * theta_a
                    theta_a[r_a == 0] = 0.0

                    # Pz angle (mrad)
                    px_a = d.px[idx_s:idx_e][ok_k]
                    py_a = d.py[idx_s:idx_e][ok_k]
                    pz_a = d.pz[idx_s:idx_e][ok_k]
                    ppp_a = np.sqrt(px_a**2 + py_a**2 + pz_a**2)
                    safe_p = np.where(ppp_a > 0, ppp_a, 1.0)
                    tpz_a = np.arcsin(np.clip(pz_a / safe_p, -1.0, 1.0)) * 1000.0
                    tpz_a[ppp_a == 0] = 0.0

                    # Energie (MeV/A)
                    ga_a = d.gamma[idx_s:idx_e][ok_k]
                    am0c2_a = d.rest_energy[idx_s:idx_e][ok_k]
                    wmevpa_a = am0c2_a * (ga_a - 1.0) / MEV / d.mass_number[k]

                    # dphi (deg)
                    tetall_a = theta_a + 360.0 * (d.turn - 1)
                    dph_a = phihf - d.harmonic * tetall_a - 90.0

                    data['scatter_x'] = x_a
                    data['scatter_y'] = y_a
                    data['scatter_z'] = z_a
                    data['scatter_r'] = r_a
                    data['scatter_theta'] = theta_a
                    data['scatter_pz'] = tpz_a
                    data['scatter_dph'] = dph_a
                    data['scatter_wmevpa'] = wmevpa_a
                    data['scatter_turn'] = d.turn - 1.0 + theta_a / 360.0

                self.iteration_data.emit(data)
        except Exception:
            pass  # Ne pas bloquer la simulation pour un probleme de plot
