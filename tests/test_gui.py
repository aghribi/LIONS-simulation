"""Headless GUI: main window, deck round trip through the form, worker run."""
import os
from pathlib import Path

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
QtWidgets = pytest.importorskip('PyQt5.QtWidgets')

import synthetic as syn  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
EXAMPLE = REPO / 'don.120.ref_phi10_mars_2017'


@pytest.fixture(scope='module')
def app():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_main_window_builds(app):
    from lions_gui import LionsMainWindow
    window = LionsMainWindow()
    assert window.input_panel is not None
    window.close()


def test_deck_round_trip(app, tmp_path):
    from lions_input_form import InputFormPanel
    first = InputFormPanel()
    first.load_data_file(str(EXAMPLE))
    first.write_data_file(str(tmp_path / 'deck.dat'))
    second = InputFormPanel()
    second.load_data_file(str(tmp_path / 'deck.dat'))
    assert normalise(second._collect_all_params()) == normalise(first._collect_all_params())


def normalise(params):
    """Free-text fields compared by value: '6990000' and '6.99d+06' are the same number."""
    def value(v):
        if isinstance(v, str):
            try:
                return tuple(float(t.replace('d', 'e').replace('D', 'e')) for t in v.split())
            except ValueError:
                return v
        if isinstance(v, (list, tuple)):
            return tuple(value(x) for x in v)
        if isinstance(v, dict):
            return {k: value(x) for k, x in v.items()}
        return v
    return value(params)


def test_worker_runs(app, workdir):
    from lions_worker import SimulationWorker
    syn.write_bmap('bmap.dat', -1.0)
    syn.write_deck('deck.dat', fhf=syn.revolution_frequency(1.0, -1.0),
                   r0=syn.larmor_radius(1.0, -1.0), w0=1.0, turns=1)
    answers = ['n', str(workdir / 'deck.dat'), '1000 2000 0 0', '0 0 0 0 0 0', '0']
    worker = SimulationWorker(answers, str(workdir))
    events = {'progress': [], 'data': 0, 'done': False, 'errors': []}
    worker.progress.connect(lambda i, n: events['progress'].append(i))
    worker.iteration_data.connect(lambda d: events.__setitem__('data', events['data'] + 1))
    worker.simulation_finished.connect(lambda: events.__setitem__('done', True))
    worker.simulation_error.connect(events['errors'].append)
    worker.run()
    assert events['errors'] == []
    assert events['done']
    assert events['progress'][0] == 0 and events['progress'][-1] == 360
    assert events['data'] == 3  # step 0, step 1, step 360
