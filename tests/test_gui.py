"""Headless GUI: main window, deck round trip through the form, worker
run."""

import os
from pathlib import Path

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
QtWidgets = pytest.importorskip('PyQt5.QtWidgets')

import synthetic as syn  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
EXAMPLE = REPO / 'examples' / 'don.120.ref_phi10_mars_2017'


@pytest.fixture(scope='module')
def app():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_main_window_builds(app):
    from lions.gui.app import LionsMainWindow

    window = LionsMainWindow()
    assert window.input_panel is not None
    window.close()


def test_deck_round_trip(app, tmp_path):
    from lions.gui.input_form import InputFormPanel

    first = InputFormPanel()
    first.load_data_file(str(EXAMPLE))
    first.write_data_file(str(tmp_path / 'deck.dat'))
    second = InputFormPanel()
    second.load_data_file(str(tmp_path / 'deck.dat'))
    assert normalise(second._collect_all_params()) == normalise(
        first._collect_all_params()
    )


def normalise(params):
    """Free-text fields compared by value: '6990000' and '6.99d+06' are
    the same number.
    """

    def value(v):
        if isinstance(v, str):
            try:
                return tuple(
                    float(t.replace('d', 'e').replace('D', 'e'))
                    for t in v.split()
                )
            except ValueError:
                return v
        if isinstance(v, (list, tuple)):
            return tuple(value(x) for x in v)
        if isinstance(v, dict):
            return {k: value(x) for k, x in v.items()}
        return v

    return value(params)


def test_worker_runs(app, workdir):
    from lions.gui.worker import SimulationWorker

    syn.write_bmap('bmap.dat', -1.0)
    syn.write_deck(
        'deck.dat',
        fhf=syn.revolution_frequency(1.0, -1.0),
        r0=syn.larmor_radius(1.0, -1.0),
        w0=1.0,
        turns=1,
    )
    answers = [
        'n',
        str(workdir / 'deck.dat'),
        '1000 2000 0 0',
        '0 0 0 0 0 0',
        '0',
    ]
    worker = SimulationWorker(answers, str(workdir))
    events = {'progress': [], 'data': 0, 'done': False, 'errors': []}
    worker.progress.connect(lambda i, n: events['progress'].append(i))
    worker.iteration_data.connect(
        lambda d: events.__setitem__('data', events['data'] + 1)
    )
    worker.simulation_finished.connect(
        lambda: events.__setitem__('done', True)
    )
    worker.simulation_error.connect(events['errors'].append)
    worker.run()
    assert events['errors'] == []
    assert events['done']
    assert events['progress'][0] == 0 and events['progress'][-1] == 360
    assert events['data'] == 3  # step 0, step 1, step 360


@pytest.mark.parametrize(
    'ntheta, theta_max',
    [(361, 2 * 3.141592653589793), (91, 3.141592653589793 / 4)],
)
def test_field_viewer_plots_engine_readers(app, workdir, ntheta, theta_max):
    from lions.gui.field_viewer import FieldViewerDialog
    from lions.io import (
        read_bfield_map,
        read_central_potential,
        read_outer_potential,
    )

    syn.write_bmap('bmap.dat', -1.0, ntet=ntheta, theta_max=theta_max)
    dialog = FieldViewerDialog()
    dialog.mag_data = read_bfield_map('bmap.dat')
    dialog.pot1_data = read_central_potential('pot1.bin')
    dialog.pot2_data = read_outer_potential('pot2.dat')
    for index in range(dialog.field_combo.count()):
        dialog.field_combo.setCurrentIndex(index)
        dialog._replot()
    dialog.close()


def test_field_viewer_opens_on_deck_maps(app):
    from lions.gui.app import LionsMainWindow
    from lions.gui.field_viewer import FieldViewerDialog

    window = LionsMainWindow()
    window.input_panel.load_data_file(str(EXAMPLE))
    window.input_panel.ficin_edit.setText(str(EXAMPLE))
    window.on_view_fields()
    dialog = window.findChildren(FieldViewerDialog)[-1]
    names = [
        Path(e.text()).name
        for e in (dialog.ficb_edit, dialog.ficv1_edit, dialog.ficv2_edit)
    ]
    assert names == [
        'cime700rec.1804.pbob1.lio.dat',
        'pot3d.impair.30pilbleu30.620.dat',
        'newgap612.dat',
    ]
    dialog.close()
    window.close()


def test_form_keeps_every_deck_entry(app, tmp_path):
    """Load and save through the form: gap offset, bunch files and all
    sections survive.
    """
    from lions.deck import read_deck, write_deck
    from lions.gui.input_form import InputFormPanel
    from test_deck import full_deck

    deck = full_deck()
    write_deck(deck, tmp_path / 'deck.dat')
    form = InputFormPanel()
    form.load_data_file(str(tmp_path / 'deck.dat'))
    form.write_data_file(str(tmp_path / 'saved.dat'))
    assert read_deck(tmp_path / 'saved.dat') == deck


def test_form_keeps_deck_options(app, tmp_path):
    """interpolation shows in the form; azimuth, without a widget, is
    kept from the loaded deck.
    """
    from lions.deck import read_deck, write_deck
    from lions.gui.input_form import InputFormPanel
    from test_deck import full_deck

    deck = full_deck()
    deck.field_interpolation, deck.field_azimuth = 'hermite', 'atan2'
    write_deck(deck, tmp_path / 'deck.dat')
    form = InputFormPanel()
    form.load_data_file(str(tmp_path / 'deck.dat'))
    assert form.field_interp_combo.currentText() == 'hermite'
    form.write_data_file(str(tmp_path / 'saved.dat'))
    assert read_deck(tmp_path / 'saved.dat') == deck
    form.field_interp_combo.setCurrentText('bilinear')
    assert form.to_deck().field_interpolation == 'bilinear'
