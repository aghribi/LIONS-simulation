"""lions.run on a uniform magnetic field with zero electric fields (1 MeV proton)."""
import os

import numpy as np
import pytest

import synthetic as syn

W0 = 1.0  # MeV/u


def uniform_case(bmap_bz, b_orbit, tail=(), turns=10, npaspt=360):
    syn.write_bmap('bmap.dat', bmap_bz)
    syn.write_deck('deck.dat', fhf=syn.revolution_frequency(W0, b_orbit),
                   r0=syn.larmor_radius(W0, b_orbit), w0=W0,
                   turns=turns, npaspt=npaspt, tail=tail)
    return syn.quiet_run('deck.dat')


@pytest.mark.parametrize('b', [-1.0, -0.5])
def test_closed_orbit(workdir, b):
    d = uniform_case(b, b)
    rho = syn.larmor_radius(W0, b)
    assert np.hypot(d.x[0], d.y[0]) == pytest.approx(rho, rel=1e-4)
    assert d.gamma[0] == pytest.approx(syn.proton_gamma(W0)[0], rel=1e-12)
    assert d.bz[0] == pytest.approx(b, rel=1e-12)


def test_boris_phase_lag(workdir):
    """Boris rotates by 2 atan(w dt / 2) per step instead of w dt."""
    d = uniform_case(-1.0, -1.0)
    wdt = 2 * np.pi / 360
    predicted = 3600 * (wdt - 2 * np.arctan(wdt / 2))
    assert -np.arctan2(d.y[0], d.x[0]) == pytest.approx(predicted, rel=1e-2)
    assert d.turn == 9  # the lag leaves the particle just before its 10th crossing


def test_map_used_as_read_without_b_factor(workdir):
    d = uniform_case(-0.5, -0.5)
    assert d.b_factor is None
    assert d.bz[0] == pytest.approx(-0.5, rel=1e-12)


def test_b_factor_rescales_map(workdir):
    d = uniform_case(-0.5, -1.0, tail=['1.0'])
    assert d.b_factor == 1.0
    assert d.bz[0] == pytest.approx(-1.0, rel=1e-12)
    assert np.hypot(d.x[0], d.y[0]) == pytest.approx(syn.larmor_radius(W0, -1.0), rel=1e-4)


@pytest.mark.parametrize('tail, expected', [
    ([], (None, 'no', 'no')),
    (['1.0'], (1.0, 'no', 'no')),
    (['no saved.bin'], (None, 'no', 'saved.bin')),
    (['1.0', 'no saved.bin'], (1.0, 'no', 'saved.bin')),
])
def test_deck_tail(workdir, tail, expected):
    d = uniform_case(-1.0, -1.0, tail=tail, turns=1)
    assert (d.b_factor, d.bunch_load_file, d.bunch_save_file) == expected
    if expected[2] != 'no':
        assert os.path.getsize(expected[2]) == 6 * 8


def test_on_step_stops_and_closes_files(workdir):
    syn.write_bmap('bmap.dat', -1.0)
    syn.write_deck('deck.dat', fhf=syn.revolution_frequency(W0, -1.0),
                   r0=syn.larmor_radius(W0, -1.0), w0=W0)
    d = syn.quiet_run('deck.dat', on_step=lambda d: d.step < 100)
    assert d.step == 100
    assert not d.outputs
