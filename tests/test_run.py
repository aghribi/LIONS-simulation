"""lions.run on a uniform magnetic field with zero electric fields (1 MeV
proton).
"""

import contextlib
import io
import os

import numpy as np
import pytest

import lions
import synthetic as syn
from lions.apertures import check_box
from lions.fields import compute_bfield, compute_fields

W0 = 1.0  # MeV/u


def uniform_case(bmap_bz, b_orbit, tail=(), turns=10, npaspt=360):
    syn.write_bmap('bmap.dat', bmap_bz)
    syn.write_deck(
        'deck.dat',
        fhf=syn.revolution_frequency(W0, b_orbit),
        r0=syn.larmor_radius(W0, b_orbit),
        w0=W0,
        turns=turns,
        npaspt=npaspt,
        tail=tail,
    )
    return syn.quiet_run('deck.dat')


@pytest.mark.parametrize('b', [-1.0, -0.5])
def test_closed_orbit(workdir, b):
    d = uniform_case(b, b)
    rho = syn.larmor_radius(W0, b)
    assert np.hypot(d.bunch.x[0], d.bunch.y[0]) == pytest.approx(rho, rel=1e-4)
    assert d.bunch.gamma[0] == pytest.approx(
        syn.proton_gamma(W0)[0], rel=1e-12
    )
    assert d.bunch.bz[0] == pytest.approx(b, rel=1e-12)


def test_boris_phase_lag(workdir):
    """Boris rotates by 2 atan(w dt / 2) per step instead of w dt."""
    d = uniform_case(-1.0, -1.0)
    wdt = 2 * np.pi / 360
    predicted = 3600 * (wdt - 2 * np.arctan(wdt / 2))
    assert -np.arctan2(d.bunch.y[0], d.bunch.x[0]) == pytest.approx(
        predicted, rel=1e-2
    )
    assert (
        d.run.turn == 9
    )  # the lag leaves the particle just before its 10th crossing


def test_map_used_as_read_without_b_factor(workdir):
    d = uniform_case(-0.5, -0.5)
    assert d.field.b_factor is None
    assert d.bunch.bz[0] == pytest.approx(-0.5, rel=1e-12)


def test_b_factor_rescales_map(workdir):
    d = uniform_case(-0.5, -1.0, tail=['1.0'])
    assert d.field.b_factor == 1.0
    assert d.bunch.bz[0] == pytest.approx(-1.0, rel=1e-12)
    assert np.hypot(d.bunch.x[0], d.bunch.y[0]) == pytest.approx(
        syn.larmor_radius(W0, -1.0), rel=1e-4
    )


@pytest.mark.parametrize(
    'tail, expected',
    [
        ([], (None, 'no', 'no')),
        (['1.0'], (1.0, 'no', 'no')),
        (['no saved.bin'], (None, 'no', 'saved.bin')),
        (['1.0', 'no saved.bin'], (1.0, 'no', 'saved.bin')),
        (['no no', 'interpolation=hermite'], (None, 'no', 'no')),
    ],
)
def test_deck_tail(workdir, tail, expected):
    d = uniform_case(-1.0, -1.0, tail=tail, turns=1)
    assert (
        d.field.b_factor,
        d.files.bunch_load,
        d.files.bunch_save,
    ) == expected
    assert d.field.interpolation == (
        'hermite' if 'interpolation=hermite' in tail else 'bilinear'
    )
    if expected[2] != 'no':
        assert os.path.getsize(expected[2]) == 6 * 8


def test_on_step_stops_and_closes_files(workdir):
    syn.write_bmap('bmap.dat', -1.0)
    syn.write_deck(
        'deck.dat',
        fhf=syn.revolution_frequency(W0, -1.0),
        r0=syn.larmor_radius(W0, -1.0),
        w0=W0,
    )
    d = syn.quiet_run('deck.dat', on_step=lambda d: d.run.step < 100)
    assert d.run.step == 100
    assert not d.outputs.open


@pytest.mark.parametrize('ntet', [361, 721, 1081])
def test_field_just_below_the_x_axis(workdir, ntet):
    """theta = 2 pi - arccos(1) = 2 pi for a tiny negative y must read map
    column 0, not ntet.
    """
    syn.write_bmap('bmap.dat', lambda r, t: -1.0 - 0.1 * np.sin(t), ntet=ntet)
    syn.write_deck(
        'deck.dat',
        fhf=syn.revolution_frequency(W0, -1.0),
        r0=syn.larmor_radius(W0, -1.0),
        w0=W0,
        turns=1,
    )
    d = lions.State()
    d.run.swap_bytes = 'n'
    with contextlib.redirect_stdout(io.StringIO()):
        lions.load(d, 'deck.dat')
    d.bunch.y[:] = -1e-300
    check_box(d.bunch, d.field, d.box)
    compute_fields(d)
    assert d.bunch.bz[0] == pytest.approx(-1.0, rel=1e-12)
    d.bunch.bz[:] = 0.0
    compute_bfield(d.bunch, d.field)
    assert d.bunch.bz[0] == pytest.approx(-1.0, rel=1e-12)


HERMITE = ['no no', 'interpolation=hermite']


@pytest.mark.parametrize('b', [-1.0, -0.5])
def test_hermite_closed_orbit(workdir, b):
    """Hermite is exact on a uniform map: same orbit as the bilinear model."""
    d = uniform_case(b, b, tail=HERMITE)
    assert d.field.interpolation == 'hermite'
    assert np.hypot(d.bunch.x[0], d.bunch.y[0]) == pytest.approx(
        syn.larmor_radius(W0, b), rel=1e-4
    )
    assert d.bunch.bz[0] == pytest.approx(b, rel=1e-12)


def smooth_field(r, t):
    """Analytic Bz with azimuthal and radial variation, its dBz/dx and
    dBz/dy.
    """
    bz = -1.0 * (1 + 0.05 * np.cos(4 * t)) - 0.4 * (r - 0.25) ** 2
    d_r, d_t = -0.8 * (r - 0.25), 0.2 * np.sin(4 * t)
    return (
        bz,
        np.cos(t) * d_r - np.sin(t) / r * d_t,
        np.sin(t) * d_r + np.cos(t) / r * d_t,
    )


def field_at(tail, r, t, z):
    """Bz and (Bx, By) / z from compute_fields at the points (r, t, z) on
    the smooth_field map.
    """
    syn.write_bmap(
        'bmap.dat',
        bz=lambda r_, t_: smooth_field(r_, t_)[0],
        dbzx=lambda r_, t_: smooth_field(r_, t_)[1],
        dbzy=lambda r_, t_: smooth_field(r_, t_)[2],
    )
    syn.write_deck(
        'deck.dat',
        fhf=syn.revolution_frequency(W0, -1.0),
        r0=0.25,
        w0=W0,
        turns=1,
        nbrion=len(r),
        tail=tail,
    )
    d = lions.State()
    d.run.swap_bytes = 'n'
    with contextlib.redirect_stdout(io.StringIO()):
        lions.load(d, 'deck.dat')
    b = d.bunch
    b.x[:], b.y[:], b.z[:] = r * np.cos(t), r * np.sin(t), z
    check_box(b, d.field, d.box)
    compute_fields(d)
    return b.bz.copy(), b.bx / z, b.by / z


def test_hermite_nodes_and_accuracy(workdir):
    """At the nodes the option returns the map; between them it is closer
    to the analytic field.
    """
    nodes_r = np.linspace(0.05, 0.5, 46)[5:40:7]
    nodes_t = np.linspace(0, 2 * np.pi, 361)[3:350:37]
    r, t = (a.ravel() for a in np.meshgrid(nodes_r, nodes_t))
    bz, gx, gy = field_at(HERMITE, r, t, 1e-4)
    exact = smooth_field(r, t)
    np.testing.assert_allclose(bz, exact[0], rtol=1e-12)
    np.testing.assert_allclose(gx, exact[1], rtol=0, atol=1e-10)
    np.testing.assert_allclose(gy, exact[2], rtol=0, atol=1e-10)

    rm, tm = r + 0.37 * 0.01, t + 0.61 * np.radians(1)
    exact = smooth_field(rm, tm)
    err = {
        name: [
            np.max(np.abs(a - e))
            for a, e in zip(field_at(tail, rm, tm, 1e-4), exact)
        ]
        for name, tail in (('bilinear', []), ('hermite', HERMITE))
    }
    for k in range(3):
        assert err['hermite'][k] < err['bilinear'][k] / 10, (k, err)


def test_atan2_azimuth_near_the_x_axis(workdir):
    """At theta = 1e-9 rad atan2 gives the azimuth to rounding; arccos(x /
    r) loses it (error ~ eps / sin theta).
    """
    theta = 1e-9 * np.arange(1, 6)
    r = np.full(5, 0.25)
    err = {}
    for name, tail in (('atan2', ['no no', 'azimuth=atan2']), ('arccos', [])):
        syn.write_bmap('bmap.dat', -1.0)
        syn.write_deck(
            'deck.dat',
            fhf=syn.revolution_frequency(W0, -1.0),
            r0=0.25,
            w0=W0,
            turns=1,
            nbrion=5,
            tail=tail,
        )
        d = lions.State()
        d.run.swap_bytes = 'n'
        with contextlib.redirect_stdout(io.StringIO()):
            lions.load(d, 'deck.dat')
        b = d.bunch
        b.x[:], b.y[:], b.z[:] = r * np.cos(theta), r * np.sin(theta), 0.0
        check_box(b, d.field, d.box)
        compute_fields(d)
        err[name] = np.max(np.abs(b.theta - theta))
    assert err['atan2'] < 1e-24 and err['arccos'] > 1e-10, err
