"""Python port against the Fortran reference on synthetic inputs.

Both programs run as executables with the same stdin, each in its own
directory with the same input files. Compared: final bunch (ficsto, full
precision), the numeric rows of every text output, and the case-specific
output files.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import synthetic as syn
from lions.deck import read_deck, write_deck

REPO = Path(__file__).resolve().parent.parent
W0 = 1.0  # MeV/u
B0 = -1.0  # T
V_DEE = 50e3  # V
NOVERT = ['1000 2000 0 0', '0 0 0 0 0 0']
N_BUNCH = 20


@pytest.fixture(scope='session')
def fortran(tmp_path_factory):
    if shutil.which('gfortran') is None:
        pytest.skip('gfortran not found')
    build = tmp_path_factory.mktemp('fortran')
    exe = build / 'lions_f'
    # No FMA contraction: NumPy does not fuse, and fusing alone shifts
    # the last bits.
    subprocess.run(
        [
            'gfortran',
            '-O2',
            '-ffp-contract=off',
            '-J',
            str(build),
            '-o',
            str(exe),
            str(REPO / 'fortran' / 'lion_m.f90'),
            str(REPO / 'fortran' / 'lions_2.0_ep.f90'),
        ],
        check=True,
        capture_output=True,
    )
    return exe


def orbit(**kwargs):
    return dict(
        fhf=syn.revolution_frequency(W0, B0),
        r0=syn.larmor_radius(W0, B0),
        w0=W0,
        **kwargs,
    )


def bunch(r0):
    """Deterministic 20-particle bunch around (r0, 0, 0), moving along
    +y."""
    p0 = np.sqrt(syn.proton_gamma(W0)[0] ** 2 - 1.0)
    u = np.linspace(-1.0, 1.0, N_BUNCH)
    v = np.roll(u, 7)
    pos = np.column_stack([r0 + 1e-3 * u, 1e-3 * v, 5e-4 * np.roll(u, 3)])
    mom = np.column_stack(
        [2e-3 * p0 * v, p0 * (1 + 1e-3 * np.roll(u, 11)), 1e-3 * p0 * u]
    )
    pos[0], mom[0] = (r0, 0.0, 0.0), (0.0, p0, 0.0)
    return pos, mom


def case_outer(bunch_file):
    """4-sector Bz with off-plane terms, vertical offset and angle, dees
    via the outer potential.
    """
    syn.write_bmap(
        'bmap.dat',
        bz=lambda r, t: B0 * (1 + 0.05 * np.cos(4 * t)),
        d2bz=-2.0,
        dbzx=lambda r, t: 0.3 * np.sin(4 * t),
        dbzy=lambda r, t: 0.3 * np.cos(4 * t),
    )
    syn.write_pot_outer(
        'pot2.dat', pot=lambda x, y: 0.5 * (1 + np.tanh(y / 0.01))
    )
    syn.write_deck(
        'deck.dat',
        **orbit(v_dee=V_DEE, z0=0.002, pr0=10.0, pz0=5.0),
        tail=[f'no {bunch_file}'],
    )
    return NOVERT + ['0'], []


def case_eighth_map(bunch_file):
    """Same field as case_outer, given as a 1/8 map (0 to 45 deg, 91
    columns).
    """
    case_outer(bunch_file)
    syn.write_bmap(
        'bmap.dat',
        bz=lambda r, t: B0 * (1 + 0.05 * np.cos(4 * t)),
        d2bz=-2.0,
        dbzx=lambda r, t: 0.3 * np.sin(4 * t),
        dbzy=lambda r, t: 0.3 * np.cos(4 * t),
        ntet=91,
        theta_max=np.pi / 4,
    )
    return NOVERT + ['0'], []


def case_central(bunch_file):
    """Whole orbit in the central region: pot1 and its d2pot (Ez) on a 3
    mm grid.

    The grid must cover r < rmax1: neither code checks the grid bounds
    there.
    """
    syn.write_bmap('bmap.dat', B0)
    syn.write_pot_central(
        'pot1.bin',
        pot=lambda x, y: np.tanh(y / 0.02),
        d2pot=lambda x, y: -100.0 * np.tanh(y / 0.02),
        half_width_mm=300.0,
        nx=201,
        ny=201,
    )
    syn.write_deck(
        'deck.dat',
        **orbit(v_dee=V_DEE, z0=0.002, pz0=5.0),
        rmax1=0.25,
        tail=[f'no {bunch_file}'],
    )
    return NOVERT + ['0'], []


def case_extraction(bunch_file):
    """Electrostatic deflector (30-90 deg), magnetic dipole (120-180
    deg), grids at 200 and 250 deg.
    """
    syn.write_bmap('bmap.dat', B0)
    grids = [
        f'{a} {1000 * np.cos(np.radians(a)):.12g}'
        f' {1000 * np.sin(np.radians(a)):.12g} 0 0'
        for a in (200, 250)
    ]
    extraction = [
        '-1 -1 1 0',
        '1',
        '0.14 2e5 0 0 0 0 30 90',
        '1',
        '0.15 0.5 0 0 120 180',
        '2',
        *grids,
        'grid',
    ]
    syn.write_deck(
        'deck.dat',
        **orbit(turns=3),
        extraction=extraction,
        tail=[f'no {bunch_file}'],
    )
    return NOVERT + ['0'], ['grid.1', 'grid.2']


def write_correlated(bunch_file, loaded):
    syn.write_bmap('bmap.dat', B0)
    with open('sigma.dat', 'w') as f:
        for i in range(6):
            f.write(
                ' '.join('1e-8' if i == j else '0' for j in range(6)) + '\n'
            )
    pos, mom = bunch(syn.larmor_radius(W0, B0))
    syn.write_bunch('bunch_py.bin', 'bunch_f.bin', pos, mom)
    syn.write_deck(
        'deck.dat',
        **orbit(turns=3),
        nbrion=N_BUNCH,
        correlations=['1', 'sigma.dat', '6d', 'cor.out no no'],
        tail=[f'{loaded} {bunch_file}'],
    )


def case_correlations(bunch_file, loaded):
    """20 particles loaded from the same bunch, ficcor written every 30
    steps.
    """
    write_correlated(bunch_file, loaded)
    return NOVERT + ['0'], ['cor.out']


def case_traps(bunch_file, loaded):
    """Same bunch with bidims on and one trap circle next to particle 0
    that catches part of it.
    """
    write_correlated(bunch_file, loaded)
    r0_mm = 1000 * syn.larmor_radius(W0, B0)
    syn.write_markers('markers.dat', [(r0_mm + 0.6, 0.0, 0.5)])
    bidims = [
        '10',
        '1',
        'c',
        '-0.3 0.3 -0.3 0.3',
        '1 0.01',
        '0.01 0.01',
        '0.01 0.01',
        '1',
        '-0.3 0.3 -0.3 0.3',
        'markers.dat',
        '100000',
        '0',
    ]
    return NOVERT + bidims + ['o', 'trapped.out'], ['cor.out', 'trapped.out']


CASES = {
    'outer': case_outer,
    'eighth_map': case_eighth_map,
    'central': case_central,
    'extraction': case_extraction,
    'correlations': case_correlations,
    'traps': case_traps,
}


def run_in(directory, cmd, setup, loaded, rewrite):
    directory.mkdir()
    with pytest.MonkeyPatch.context() as mp:
        mp.chdir(directory)
        syn.write_pot_central('pot1.bin')
        syn.write_pot_outer('pot2.dat')
        args = ('final.bin',) if loaded is None else ('final.bin', loaded)
        answers, extra = setup(*args)
        if rewrite:
            write_deck(read_deck('deck.dat'), 'deck.dat')
    stdin = '\n'.join(['n', 'deck.dat', *answers]) + '\n'
    subprocess.run(
        cmd,
        input=stdin,
        text=True,
        cwd=directory,
        check=True,
        capture_output=True,
    )
    return extra


def assert_rows_match(path_f, path_p):
    rf, rp = syn.numeric_rows(path_f), syn.numeric_rows(path_p)
    assert len(rf) == len(rp), (
        f'{path_p.name}: {len(rf)} Fortran rows, {len(rp)} Python rows'
    )
    for a, b in zip(rf, rp):
        # 1e-6: Fortran e12.7 fields carry 7 significant digits, and the
        # end-of-run energy uses a single-precision MeV constant
        np.testing.assert_allclose(
            b, a, rtol=1e-6, atol=1e-12, err_msg=path_p.name
        )


@pytest.mark.parametrize('rewrite', [False, True], ids=['deck', 'rewritten'])
@pytest.mark.parametrize('name', CASES)
def test_python_matches_fortran(tmp_path, fortran, name, rewrite):
    """rewritten: the deck goes through read_deck and write_deck before
    both runs.
    """
    setup = CASES[name]
    loaded = {'correlations': True, 'traps': True}.get(name)
    extra = run_in(
        tmp_path / 'f',
        [str(fortran)],
        setup,
        'bunch_f.bin' if loaded else None,
        rewrite,
    )
    run_in(
        tmp_path / 'p',
        [sys.executable, '-m', 'lions'],
        setup,
        'bunch_py.bin' if loaded else None,
        rewrite,
    )

    n = N_BUNCH if loaded else 1
    xf, pf = syn.read_fortran_bunch(tmp_path / 'f' / 'final.bin', n)
    xp, pp = np.fromfile(tmp_path / 'p' / 'final.bin', '<f8').reshape(2, n, 3)
    np.testing.assert_allclose(xp, xf, rtol=0, atol=1e-12)
    np.testing.assert_allclose(pp, pf, rtol=1e-12, atol=1e-15)
    for out in ['z.out', 'xy.out', *extra]:
        assert_rows_match(tmp_path / 'f' / out, tmp_path / 'p' / out)
        assert syn.numeric_rows(tmp_path / 'p' / out), (
            f'{out} empty: case not exercised'
        )


def test_fortran_ignores_deck_options(tmp_path, fortran):
    """A deck with interpolation=hermite and azimuth=atan2 gives the
    Fortran the same run as the plain deck.
    """

    def with_option(bunch_file):
        out = case_outer(bunch_file)
        deck = read_deck('deck.dat')
        deck.field_interpolation, deck.field_azimuth = 'hermite', 'atan2'
        write_deck(deck, 'deck.dat')
        return out

    extra = run_in(tmp_path / 'plain', [str(fortran)], case_outer, None, False)
    run_in(tmp_path / 'option', [str(fortran)], with_option, None, False)
    assert (tmp_path / 'plain' / 'final.bin').read_bytes() == (
        tmp_path / 'option' / 'final.bin'
    ).read_bytes()
    for out in ['z.out', 'xy.out', *extra]:
        assert (tmp_path / 'plain' / out).read_text() == (
            tmp_path / 'option' / out
        ).read_text()
