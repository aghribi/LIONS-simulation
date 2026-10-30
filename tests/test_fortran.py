"""Python port against the Fortran reference, single proton, same synthetic inputs.

Non-uniform case: 4-sector Bz with off-plane terms, vertical offset and angle,
accelerating dees through the outer-region potential. The central region is
not exercised (rmax1 = 0.1 mm).
"""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import synthetic as syn

REPO = Path(__file__).resolve().parent.parent
STDIN = 'n\ndeck.dat\n1000 2000 0 0\n0 0 0 0 0 0\n0\n'
W0 = 1.0       # MeV/u
B0 = -1.0      # T
V_DEE = 50e3   # V


@pytest.fixture(scope='session')
def fortran(tmp_path_factory):
    if shutil.which('gfortran') is None:
        pytest.skip('gfortran not found')
    build = tmp_path_factory.mktemp('fortran')
    exe = build / 'lions_f'
    # No FMA contraction: NumPy does not fuse, and fusing alone shifts the last bits.
    subprocess.run(['gfortran', '-O2', '-ffp-contract=off', '-J', str(build), '-o', str(exe),
                    str(REPO / 'lion_m.f90'), str(REPO / 'lions_2.0_ep.f90')],
                   check=True, capture_output=True)
    return exe


def write_case(tail):
    syn.write_bmap('bmap.dat',
                   bz=lambda r, t: B0 * (1 + 0.05 * np.cos(4 * t)),
                   d2bz=-2.0,
                   dbzx=lambda r, t: 0.3 * np.sin(4 * t),
                   dbzy=lambda r, t: 0.3 * np.cos(4 * t))
    syn.write_pot_outer('pot2.dat', pot=lambda x, y: 0.5 * (1 + np.tanh(y / 0.01)))
    syn.write_deck('deck.dat', fhf=syn.revolution_frequency(W0, B0),
                   r0=syn.larmor_radius(W0, B0), w0=W0, v_dee=V_DEE,
                   z0=0.002, pr0=10.0, pz0=5.0, tail=tail)


def test_python_matches_fortran(workdir, fortran):
    write_case(['no fortran.bin'])
    subprocess.run([str(fortran)], input=STDIN, text=True, check=True, capture_output=True)
    xf, pf = syn.read_fortran_bunch('fortran.bin')

    write_case(['no python.bin'])
    d = syn.quiet_run('deck.dat')
    xp, pp = np.fromfile('python.bin', '<f8').reshape(2, 1, 3)

    ga0 = syn.proton_gamma(W0)[0]
    assert abs((d.ga[0] - 1) / (ga0 - 1) - 1) > 1e-2, 'dee field not exercised'
    assert abs(d.z[0]) > 0 and abs(d.bx[0]) > 0, 'off-plane terms not exercised'
    np.testing.assert_allclose(xp, xf, rtol=0, atol=1e-12)
    np.testing.assert_allclose(pp, pf, rtol=1e-12, atol=0)
