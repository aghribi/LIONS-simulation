"""Synthetic fields of cases B0 and B (DESIGN.md sections 5 and 12). No
RF-Track or OPAL dependency.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tests'))

import synthetic as syn  # noqa: E402 (after sys.path)

B0 = -1.0  # T, LIONS sign for protons
K, R0 = 0.5, 0.15  # 1/m, m: case B0 gradient
EPS = 0.3  # case B flutter
MAP = dict(rmin=0.05, rmax=0.6, nr=111)


def bz_b0(r, t):
    """Case B0 field Bz (T) at radius r (m) and angle t."""
    return B0 * (1 - K * (r - R0)) + 0 * t


def write_field(case, ntet):
    """Write the potential and bmap files of case B0 or B with ntet
    columns.
    """
    syn.write_pot_central('pot1.bin')
    syn.write_pot_outer('pot2.dat')
    if case == 'B0':
        syn.write_bmap(
            'bmap.dat',
            bz=bz_b0,
            dbzx=lambda r, t: -K * B0 * np.cos(t) + 0 * r,
            dbzy=lambda r, t: -K * B0 * np.sin(t) + 0 * r,
            ntet=ntet,
            **MAP,
        )
    else:
        syn.write_bmap(
            'bmap.dat',
            bz=lambda r, t: B0 * (1 + EPS * np.cos(4 * t)),
            dbzx=lambda r, t: 4 * EPS * B0 * np.sin(t) * np.sin(4 * t) / r,
            dbzy=lambda r, t: -4 * EPS * B0 * np.cos(t) * np.sin(4 * t) / r,
            ntet=ntet,
            **MAP,
        )


def b_b0(x, y, z):
    """Case B0 field in the LIONS frame (m, T), first order in z."""
    r = np.hypot(x, y)
    return z * -K * B0 * x / r, z * -K * B0 * y / r, bz_b0(r, 0.0)
