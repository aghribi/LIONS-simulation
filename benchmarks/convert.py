"""LIONS to OPAL-cycl and RF-Track conversions. One function per handoff.

Frames. LIONS: lab frame (x, y, z), z vertical, positions in m, momenta
p = beta gamma. OPAL-cycl: same axes, m and beta gamma
(`-trackOrbit.dat`), map Bz stored with the opposite sign. RF-Track: Y
vertical, positions in mm, momenta in MeV/c. The proper rotation (x, y,
z) -> (X, Y, Z) = (x, z, -y) keeps both frames right-handed; fields map
the same way.
"""

import numpy as np

C = 299792458.0


def write_opal_carbon(path, field):
    """OPAL CARBONCYCL map (mm, deg, kG) from a LIONS FieldMap.

    OPAL computes Bz = -(map value), so the map holds -Bz. It wants the
    azimuth from 0 without the 360 deg endpoint, and reads every value
    with fscanf("%16lE"): 16 characters at most.
    """
    ntheta = field.ntheta
    span_deg = np.degrees(field.dtheta) * (ntheta - 1)
    if not np.isclose(span_deg, 360.0):
        raise ValueError(
            f'OPAL CARBONCYCL needs a full 360 deg map, got {span_deg} deg'
        )
    bz_kg = -10.0 * field.bz[: field.nr, : ntheta - 1]
    with open(path, 'w') as f:
        f.write(
            f'{field.rmin * 1e3:.15g}\n{field.dr * 1e3:.15g}\n0.0\n'
            f'{np.degrees(field.dtheta):.15g}\n{ntheta - 1}\n{field.nr}\n'
        )
        for row in bz_kg:
            f.write(' '.join(f'{v:.9E}' for v in row) + '\n')


def read_opal_orbit(path):
    """Single-particle `-trackOrbit.dat`: positions (n, 3) in m and momenta
    (n, 3) in beta gamma.
    """
    d = np.loadtxt(path, comments='#', usecols=range(1, 7))
    return d[:, [0, 2, 4]], d[:, [1, 3, 5]]


def to_rftrack(pos, mom, mass_mev):
    """LIONS positions (m) and momenta (beta gamma), shape (n, 3), to
    RF-Track [X Px Y Py Z Pz].
    """
    pos_mm = np.asarray(pos) * 1e3
    p_mev = np.asarray(mom) * mass_mev
    return np.column_stack(
        [
            pos_mm[:, 0],
            p_mev[:, 0],
            pos_mm[:, 2],
            p_mev[:, 2],
            -pos_mm[:, 1],
            -p_mev[:, 1],
        ]
    )


def from_rftrack(ps, mass_mev):
    """RF-Track [X Px Y Py Z Pz] (mm, MeV/c) to LIONS positions (m) and
    momenta (beta gamma).
    """
    ps = np.atleast_2d(ps)
    pos = np.column_stack([ps[:, 0], -ps[:, 4], ps[:, 2]]) * 1e-3
    mom = np.column_stack([ps[:, 1], -ps[:, 5], ps[:, 3]]) / mass_mev
    return pos, mom


def to_rftrack_field(b):
    """LIONS field vector (Bx, By, Bz) to RF-Track (BX, BY, BZ)."""
    return b[0], b[2], -b[1]
