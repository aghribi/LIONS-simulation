"""Synthetic LIONS inputs: field maps, potentials and decks written on the
fly.
"""

import contextlib
import io
import struct

import numpy as np

import lions

# Constants as used by lions.init
U = 1.6605402e-27
QE = 1.60217733e-19
C = 299792458.0
MEV = 1.60217733e-13
ME_U = 0.0005486

HYDROGEN = 1.00782503  # atomic mass (u)


def proton_gamma(w_mev):
    m0 = U * (HYDROGEN - ME_U)
    return 1.0 + w_mev * MEV / (m0 * C**2), m0


def larmor_radius(w_mev, b):
    ga, m0 = proton_gamma(w_mev)
    return m0 * C * np.sqrt(ga**2 - 1.0) / (QE * abs(b))


def revolution_frequency(w_mev, b):
    ga, m0 = proton_gamma(w_mev)
    return QE * abs(b) / (2 * np.pi * ga * m0)


def write_bmap(
    path,
    bz,
    d2bz=0.0,
    dbzx=0.0,
    dbzy=0.0,
    rmin=0.05,
    rmax=0.5,
    nr=46,
    ntet=361,
    theta_max=2 * np.pi,
):
    """Field map over 0..theta_max (360 deg, or 45 deg with ntet = 91 for a
    1/8 map).

    Each field is a constant or a function f(r, theta).
    """
    r = np.linspace(rmin, rmax, nr)[:, None]
    teta = np.linspace(0.0, theta_max, ntet)[None, :]
    with open(path, 'w') as f:
        f.write('synthetic\n')
        f.write(f'{rmin} {rmax} {nr} {ntet}\n')
        for field in (bz, d2bz, dbzx, dbzy):
            values = np.broadcast_to(
                field(r, teta) if callable(field) else field, (nr, ntet)
            )
            for row in values:
                f.write(' '.join(f'{v:.15e}' for v in row) + '\n')


def fortran_record(payload):
    """Fortran unformatted sequential record, little-endian 4-byte markers."""
    return (
        struct.pack('<i', len(payload))
        + payload
        + struct.pack('<i', len(payload))
    )


def write_pot_central(
    path, pot=None, d2pot=None, half_width_mm=10.0, nx=5, ny=5
):
    """Central-region potential pot(x, y) and d2pot(x, y), x and y in m.
    Zero by default.
    """
    h = half_width_mm
    x = np.linspace(-h, h, nx)[:, None] / 1000.0
    y = np.linspace(-h, h, ny)[None, :] / 1000.0
    with open(path, 'wb') as f:
        f.write(fortran_record(struct.pack('<4d2i', -h, h, -h, h, nx, ny)))
        for field in (pot, d2pot):
            values = np.broadcast_to(
                field(x, y) if field is not None else 0.0, (nx, ny)
            )
            for row in values:
                f.write(
                    fortran_record(np.ascontiguousarray(row, '<f8').tobytes())
                )


def write_bunch(python_path, fortran_path, pos, mom):
    """Same bunch in the Python (raw doubles) and Fortran (records)
    formats. pos, mom: (n, 3).
    """
    pos = np.ascontiguousarray(pos, '<f8')
    mom = np.ascontiguousarray(mom, '<f8')
    with open(python_path, 'wb') as f:
        f.write(pos.tobytes() + mom.tobytes())
    with open(fortran_path, 'wb') as f:
        f.write(fortran_record(pos.tobytes()) + fortran_record(mom.tobytes()))


def write_markers(path, traps):
    """Marker file with trap circles (x, y, r in mm)."""
    with open(path, 'w') as f:
        f.write('0 0 0 0\n0 0\n0 0\n')
        f.write('cercle\n')
        f.write(f'{len(traps)}\n')
        for x, y, r in traps:
            f.write(f'{x} {y} {r} 0 360 360 -100\n')
        f.write('fintyp\n')


def write_pot_outer(
    path, pot=None, xmin=-1.0, xmax=1.0, ymin=-1.0, ymax=1.0, nx=41, ny=201
):
    """Outer-region potential pot(x, y) in the rotated gap frame (m), and
    d2pot = 0.
    """
    x = np.linspace(xmin, xmax, nx)[:, None]
    y = np.linspace(ymin, ymax, ny)[None, :]
    values = np.broadcast_to(pot(x, y) if pot is not None else 0.0, (nx, ny))
    with open(path, 'w') as f:
        f.write(f'{xmin} {xmax} {ymin} {ymax} {nx} {ny}\n')
        for row in values:
            f.write(' '.join(f'{v:.15e}' for v in row) + '\n')
        for _ in range(nx):
            f.write(' '.join(['0.0'] * ny) + '\n')


def write_deck(
    path,
    *,
    fhf,
    r0,
    w0,
    npaspt=360,
    turns=10,
    v_dee=0.0,
    teta1=20.0,
    harm=1.0,
    z0=0.0,
    pr0=0.0,
    pz0=0.0,
    phi0=0.0,
    rmax1=0.0001,
    nbrion=1,
    extraction=('1 0 0 0',),
    correlations=('0',),
    tail=(),
    species=(HYDROGEN, 1, 1),
):
    """One-family deck, protons by default. Outputs z.out, xy.out.

    species: (atomic mass in u, mass number, charge state).

    extraction, correlations: raw deck lines from 'rex tet1ex tet2ex
    eex' and from 'iflcor'.
    """
    lines = [
        "'synthetic'",
        f'{npaspt}',
        f'0.25 {turns}',
        f'{fhf:.17g} {-v_dee:.17g} {teta1} {harm} 0.0',
        '0 0',
        f'{v_dee:.17g} {v_dee:.17g}',
        '0',
        '0',
        f'-0.1 0.1 {rmax1}',
        *extraction,
        '1',
        f'{nbrion} {species[0]} {species[1]} {species[2]}',
        f'{r0:.17g} 0 {z0:.17g} {pr0:.17g} {pz0:.17g} {w0:.17g} {phi0:.17g}',
        *correlations,
        f'{turns * npaspt} {npaspt}',
        'z.out xy.out',
        "'bmap.dat' 'pot1.bin' 'pot2.dat'",
        *tail,
    ]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def quiet_run(deck, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return lions.run(deck, **kwargs)


def numeric_rows(path):
    """Rows of a text output that hold only numbers (headers and labels
    dropped).
    """
    rows = []
    for line in open(path):
        try:
            rows.append([float(t.replace('D', 'E')) for t in line.split()])
        except ValueError:
            continue
    return [r for r in rows if r]


def read_fortran_bunch(path, n=1):
    """Bunch file written by the Fortran (two unformatted records): x,y,z
    then px,py,pz.
    """
    raw = open(path, 'rb').read()
    out, pos = [], 0
    for _ in range(2):
        size = struct.unpack_from('<i', raw, pos)[0]
        out.append(np.frombuffer(raw, '<f8', size // 8, pos + 4).reshape(n, 3))
        pos += size + 8
    return out[0], out[1]
