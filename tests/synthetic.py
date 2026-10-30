"""Synthetic LIONS inputs: field maps, potentials and decks written on the fly."""
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


def write_bmap(path, bz, d2bz=0.0, dbzx=0.0, dbzy=0.0,
               rmin=0.05, rmax=0.5, nr=46, ntet=361):
    """Full 360 deg map. Each field is a constant or a function f(r, teta)."""
    r = np.linspace(rmin, rmax, nr)[:, None]
    teta = np.linspace(0.0, 2 * np.pi, ntet)[None, :]
    with open(path, 'w') as f:
        f.write('synthetic\n')
        f.write(f'{rmin} {rmax} {nr} {ntet}\n')
        for field in (bz, d2bz, dbzx, dbzy):
            values = np.broadcast_to(field(r, teta) if callable(field) else field, (nr, ntet))
            for row in values:
                f.write(' '.join(f'{v:.15e}' for v in row) + '\n')


def write_pot_central(path, nx=5, ny=5):
    """Zero central potential, Fortran unformatted sequential, native byte order."""
    def record(payload):
        return struct.pack('<i', len(payload)) + payload + struct.pack('<i', len(payload))

    with open(path, 'wb') as f:
        f.write(record(struct.pack('<4d2i', -10.0, 10.0, -10.0, 10.0, nx, ny)))
        row = np.zeros(ny).astype('<f8').tobytes()
        for _ in range(2 * nx):
            f.write(record(row))


def write_pot_outer(path, pot=None, xmin=-1.0, xmax=1.0, ymin=-1.0, ymax=1.0, nx=41, ny=201):
    """Outer-region potential pot(x, y) in the rotated gap frame (m), and d2pot = 0."""
    x = np.linspace(xmin, xmax, nx)[:, None]
    y = np.linspace(ymin, ymax, ny)[None, :]
    values = np.broadcast_to(pot(x, y) if pot is not None else 0.0, (nx, ny))
    with open(path, 'w') as f:
        f.write(f'{xmin} {xmax} {ymin} {ymax} {nx} {ny}\n')
        for row in values:
            f.write(' '.join(f'{v:.15e}' for v in row) + '\n')
        for _ in range(nx):
            f.write(' '.join(['0.0'] * ny) + '\n')


def write_deck(path, *, fhf, r0, w0, npaspt=360, turns=10, v_dee=0.0, teta1=20.0,
               harm=1.0, z0=0.0, pr0=0.0, pz0=0.0, phi0=0.0, tail=()):
    """Single-proton deck, no extraction, no correlations. Outputs z.out, xy.out."""
    lines = [
        "'synthetic'", f'{npaspt}', f'0.25 {turns}',
        f'{fhf:.17g} {-v_dee:.17g} {teta1} {harm} 0.0', '0 0', f'{v_dee:.17g} {v_dee:.17g}',
        '0', '0', '-0.1 0.1 0.0001', '1 0 0 0', '1',
        f'1 {HYDROGEN} 1 1', f'{r0:.17g} 0 {z0:.17g} {pr0:.17g} {pz0:.17g} {w0:.17g} {phi0:.17g}',
        '0', f'{turns * npaspt} {npaspt}', 'z.out xy.out',
        "'bmap.dat' 'pot1.bin' 'pot2.dat'",
        *tail,
    ]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def quiet_run(deck, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return lions.run(deck, **kwargs)


def read_fortran_bunch(path, n=1):
    """Bunch file written by the Fortran (two unformatted records): x,y,z then px,py,pz."""
    raw = open(path, 'rb').read()
    out, pos = [], 0
    for _ in range(2):
        size = struct.unpack_from('<i', raw, pos)[0]
        out.append(np.frombuffer(raw, '<f8', size // 8, pos + 4).reshape(n, 3))
        pos += size + 8
    return out[0], out[1]
