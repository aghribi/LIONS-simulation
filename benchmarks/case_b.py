"""Case B0 (smooth radial gradient, exact references) and case B (4-sector field on the bilinear map).

usage: uv run --group bench python benchmarks/case_b.py [outdir]
Writes case_b.parquet. OPAL runs in Docker; outdir must be visible to Docker.
"""
import contextlib
import importlib.metadata
import io
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'benchmarks'))

import lions  # noqa: E402
import lions.apertures  # noqa: E402
import lions.fields  # noqa: E402
import RF_Track as rft  # noqa: E402
import synthetic as syn  # noqa: E402
import orbits  # noqa: E402
from case_a import opal_version  # noqa: E402

B0 = -1.0                               # T, LIONS sign for protons
K, R0 = 0.5, 0.15                       # 1/m, m: case B0 gradient
EPS = 0.3                               # case B flutter
ENERGIES = [1.0, 2.0, 3.0, 4.0, 5.0]    # MeV
STEPS = [360, 720, 1440]
RFT_TOLERANCES = [1e-9, 1e-12]
NTHETA = {'B0': [361], 'B': [361, 721]}
MAP = dict(rmin=0.05, rmax=0.6, nr=111)
H = 1e-6                                # m, finite-difference step of the Jacobian
FFT_TURNS, FFT_OFFSET = 256, 1e-4       # LIONS second tune instrument, case B0


def bz_b0(r, t):
    return B0 * (1 - K * (r - R0)) + 0 * t


def write_field(case, ntet):
    syn.write_pot_central('pot1.bin')
    syn.write_pot_outer('pot2.dat')
    if case == 'B0':
        syn.write_bmap('bmap.dat', bz=bz_b0, dbzx=lambda r, t: -K * B0 * np.cos(t) + 0 * r,
                       dbzy=lambda r, t: -K * B0 * np.sin(t) + 0 * r, ntet=ntet, **MAP)
    else:
        syn.write_bmap('bmap.dat', bz=lambda r, t: B0 * (1 + EPS * np.cos(4 * t)),
                       dbzx=lambda r, t: 4 * EPS * B0 * np.sin(t) * np.sin(4 * t) / r,
                       dbzy=lambda r, t: -4 * EPS * B0 * np.cos(t) * np.sin(4 * t) / r, ntet=ntet, **MAP)


def b_b0(x, y, z):
    """Case B0 field in the LIONS frame (m, T), first order in z."""
    r = np.hypot(x, y)
    return z * -K * B0 * x / r, z * -K * B0 * y / r, bz_b0(r, 0.0)


def lions_state(w, n, nbrion=1):
    d = lions.State()
    d.run.swap_bytes = 'n'
    with contextlib.redirect_stdout(io.StringIO()):
        syn.write_deck('deck.dat', fhf=syn.revolution_frequency(w, B0), r0=syn.larmor_radius(w, B0),
                       w0=w, turns=1, npaspt=n, nbrion=nbrion)
        lions.load(d, 'deck.dat')
    return d


def check_parity(points=200, z_max=5e-3, seed=0):
    """RF-Track closed form against LIONS compute_bfield (DESIGN.md section 12). Raises on failure."""
    rng = np.random.default_rng(seed)
    d = lions_state(1.0, STEPS[0], nbrion=points)
    b, f = d.bunch, d.field
    r = rng.uniform(f.rmin + f.dr, f.rmax - 0.03, points)
    t = rng.uniform(0.0, 2 * np.pi, points)
    b.x[:], b.y[:], b.z[:] = r * np.cos(t), r * np.sin(t), rng.uniform(-z_max, z_max, points)
    lions.apertures.check_box(b, f, d.box)
    lions.fields.compute_bfield(b, f)
    bx, by, bz = b_b0(b.x, b.y, b.z)
    chord = K * abs(B0) * np.abs(b.z) * f.dtheta**2 / 8
    if not (np.all(np.abs(b.bz - bz) < 1e-13) and np.all(np.abs(b.bx - bx) <= chord + 1e-13)
            and np.all(np.abs(b.by - by) <= chord + 1e-13)):
        raise RuntimeError('case B0: RF-Track field differs from LIONS beyond the section 12 bounds')
    return {'parity_bz_max_abs': (np.max(np.abs(b.bz - bz)), 'T'),
            'parity_bxy_max_abs': (max(np.max(np.abs(b.bx - bx)), np.max(np.abs(b.by - by))), 'T')}


def exact_b0(p, mass_kg, charge):
    """Closed orbit radius and fractional linear tunes of case B0 (DESIGN.md section 12)."""
    rigidity = p * mass_kg * syn.C / charge
    a = 1 + K * R0
    r = (a - np.sqrt(a**2 - 4 * K * rigidity / abs(B0))) / (2 * K)
    n = K * r / (1 - K * (r - R0))
    fold = lambda q: np.arccos(np.cos(2 * np.pi * q)) / (2 * np.pi)   # range of the trace and the FFT
    return {'r_eo': (r, 'm'), 'qr': (fold(np.sqrt(1 - n)), '1'), 'qz': (fold(np.sqrt(n)), '1')}


def orbit_metrics(track, guess, p, rho, iters=8):
    """Closed orbit, then tunes and symplectic error from the one-turn Jacobian at it."""
    s, res = orbits.closed_orbit(track, guess, iters=iters)
    m = orbits.jacobian(track, s, H, H * p / rho)
    qr, qz = orbits.linear_tunes(m)
    eps_s, det_err = orbits.symplectic_error(m)
    return s, {'r_eo': (s[0], 'm'), 'pr_eo': (s[1], 'beta gamma'), 'residual': (res, 'm'),
               'qr': (qr, '1'), 'qz': (qz, '1'), 'symplectic_error': (eps_s, '1'),
               'det_error': (det_err, '1')}


def fft_tunes(track, s):
    """Tunes of the turn-by-turn r and z signals of two particles offset from the closed orbit s."""
    states = np.array([s + [FFT_OFFSET, 0, 0, 0], s + [0, 0, FFT_OFFSET, 0]])
    series = orbits.section_series(*track(states, 1.2 * FFT_TURNS))
    if min(len(x) for x in series) < FFT_TURNS:
        raise RuntimeError(f'fewer than {FFT_TURNS} section crossings')
    return {'qr_fft': (orbits.tune(series[0][:FFT_TURNS, 1]), '1'),
            'qz_fft': (orbits.tune(series[1][:FFT_TURNS, 3]), '1')}


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    versions = {'lions': importlib.metadata.version('lions'), 'opal': opal_version(),
                'rf_track': str(rft.RF_TRACK_VERSION), 'exact': None}
    commit = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True,
                            text=True).stdout.strip()
    records = []

    def add(case, ntet, code, integrator, tol, w, n, values):
        for metric, (value, unit) in values.items():
            records.append({'case': case, 'map_dtheta_deg': 360 / (ntet - 1), 'code': code,
                            'version': versions[code], 'integrator': integrator, 'tolerance': tol,
                            'energy_mev': w, 'steps_per_turn': n, 'metric': metric,
                            'value': float(value), 'unit': unit})

    for case, ntets in NTHETA.items():
        for ntet in ntets:
            work = outdir / f'{case}_{ntet}'
            work.mkdir(exist_ok=True)
            with contextlib.chdir(work):
                write_field(case, ntet)
                if case == 'B0':
                    add(case, ntet, 'exact', None, None, None, None, check_parity())
                for w in ENERGIES:
                    rho = syn.larmor_radius(w, B0)
                    d = lions_state(w, STEPS[0])
                    p, mass, charge = d.bunch.p[0], d.bunch.mass[0], d.bunch.charge[0]
                    f_rev = d.rf.frequency / d.rf.harmonic
                    guess = (rho, 0.0)
                    if case == 'B0':
                        exact = exact_b0(p, mass, charge)
                        add(case, ntet, 'exact', None, None, w, None, exact)
                        guess = (exact['r_eo'][0], 0.0)
                        for tol in RFT_TOLERANCES:
                            rt = orbits.RFTrackTracker(b_b0, d.field.rmax, mass, p, 1 / f_rev,
                                                       STEPS[-1], tol)
                            add(case, ntet, 'rf_track', 'rkf45', tol, w, None,
                                orbit_metrics(rt.track, guess, p, rho)[1])
                    for n in STEPS:
                        lt = orbits.LionsTracker('.', w, B0, n)
                        s, m = orbit_metrics(lt.track, guess, p, rho)
                        if case == 'B0' and n == STEPS[-1]:
                            m |= fft_tunes(lt.track, s)
                        add(case, ntet, 'lions', 'boris', None, w, n, m)
                        if ntet != ntets[0]:
                            continue
                        for integ in ('LF2', 'RK4'):
                            ot = orbits.OpalTracker(work / 'opal', d.field, mass, p, f_rev, n, integ)
                            add(case, ntet, 'opal', integ, None, w, n,
                                orbit_metrics(ot.track, (s[0], s[1]), p, rho, iters=4)[1])
                    print(f'{case} ntheta={ntet} {w} MeV done', file=sys.stderr)

    table = pa.Table.from_pylist(records).replace_schema_metadata({
        'case': f'B0: Bz = {B0} (1 - {K} (r - {R0})) T; B: Bz = {B0} (1 + {EPS} cos 4 theta) T; '
                'protons, theta = 0 section',
        'map': f"r = {MAP['rmin']} to {MAP['rmax']} m, {MAP['nr']} rows, 0 to 360 deg",
        'fd_step_m': str(H), 'fft_turns': str(FFT_TURNS), 'fft_offset_m': str(FFT_OFFSET),
        'git_commit': commit, 'machine': platform.platform(), 'python': platform.python_version(),
        **{f'version_{k}': v for k, v in versions.items() if v is not None}})
    pq.write_table(table, outdir / 'case_b.parquet')
    return table


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'benchmarks' / 'runs' / 'case_b'
    main(out.resolve())
