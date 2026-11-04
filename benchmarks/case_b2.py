"""Case B discriminating runs (DESIGN.md section 14): B-smooth (RF-Track,
no map) and B-gradient.

usage: uv run --group bench python benchmarks/case_b2.py [outdir] Writes
case_b2.parquet. OPAL runs in Docker; outdir must be visible to Docker.
"""

import contextlib
import importlib.metadata
import platform
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
import RF_Track as rft  # noqa: N813 (external name)
import synthetic as syn  # noqa: E402
from provenance import provenance  # noqa: E402
import orbits  # noqa: E402
from case_a import opal_version  # noqa: E402
from case_b import ENERGIES, STEPS, lions_state, exact_b0, orbit_metrics  # noqa: E402
from fields_b import B0, R0, EPS, MAP  # noqa: E402

K_GRADIENT = 0.5  # 1/m
TOL = 1e-12
NTHETA = 361
RUNS = {
    'B-smooth': dict(k=0.0, energies=[1.0], lions_opal=False),
    'B-gradient': dict(k=K_GRADIENT, energies=ENERGIES, lions_opal=True),
}


def fields(k):
    """Bz, dBz/dx, dBz/dy as functions of (r, theta) for B0 g(r) f(theta)
    (DESIGN.md section 14).
    """

    def g(r):
        return 1 - k * (r - R0)

    def f(t):
        return 1 + EPS * np.cos(4 * t)

    return (
        lambda r, t: B0 * g(r) * f(t),
        lambda r, t: (
            B0
            * (
                -k * f(t) * np.cos(t)
                + 4 * EPS * g(r) * np.sin(4 * t) * np.sin(t) / r
            )
        ),
        lambda r, t: (
            B0
            * (
                -k * f(t) * np.sin(t)
                - 4 * EPS * g(r) * np.sin(4 * t) * np.cos(t) / r
            )
        ),
    )


def closed_form(k):
    """Field b(x, y, z) in the LIONS frame (m, T), first order in z."""
    bz, dx, dy = fields(k)

    def b(x, y, z):
        r, t = np.hypot(x, y), np.arctan2(y, x)
        return z * dx(r, t), z * dy(r, t), bz(r, t)

    return b


def check_parity(b, k, points=200, z_max=5e-3, seed=0):
    """Closed form against LIONS compute_bfield at the map nodes in theta.
    Raises on failure.

    Bz is exact there. Bx, By carry the linear interpolation in r of
    g(r) / r, at most
    |z| dr^2 / 8 * 8 eps |B0| (1 + k r0) / r_i^3, r_i the inner node of the cell.
    """
    rng = np.random.default_rng(seed)
    d = lions_state(1.0, STEPS[0], nbrion=points)
    bunch, f = d.bunch, d.field
    r = rng.uniform(f.rmin + f.dr, f.rmax - 0.03, points)
    t = rng.integers(0, f.ntheta - 1, points) * f.dtheta
    bunch.x[:], bunch.y[:], bunch.z[:] = (
        r * np.cos(t),
        r * np.sin(t),
        rng.uniform(-z_max, z_max, points),
    )
    lions.apertures.check_box(bunch, f, d.box)
    lions.fields.compute_bfield(bunch, f)
    bx, by, bz = b(bunch.x, bunch.y, bunch.z)
    r_i = f.rmin + np.floor((r - f.rmin) / f.dr) * f.dr
    bound = (
        np.abs(bunch.z) * f.dr**2 * EPS * abs(B0) * (1 + k * R0) / r_i**3
        + 1e-13
    )
    dz, dxy = (
        np.abs(bz - bunch.bz),
        np.maximum(np.abs(bx - bunch.bx), np.abs(by - bunch.by)),
    )
    if np.max(dz) > 1e-13 or np.any(dxy > bound):
        raise RuntimeError(
            f'closed form differs from LIONS at the map nodes: Bz {np.max(dz)} T, '
            f'Bx, By {np.max(dxy / bound)} of the bound'
        )
    return {
        'parity_bz_max_abs': (np.max(dz), 'T'),
        'parity_bxy_max_abs': (np.max(dxy), 'T'),
        'parity_bxy_max_over_bound': (np.max(dxy / bound), '1'),
    }


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    versions = {
        'lions': importlib.metadata.version('lions'),
        'opal': opal_version(),
        'rf_track': str(rft.RF_TRACK_VERSION),
        'exact': None,
    }
    records = []

    def add(case, code, integrator, w, n, values):
        for metric, (value, unit) in values.items():
            records.append(
                {
                    'case': case,
                    'code': code,
                    'version': versions[code],
                    'integrator': integrator,
                    'energy_mev': w,
                    'steps_per_turn': n,
                    'metric': metric,
                    'value': float(value),
                    'unit': unit,
                }
            )

    for case, run in RUNS.items():
        work = outdir / case
        work.mkdir(exist_ok=True)
        bz, dx, dy = fields(run['k'])
        b = closed_form(run['k'])
        with contextlib.chdir(work):
            syn.write_pot_central('pot1.bin')
            syn.write_pot_outer('pot2.dat')
            syn.write_bmap(
                'bmap.dat', bz=bz, dbzx=dx, dbzy=dy, ntet=NTHETA, **MAP
            )
            add(case, 'exact', None, None, None, check_parity(b, run['k']))
            for w in run['energies']:
                rho = syn.larmor_radius(w, B0)
                d = lions_state(w, STEPS[0])
                p, mass, charge = (
                    d.bunch.p[0],
                    d.bunch.mass[0],
                    d.bunch.charge[0],
                )
                f_rev = d.rf.frequency / d.rf.harmonic
                guess = (
                    exact_b0(p, mass, charge)['r_eo'][0] if run['k'] else rho,
                    0.0,
                )
                for n in STEPS:
                    rt = orbits.RFTrackTracker(
                        b, d.field.rmax, mass, p, 1 / f_rev, n, TOL
                    )
                    s, m = orbit_metrics(rt.track, guess, p, rho)
                    add(case, 'rf_track', 'rkf45', w, n, m)
                    if not run['lions_opal']:
                        continue
                    lt = orbits.LionsTracker(
                        '.', w, syn.revolution_frequency(w, B0), n
                    )
                    s, m = orbit_metrics(lt.track, guess, p, rho)
                    add(case, 'lions', 'boris', w, n, m)
                    for integ in ('LF2', 'RK4'):
                        ot = orbits.OpalTracker(
                            work / 'opal', d.field, mass, p, f_rev, n, integ
                        )
                        add(
                            case,
                            'opal',
                            integ,
                            w,
                            n,
                            orbit_metrics(
                                ot.track, (s[0], s[1]), p, rho, iters=4
                            )[1],
                        )
                print(f'{case} {w} MeV done', file=sys.stderr)

    table = pa.Table.from_pylist(records).replace_schema_metadata(
        {
            'case': f'Bz = {B0} (1 - k (r - {R0})) (1 + {EPS} cos 4 theta) T, protons, theta = 0 section; '
            f'B-smooth k = 0 (RF-Track, closed form), B-gradient k = {K_GRADIENT} /m',
            'map': f"r = {MAP['rmin']} to {MAP['rmax']} m, {MAP['nr']} rows, {NTHETA} columns over 360 deg",
            'rf_track_tolerance': str(TOL),
            **provenance(),
            'machine': platform.platform(),
            'python': platform.python_version(),
            **{
                f'version_{k}': v for k, v in versions.items() if v is not None
            },
        }
    )
    pq.write_table(table, outdir / 'case_b2.parquet')
    return table


if __name__ == '__main__':
    out = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else ROOT / 'benchmarks' / 'runs' / 'case_b2'
    )
    main(out.resolve())
