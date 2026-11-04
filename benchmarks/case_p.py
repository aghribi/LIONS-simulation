"""Case P: synthetic machine for LIONS, F18 (docs/training/PHYSICAL.md) or NEON (docs/training/PHYSICAL_NEON.md). Synthetic: says nothing about a real machine.

usage: [CASE_P_MISSION=NEON] uv run --extra game --extra rl --group bench --group dev python benchmarks/case_p.py {design,dees,slip,offaxis,published} [--outroot DIR]
Energies are per nucleon unless named total; gains per turn are total, as T0's. CASE_P_MISSION is read at import, since the machine sets the module constants.
design: isochronous map from closed orbits, tunes (analytic field, Boris; LIONS Hermite at a few energies).
dees: gain per turn against T0 energy_gain. slip: RF phase slip from injection to extraction at the design setting.
offaxis: from the verdict start energy, centred and off the median plane. published: the design against CIME's published ranges (NEON).
"""

import argparse
import concurrent.futures
import contextlib
import io
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
sys.path.insert(0, str(ROOT / 'tests'))

import lions  # noqa: E402
import orbits  # noqa: E402
import synthetic as syn  # noqa: E402
from lions.game import machine as mc  # noqa: E402
from lions.game import physics  # noqa: E402
from lions.game.env import mission_by_id  # noqa: E402
from lions.state import MAP_NR, OUTER_NX, OUTER_NY  # noqa: E402
from provenance import provenance  # noqa: E402

STEPS = 360  # per turn, as in the timings of PHYSICAL.md
OPTIONS = {'interpolation': 'hermite', 'azimuth': 'atan2'}
GAP_DEG = 1.5  # tanh half width of a dee gap, seen from the centre: invented, as the whole central region
GAP_MIN = 5e-4  # m, smallest gap half width, resolved by the potential grid
WORKERS = 2
MISSION = os.environ.get('CASE_P_MISSION', 'F18')
REGISTRATION = (
    'docs/training/PHYSICAL_NEON.md'
    if MISSION == 'NEON'
    else 'docs/training/PHYSICAL.md'
)
M = mission_by_id(MISSION)
M.reset(0)
MACHINE, ION, HARMONIC = M.machine, M.ion, M.harmonic
A = ION.mass_number
W_EXT, DEE_KV = M.p['w_target'] / A, M.p['dee_kv']
F0 = mc.nominal(MACHINE, ION, W_EXT * A)[
    1
]  # Hz, revolution frequency at the design setting of T0
EPS = np.sqrt(
    2 * MACHINE.flutter
)  # F = <(B - <B>)^2> / <B>^2 = eps^2 / 2 for B = <B> (1 + eps cos N theta), lecture eq. for F
SPECIES = (
    syn.HYDROGEN if MISSION == 'F18' else A,
    A,
    ION.charge,
)  # atomic mass in the deck: T0's for an ion, A u
MASS_KG = syn.U * (SPECIES[0] - ION.charge * syn.ME_U)
CHARGE = ION.charge * syn.QE
BC = (
    2 * np.pi * F0 * MASS_KG / CHARGE
)  # T, field that turns a slow ion at F0, with LIONS constants


def gamma(w):
    return 1.0 + w * A * syn.MEV / (MASS_KG * syn.C**2)


def momentum(w):
    return np.sqrt(gamma(w) ** 2 - 1)


class Field:
    """B = -Bc gamma_iso(r) (1 + g(r)) (1 + eps cos N theta), LIONS sign. g: corrections at r_nodes, from closed orbits."""

    def __init__(self, r_nodes=None, g=None):
        self.n, self.r_nodes = MACHINE.sectors, r_nodes
        self.spline = (
            CubicSpline(r_nodes, g, bc_type='natural')
            if r_nodes is not None
            else None
        )

    def mean(self, r):
        """<B> and its first two radial derivatives."""
        a = 2 * np.pi * F0 / syn.C
        s = 1 - (a * r) ** 2
        gam = s**-0.5
        d1, d2 = (
            a * a * r * s**-1.5,
            a * a * s**-1.5 + 3 * a**4 * r * r * s**-2.5,
        )
        g, g1, g2 = 0.0, 0.0, 0.0
        if self.spline is not None:
            c = np.clip(r, self.r_nodes[0], self.r_nodes[-1])
            g, g1, g2 = self.spline(c), self.spline(c, 1), self.spline(c, 2)
        return (
            -BC * gam * (1 + g),
            -BC * (d1 * (1 + g) + gam * g1),
            -BC * (d2 * (1 + g) + 2 * d1 * g1 + gam * g2),
        )

    def nodes(self, r, t):
        """On the median plane: Bz, its Laplacian in (x, y), dBz/dx, dBz/dy."""
        b, b1, b2 = self.mean(r)
        c, s = np.cos(self.n * t), np.sin(self.n * t)
        bz, dr, dt = (
            b * (1 + EPS * c),
            b1 * (1 + EPS * c),
            -b * EPS * self.n * s,
        )
        lap = (b2 + b1 / r) * (1 + EPS * c) - b * EPS * self.n**2 * c / r**2
        return (
            bz,
            lap,
            np.cos(t) * dr - np.sin(t) / r * dt,
            np.sin(t) * dr + np.cos(t) / r * dt,
        )

    def __call__(self, x, y, z):
        r, t = np.hypot(x, y), np.arctan2(y, x)
        bz, lap, gx, gy = self.nodes(r, t)
        return (
            z * gx,
            z * gy,
            bz - z * z / 2 * lap,
        )  # lowest order off the plane, as LIONS


def closed_orbit(field, w, guess):
    """State, frequency [Hz], mean radius [m] of the closed orbit at kinetic energy w, Boris at the deck step."""
    p = momentum(w)
    track = orbits.BorisTracker(field, MASS_KG, CHARGE, p, F0, STEPS).track
    s, _ = orbits.closed_orbit(track, (guess, 0.0))
    t, pos, mom = track(s, 3)
    sec = orbits.section_series(t, pos, mom)[0]
    k = pos[:STEPS, 0]
    return (
        s,
        (len(sec) - 1) / (sec[-1, 0] - sec[0, 0]),
        float(np.hypot(k[:, 0], k[:, 1]).mean()),
        track,
        p,
    )


def isochronise(energies, iterations=8, tol=2e-6):
    field = Field()
    for it in range(iterations):
        rows = []
        for w in energies:
            p = momentum(w)
            guess = (
                rows[-1][0][0] * p / rows[-1][-1]
                if rows
                else p * MASS_KG * syn.C / (CHARGE * BC)
            )
            rows.append(closed_orbit(field, w, guess))
        f = np.array([r[1] for r in rows])
        rbar = np.array([r[2] for r in rows])
        worst = np.max(np.abs(f / F0 - 1))
        print(f'iteration {it}: max |f / F0 - 1| = {worst:.2e}', flush=True)
        if worst < tol:
            break
        # A local relative field change delta moves the orbit by -delta / gamma^2 (k = gamma^2 - 1): f changes by delta / gamma^2.
        delta = gamma(energies) ** 2 * (F0 / f - 1)
        old = field.spline(rbar) if field.spline is not None else 0.0
        field = Field(rbar, (1 + old) * (1 + delta) - 1)
    return field, rows


def rows_of(energies, field, rows):
    out = []
    for w, (s, f, rbar, track, p) in zip(energies, rows):
        m = orbits.jacobian(track, s, 1e-6, 1e-6 * p)
        tr_r, tr_z = np.trace(m[:2, :2]) / 2, np.trace(m[2:, 2:]) / 2
        k = gamma(w) ** 2 - 1
        out.append(
            {
                'code': 'boris_analytic',
                'w_mev': float(w),
                'r_co_m': float(s[0]),
                'rbar_m': rbar,
                'f_rel_error': float(f / F0 - 1),
                'half_trace_r': float(tr_r),
                'half_trace_z': float(tr_z),
                'qr_frac': float(
                    np.arccos(np.clip(tr_r, -1, 1)) / (2 * np.pi)
                ),
                'qz_frac': float(
                    np.arccos(np.clip(tr_z, -1, 1)) / (2 * np.pi)
                ),
                'qz_smooth': float(
                    np.sqrt(
                        physics.qz_squared(
                            k,
                            MACHINE.sectors,
                            MACHINE.flutter,
                            np.radians(MACHINE.spiral_deg),
                        )
                    )
                ),
                'qr_smooth': float(
                    np.sqrt(
                        physics.qr_squared(
                            k,
                            MACHINE.sectors,
                            MACHINE.flutter,
                            np.radians(MACHINE.spiral_deg),
                        )
                    )
                ),
            }
        )
    return out


def lions_check(job):
    """Closed orbit, frequency and tunes with LIONS Hermite on the written map, no E field."""
    workdir, w, guess = job
    track = orbits.LionsTracker(
        workdir, w, F0, STEPS, species=SPECIES, options=OPTIONS
    ).track
    s, _ = orbits.closed_orbit(track, (guess, 0.0))
    p = momentum(w)
    t, pos, mom = track(s, 3)
    sec = orbits.section_series(t, pos, mom)[0]
    m = orbits.jacobian(track, s, 1e-6, 1e-6 * p)
    tr_r, tr_z = np.trace(m[:2, :2]) / 2, np.trace(m[2:, 2:]) / 2
    return {
        'code': 'lions_hermite',
        'w_mev': float(w),
        'r_co_m': float(s[0]),
        'rbar_m': float('nan'),
        'f_rel_error': float(
            (len(sec) - 1) / (sec[-1, 0] - sec[0, 0]) / F0 - 1
        ),
        'half_trace_r': float(tr_r),
        'half_trace_z': float(tr_z),
        'qr_frac': float(np.arccos(np.clip(tr_r, -1, 1)) / (2 * np.pi)),
        'qz_frac': float(np.arccos(np.clip(tr_z, -1, 1)) / (2 * np.pi)),
        'qz_smooth': float('nan'),
        'qr_smooth': float('nan'),
    }


def dee_potential(x, y):
    """1 on the wedge around the x axis (y < 0 in the gap frame), 0 outside. The gap width follows the radius, so transit stays the same."""
    return 0.5 * (
        1 - np.tanh(y / np.maximum(x * np.tan(np.radians(GAP_DEG)), GAP_MIN))
    )


def write_machine(outdir, field):
    outdir.mkdir(parents=True, exist_ok=True)
    rmax = MACHINE.radius * 1.15
    with contextlib.chdir(outdir):
        nodes = lambda i: lambda r, t: field.nodes(r, t)[i]  # noqa: E731
        syn.write_bmap(
            'bmap.dat',
            nodes(0),
            nodes(1),
            nodes(2),
            nodes(3),
            rmin=0.005,
            rmax=rmax,
            nr=MAP_NR,
            ntet=361,
        )
        syn.write_pot_central('pot1.bin')
        width = 3 * rmax * np.tan(np.radians(GAP_DEG))
        syn.write_pot_outer(
            'pot2.dat',
            dee_potential,
            xmin=0.0,
            xmax=rmax,
            ymin=-width,
            ymax=width,
            nx=OUTER_NX,
            ny=OUTER_NY,
        )


def accelerate(
    workdir, w0, r0, pr0_mrad, phi0_deg, turns, one_turn=False, z0=0.0
):
    """Design setting, one ion from (r0, theta = 0): crossings of theta = 0 as (t [s], energy [MeV/u], r [m], z [m]), or the total gain [MeV] over one nominal turn."""
    with contextlib.chdir(workdir):
        syn.write_deck(
            'deck.dat',
            fhf=HARMONIC * F0,
            r0=r0,
            w0=w0,
            npaspt=STEPS,
            turns=turns,
            v_dee=DEE_KV * 1e3,
            teta1=MACHINE.dee_angle_deg / 2,
            harm=float(HARMONIC),
            z0=z0,
            pr0=pr0_mrad,
            phi0=phi0_deg,
            species=SPECIES,
            tail=orbits.deck_tail(OPTIONS),
        )
        rows = []

        def record(d):
            b = d.bunch
            rows.append(
                (
                    d.run.t,
                    b.x[0],
                    b.y[0],
                    b.z[0],
                    float(np.sqrt(b.px[0] ** 2 + b.py[0] ** 2 + b.pz[0] ** 2)),
                )
            )
            return bool(np.hypot(b.x[0], b.y[0]) < MACHINE.radius * 1.12)

        with contextlib.redirect_stdout(io.StringIO()):
            lions.run('deck.dat', on_step=record)
    a = np.array(rows)
    t, x, y, z, p = a.T
    if one_turn:
        return (
            (np.sqrt(1 + p[STEPS] ** 2) - np.sqrt(1 + p[0] ** 2))
            * MASS_KG
            * syn.C**2
            / syn.MEV
        )
    k = np.nonzero((y[:-1] < 0) & (y[1:] >= 0) & (x[:-1] > 0))[0]
    frac = -y[k] / (y[k + 1] - y[k])
    pc = p[k] + frac * (p[k + 1] - p[k])
    w = (np.sqrt(1 + pc * pc) - 1) * MASS_KG * syn.C**2 / syn.MEV / A
    return np.column_stack(
        [t[k] + frac * (t[k + 1] - t[k]), w, np.hypot(x[k], y[k]), z[k]]
    )


def gain_job(job):
    workdir, w, r0, phi = job
    return {
        'w_mev': float(w),
        'phi0_deg': float(phi),
        'gain_mev': float(
            accelerate(workdir, w, r0, 0.0, phi, 2, one_turn=True)
        ),
    }


def workdirs(outdir, n):
    dirs = []
    for i in range(n):
        d = outdir / f'work{i}'
        d.mkdir(exist_ok=True)
        for name in ('bmap.dat', 'pot1.bin', 'pot2.dat'):
            if not (d / name).exists():
                os.symlink(outdir / name, d / name)
        dirs.append(d)
    return dirs


def write(records, path, step, *inputs):
    meta = {
        'case': f'synthetic {MISSION}, {step}',
        'registration': REGISTRATION,
        'label': 'SYNTHETIC',
        **provenance(*inputs),
        'machine': platform.platform(),
        'python': platform.python_version(),
        'flutter_formula': 'F = eps^2 / 2',
        'eps': str(EPS),
        'gap_half_width_deg': str(GAP_DEG),
        'steps_per_turn': str(STEPS),
        'interpolation': 'hermite',
        'units': 'w_mev MeV, r m, gain MeV, phase deg, f_rel_error 1',
    }
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta), path
    )


def load_field(outdir):
    t = pq.read_table(outdir / 'case_p_field.parquet').to_pylist()
    return Field(
        np.array([r['r_m'] for r in t]), np.array([r['g'] for r in t])
    )


def design(outdir):
    energies = np.geomspace(MACHINE.w_injection, W_EXT * 1.03, 30)
    field, rows = isochronise(energies)
    outdir.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist(
            [
                {'r_m': float(a), 'g': float(b)}
                for a, b in zip(field.r_nodes, field.spline(field.r_nodes))
            ]
        ),
        outdir / 'case_p_field.parquet',
    )
    write_machine(outdir, field)
    records = rows_of(energies, field, rows)
    picks = energies[[0, 8, 15, 22, 29]]
    guess = {float(w): r[0][0] for w, r in zip(energies, rows)}
    dirs = workdirs(outdir, WORKERS)
    jobs = [
        (dirs[i % WORKERS], float(w), guess[float(w)])
        for i, w in enumerate(picks)
    ]
    with concurrent.futures.ThreadPoolExecutor(1) as pool:
        records += list(pool.map(lions_check, jobs))
    write(
        records,
        outdir / 'case_p_orbits.parquet',
        'closed orbits and tunes',
        outdir / 'bmap.dat',
    )
    for r in records:
        print(
            {
                k: (round(v, 7) if isinstance(v, float) else v)
                for k, v in r.items()
            },
            flush=True,
        )


def dees(outdir):
    t = {
        r['w_mev']: r['r_co_m']
        for r in pq.read_table(outdir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    }
    keys = sorted(t)
    energies = [keys[i] for i in (0, 5, 10, 15, 20, 25, 28)]
    dirs = workdirs(outdir, WORKERS)
    jobs = [
        (dirs[0], w, t[w], float(phi))
        for w in energies
        for phi in range(0, 360, 60)
    ]
    per = [[] for _ in dirs]
    for i, j in enumerate(jobs):
        per[i % WORKERS].append((dirs[i % WORKERS], *j[1:]))
    with concurrent.futures.ProcessPoolExecutor(WORKERS) as pool:
        rows = [r for part in pool.map(serial_gains, per) for r in part]
    expected = float(mc.energy_gain(MACHINE, ION, DEE_KV, HARMONIC, 0.0))
    out = []
    for w in energies:
        sel = sorted(
            (r for r in rows if r['w_mev'] == w), key=lambda r: r['phi0_deg']
        )
        phi = np.radians([r['phi0_deg'] for r in sel])
        g = np.array([r['gain_mev'] for r in sel])
        c = np.linalg.lstsq(
            np.column_stack([np.cos(phi), np.sin(phi), np.ones_like(phi)]),
            g,
            rcond=None,
        )[0]
        amp = float(np.hypot(c[0], c[1]))
        crest = float(np.degrees(np.arctan2(c[1], c[0])))
        res = float(
            np.max(
                np.abs(
                    np.column_stack(
                        [np.cos(phi), np.sin(phi), np.ones_like(phi)]
                    )
                    @ c
                    - g
                )
            )
        )
        out.append(
            {
                'w_mev': w,
                'gain_amp_mev': amp,
                't0_gain_mev': expected,
                'ratio': amp / expected,
                'crest_phi0_deg': crest,
                'offset_mev': float(c[2]),
                'fit_residual_mev': res,
            }
        )
        print(out[-1], flush=True)
    write(
        rows,
        outdir / 'case_p_gain_raw.parquet',
        'dees, raw',
        outdir / 'bmap.dat',
        outdir / 'pot2.dat',
    )
    write(
        out,
        outdir / 'case_p_gain.parquet',
        'dees',
        outdir / 'bmap.dat',
        outdir / 'pot2.dat',
    )


def serial_gains(jobs):
    return [gain_job(j) for j in jobs]


def slip(
    outdir, w0=None, r0=None, pr0=0.0, phi0=None, tag='injection', z0=0.0
):
    rec = [
        r
        for r in pq.read_table(outdir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    ]
    w0 = rec[0]['w_mev'] if w0 is None else w0
    r0 = (
        next(r['r_co_m'] for r in rec if abs(r['w_mev'] - w0) < 1e-9)
        if r0 is None
        else r0
    )
    gains = pq.read_table(outdir / 'case_p_gain.parquet').to_pylist()
    phi0 = (
        float(np.median([g['crest_phi0_deg'] for g in gains]))
        if phi0 is None
        else phi0
    )
    turns = int(
        2
        * W_EXT
        * A
        / float(mc.energy_gain(MACHINE, ION, DEE_KV, HARMONIC, 0.0))
    )
    s = accelerate(workdirs(outdir, 1)[0], w0, r0, pr0, phi0, turns, z0=z0)
    phase = 360.0 * HARMONIC * (F0 * s[:, 0] - np.arange(len(s)))
    phase -= phase[0]
    out = [
        {
            'tag': tag,
            'turn': i,
            't_s': float(a[0]),
            'w_mev': float(a[1]),
            'r_m': float(a[2]),
            'z_m': float(a[3]),
            'phase_deg': float(b),
        }
        for i, (a, b) in enumerate(zip(s, phase))
    ]
    write(
        out,
        outdir / f'case_p_slip_{tag}.parquet',
        'phase slip',
        outdir / 'bmap.dat',
        outdir / 'pot2.dat',
    )
    end = np.nonzero(s[:, 1] >= W_EXT)[0]
    k = end[0] if len(end) else len(s) - 1
    print(
        f'{tag}: {len(s)} crossings, reached r = {s[k, 2]:.4f} m at turn {k}, w = {s[k, 1]:.3f} MeV; max |slip| to there {np.max(np.abs(phase[: k + 1])):.2f} deg, end slip {phase[k]:.2f} deg',
        flush=True,
    )
    print('w', np.round(s[: k + 1 : max(k // 12, 1), 1], 3).tolist())
    print('phase', np.round(phase[: k + 1 : max(k // 12, 1)], 1).tolist())
    print(
        f'{tag}: max |z| at the section {1e3 * np.max(np.abs(s[: k + 1, 3])):.3f} mm, final energy {s[k, 1]:.3f} MeV'
    )
    gain = np.diff(s[: k + 1, 1]) * A
    print(
        'gain per turn: median',
        np.median(gain),
        'max',
        gain.max(),
        'T0 peak',
        mc.energy_gain(MACHINE, ION, DEE_KV, HARMONIC, 0.0),
    )


def start_energy(gains, rec):
    """The closed-orbit energy of the lowest gain row from which every gain ratio is within 5%."""
    keep = [g['w_mev'] for g in sorted(gains, key=lambda g: g['w_mev'])]
    for g in sorted(gains, key=lambda g: -g['w_mev']):
        if abs(g['ratio'] - 1) >= 0.05:
            keep = [w for w in keep if w > g['w_mev']]
            break
    return min(rec, key=lambda w: abs(w - keep[0]))


def published(outdir):
    """Registered checks against CIME's published ranges, PHYSICAL_NEON.md. Qr close to 1 is read from its fractional part."""
    rows = [
        r
        for r in pq.read_table(outdir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    ]
    inner = [r for r in rows if 0.2 <= r['rbar_m'] <= 1.4]
    qr, qz = [r['qr_frac'] for r in inner], [r['qz_frac'] for r in inner]
    mean = -load_field(outdir).mean(MACHINE.radius)[0]
    print(
        f'{len(inner)} orbits at 0.2 to 1.4 m. |Qr - 1|: max {max(qr):.4f}, limit 0.04. Qz: {min(qz):.3f} to {max(qz):.3f}, limits 0.2 and 0.45'
    )
    print(
        f'<B> at r_ext {mean:.4f} T, range 0.75 to 1.56. B rho {mean * MACHINE.radius:.4f} T m, limit 2.344. f_rf {HARMONIC * F0 / 1e6:.4f} MHz, range 9.6 to 14.5'
    )
    print(
        f'peak B {mean * (1 + EPS):.3f} T, published maximum on the hills 2.1 T'
    )


def offaxis(outdir):
    """Gate of PHYSICAL.md: from the start energy at the crest, z0 = 0, T0's z_injection and twice it (T0 spreads z0 by 30% rms)."""
    rec = {
        r['w_mev']: r['r_co_m']
        for r in pq.read_table(outdir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    }
    gains = pq.read_table(outdir / 'case_p_gain.parquet').to_pylist()
    w0 = start_energy(gains, rec)
    phi0 = min(gains, key=lambda g: abs(g['w_mev'] - w0))['crest_phi0_deg']
    for factor in (0.0, 1.0, 2.0):
        slip(
            outdir,
            w0,
            rec[w0],
            phi0=phi0,
            tag=f'offaxis_{factor:g}',
            z0=factor * MACHINE.z_injection * 1e-3,
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'step', choices=('design', 'dees', 'slip', 'offaxis', 'published')
    )
    parser.add_argument(
        '--outroot',
        type=Path,
        default=ROOT
        / 'benchmarks'
        / 'runs'
        / ('case_p' if MISSION == 'F18' else f'case_p_{MISSION}'),
    )
    args = parser.parse_args()
    outdir = args.outroot.resolve()
    {
        'design': design,
        'dees': dees,
        'slip': slip,
        'offaxis': offaxis,
        'published': published,
    }[args.step](outdir)
