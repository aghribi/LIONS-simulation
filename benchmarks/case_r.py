"""Case R: reduced model against LIONS and T0 on a mission's tunes,
docs/training/PHYSICAL.md (F18) and PHYSICAL_NEON.md (NEON). Synthetic
map: says nothing about a real machine.

usage: [CASE_P_MISSION=NEON] uv run --extra game --extra rl --group
bench --group dev python benchmarks/case_r.py {tables,check,compare}
[--mapdir DIR] [--verdict DIR] [--outroot DIR] [--workers N] tables:
LIONS Hermite closed orbits at the energies of case P, no RF, into the
scenario table. check: the rules of the model, each against its
reference. compare: reduced model and T0 against LIONS on the tunes of
the verdict, with pass times. Energies in the table are total, MeV.
"""

import argparse
import concurrent.futures
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

import case_p  # noqa: E402
import case_v  # noqa: E402
import synthetic as syn  # noqa: E402
from lions.game import machine as mc  # noqa: E402
from lions.game import reduced  # noqa: E402
from provenance import provenance  # noqa: E402

TABLE = (
    ROOT
    / 'src'
    / 'lions'
    / 'game'
    / 'scenarios'
    / f'{case_p.MISSION.lower()}_synthetic_orbits.parquet'
)
R_MIN = 0.005  # m, end of the map grid, as the verdict
NEON = case_p.MISSION == 'NEON'
THRESHOLD = {
    'transmission': 1 / case_v.N,
    'energy': 0.003 if NEON else 0.005,
    'spread': 0.001 if NEON else 0.005,
    'phase_deg': 6.91 if NEON else 5.06,
    'pass_s': 0.15,
}
START_MEV, EXTRACTION, PHASE_RMS_DEG = (
    (1.0, 'energy', 10.0) if NEON else (None, 'mean', None)
)


def orbit_job(job):
    """Return the LIONS closed-orbit check of (workdir, w, guess)."""
    workdir, w, guess = job
    return case_p.lions_check((workdir, w, guess))


def tables(mapdir, workers):
    """Write the scenario table from LIONS closed orbits, no RF."""
    mp = case_v.Map(mapdir, START_MEV, EXTRACTION)
    rec = [
        r
        for r in pq.read_table(mapdir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    ]
    dirs = [mp.workdir(mapdir / f'table_{i}') for i in range(len(rec))]
    with concurrent.futures.ProcessPoolExecutor(workers) as pool:
        rows = list(
            pool.map(
                orbit_job,
                [(d, r['w_mev'], r['r_co_m']) for d, r in zip(dirs, rec)],
            )
        )
    gains = [
        g
        for g in pq.read_table(mapdir / 'case_p_gain.parquet').to_pylist()
        if g['w_mev'] >= mp.w0
    ]
    lnw, per_kv = (
        np.log([g['w_mev'] for g in gains]),
        np.array([g['gain_amp_mev'] for g in gains]) / case_p.DEE_KV,
    )
    out = []
    for r, b in zip(rows, rec):
        w = r['w_mev']
        g = 1 + w * case_p.A / case_v.REST_MEV
        p = np.sqrt(g * g - 1)
        f = case_p.F0 * (1 + r['f_rel_error'])
        out.append(
            {
                'w_mev': w * case_p.A,
                'p': float(p),
                'path_m': float(p / g * syn.C / f),
                'f_hz': float(f),
                'r_section_m': r['r_co_m'],
                'qz': r['qz_frac'],
                'qr': 1 + r['qr_frac'],
                'gain_mev_per_kv': float(np.interp(np.log(w), lnw, per_kv)),
                'gain_measured': bool(w >= mp.w0),
            }
        )
        print(
            {
                k: out[-1][k]
                for k in ('w_mev', 'f_hz', 'r_section_m', 'qz', 'qr')
            },
            'boris r',
            b['r_co_m'],
            flush=True,
        )
    meta = {
        'case': (
            f'synthetic {case_p.MISSION}, '
            'closed orbits of LIONS Hermite, no RF'
        ),
        'label': 'SYNTHETIC',
        'registration': case_p.REGISTRATION,
        **provenance(*mp.inputs),
        'machine': platform.platform(),
        'units': (
            'w_mev MeV, p beta gamma, path_m m, f_hz Hz, '
            'r_section_m m at theta = 0, '
            'gain_mev_per_kv MeV/kV peak per turn; '
            'below w_start the gain is held, not measured'
        ),
        'rest_mev': str(case_v.REST_MEV),
        'w_start_mev': str(mp.w0 * case_p.A),
        'crest_deg': str(mp.crest_deg),
        'r_ext_m': str(mp.r_ext),
        'r_min_m': str(R_MIN),
        'b0_nominal_t': str(
            mc.nominal(case_p.MACHINE, case_p.ION, case_p.W_EXT * case_p.A)[0]
        ),
        'harmonic': str(case_p.HARMONIC),
    }
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist(out).replace_schema_metadata(meta), TABLE
    )


class Scaled:
    """A field multiplied by the constant s."""

    def __init__(self, field, s):
        self.field, self.s = field, s

    def __call__(self, x, y, z):
        """Return s times the field components at (x, y, z)."""
        return tuple(self.s * c for c in self.field(x, y, z))


def check(mapdir):
    """Check each rule of the reduced model against its reference."""
    t = reduced.Tables(TABLE)
    data = pq.read_table(TABLE).to_pydict()
    # Interpolation of the path length: drop every other energy, predict
    # it.
    lnp, lnl = np.log(data['p']), np.log(data['path_m'])
    from scipy.interpolate import CubicSpline

    miss = np.exp(CubicSpline(lnp[::2], lnl[::2])(lnp[1::2]) - lnl[1::2]) - 1
    print(
        'path length, every other energy held out: max relative error',
        float(np.max(np.abs(miss[:-1]))),
        '; scaled grid against the table, max',
        float(np.max(np.abs(t.path(lnp) / np.array(data['path_m']) - 1))),
    )
    # f_rev(p, s) = v(p) / L(p / s) against Boris closed orbits in a
    # scaled field.
    field = case_p.load_field(mapdir)
    for s in (0.997, 1.002):
        for w in (case_p.W_EXT / 18, case_p.W_EXT / 1.8):
            p = case_p.momentum(w)
            f = case_p.closed_orbit(
                Scaled(field, s),
                w,
                float(
                    np.exp(
                        np.interp(
                            np.log(p / s), lnp, np.log(data['r_section_m'])
                        )
                    )
                ),
            )[1]
            g = 1 + w * case_p.A / case_v.REST_MEV
            predicted = (p / g) * syn.C / t.path(np.log(p / s))
            print(
                f'scale {s}, w {w} MeV: f_rev Boris {f:.1f} Hz, '
                f'rule {predicted:.1f} Hz, '
                f'relative {predicted / f - 1:+.2e}',
                flush=True,
            )
    # Gain linear in the dee voltage: LIONS, one turn on the crest.
    mp = case_v.Map(mapdir, START_MEV, EXTRACTION)
    work = mp.workdir(mapdir / 'check_r')
    gains = [
        g
        for g in pq.read_table(mapdir / 'case_p_gain.parquet').to_pylist()
        if g['w_mev'] >= mp.w0
    ]
    orb = {
        r['w_mev']: r['r_co_m']
        for r in pq.read_table(mapdir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    }
    keep = case_p.DEE_KV
    for g in (gains[1], gains[min(3, len(gains) - 1)]):
        out = []
        for kv in (40.0, 60.0):
            case_p.DEE_KV = kv
            out.append(
                case_p.accelerate(
                    work,
                    g['w_mev'],
                    orb[g['w_mev']],
                    0.0,
                    g['crest_phi0_deg'],
                    2,
                    one_turn=True,
                )
            )
        print(
            f"gain at {g['w_mev']:.3f} MeV: 40 kV {out[0]:.5f}, "
            f'60 kV {out[1]:.5f} MeV; '
            f'ratio {out[1] / out[0]:.4f} against 1.5',
            flush=True,
        )
    case_p.DEE_KV = keep
    # Betatron action against the off-axis tracks of case P. The tracks
    # start at the start energy of the gate, not of this table.
    orbit_rows = {
        r['w_mev']: r['r_co_m']
        for r in pq.read_table(mapdir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    }
    p0 = case_p.momentum(
        case_p.start_energy(
            pq.read_table(mapdir / 'case_p_gain.parquet').to_pylist(),
            orbit_rows,
        )
    )
    for tag, z0 in (('offaxis_1', 3.0), ('offaxis_2', 6.0)):
        path = mapdir / f'case_p_slip_{tag}.parquet'
        rows = pq.read_table(path).to_pylist()
        w, z = (
            np.array([r['w_mev'] for r in rows]),
            np.abs([r['z_m'] for r in rows]) * 1e3,
        )
        for lo, hi in (
            np.array(((0.15, 0.5), (0.5, 2), (2, 6), (6, 12), (12, 18)))
            * case_p.W_EXT
            / 18
        ):
            on = (w >= lo) & (w < hi)
            p = case_p.momentum(w[on].mean())

            def j(q):
                return (
                    np.interp(np.log(q), t.qz_lnp, t.qz)
                    * q
                    / t.path(np.log(q))
                )  # noqa: E731

            print(
                f'{tag} {lo} to {hi} MeV/u: '
                f'envelope of |z| {z[on].max():.3f} mm, '
                f'action rule {z0 * np.sqrt(j(p0) / j(p)):.3f} mm',
                flush=True,
            )


def errors_of(curves_lions, phase_model, probe):
    """Return the count-weighted rms phase error, deg, against LIONS."""
    counts, phase = np.array(curves_lions[0]), np.array(curves_lions[1])
    on = (counts > 0) & (probe > 0)
    return float(
        np.sqrt(
            np.average((phase_model[on] - phase[on]) ** 2, weights=counts[on])
        )
    )


def compare(verdict, outdir):
    """Run reduced model and T0 on the verdict tunes, against LIONS."""
    model = reduced.Reduced(TABLE, PHASE_RMS_DEG)
    rows = pq.read_table(verdict / 'case_v_tunes.parquet').to_pylist()
    parts = pq.read_table(verdict / 'case_v_particles.parquet').to_pylist()
    curves = pq.read_table(verdict / 'case_v_curves.parquet').to_pylist()
    out, by_particle, times = [], [], []
    for row in rows:
        if row.get('harmonic', case_p.HARMONIC) != case_p.HARMONIC:
            continue
        m = case_v.of_tune(row)
        lions = sorted(
            (p for p in parts if p['tune'] == row['tune']),
            key=lambda p: p['particle'],
        )
        extracted = np.array([p['fate'] == case_v.EXTRACTED for p in lions])
        w_ext = np.array([p['w_extraction_mev'] for p in lions])[extracted]
        none = not len(w_ext)
        c = sorted(
            (x for x in curves if x['tune'] == row['tune']),
            key=lambda x: x['bin'],
        )
        counts = (
            np.array([x['lions_probe'] for x in c]) * case_v.N * case_v.N
        )  # the tables hold counts / N^2
        phase_lions = np.array([x['lions_phase_deg'] for x in c])
        # Beam.energies holds what the slit keeps, and the LIONS columns
        # are before it: all three are read with the slit open.
        unslit = {**m.knobs, 'slit': 0.0}
        began = time.perf_counter()
        b = model.run(m.machine, m.ion, unslit, m.errors, m.draws)
        wall = time.perf_counter() - began
        times.append(wall)
        began = time.perf_counter()
        t0 = mc.run(m.machine, m.ion, unslit, m.errors, m.draws)
        wall_t0 = time.perf_counter() - began
        converged = {
            f: model.run(
                m.machine, m.ion, unslit, m.errors, m.draws, limit_factor=f
            )
            for f in (2, 4)
        }
        rec = {
            'tune': row['tune'],
            'method': row['method'],
            'condition': row['condition'],
            'wall_s': wall,
            'wall_t0_s': wall_t0,
            'turns_reduced': b.turns,
            'turns_t0': t0.turns,
            'lions_extracted': int(extracted.sum()),
            'lions_transmission': float(extracted.mean()),
            'lions_w_mev': float('nan') if none else float(w_ext.mean()),
            'lions_spread': float('nan')
            if none
            else float(w_ext.std() / w_ext.mean()),
        }
        for name, beam in (('reduced', b), ('t0', t0)):
            e = beam.energies
            rec |= {
                f'{name}_extracted': len(e),
                f'{name}_transmission': len(e) / case_v.N,
                f'{name}_w_mev': float(e.mean()) if len(e) else float('nan'),
                f'{name}_spread': float(e.std() / e.mean())
                if len(e) > 1
                else 0.0,
                f'{name}_phase_rms_deg': errors_of(
                    (counts, phase_lions), beam.phase, beam.probe
                ),
            }
        rec |= {
            f'reduced_extracted_{f}x': len(x.energies)
            for f, x in converged.items()
        }
        out.append(rec)
        print(
            row['tune'],
            row['method'],
            f"extracted LIONS {rec['lions_extracted']} "
            f"reduced {rec['reduced_extracted']} "
            f"T0 {rec['t0_extracted']}; "
            f"spread {rec['lions_spread']:.4f} "
            f"{rec['reduced_spread']:.4f} {rec['t0_spread']:.4f}; "
            f"phase rms {rec['reduced_phase_rms_deg']:.1f} "
            f"{rec['t0_phase_rms_deg']:.1f}; {wall * 1e3:.0f} ms",
            flush=True,
        )
        by_particle += particle_rows(model, m, lions, row)
    meta = {
        'case': (
            f'synthetic {case_p.MISSION}, reduced model against LIONS and T0'
        ),
        'label': 'SYNTHETIC',
        'registration': case_p.REGISTRATION,
        'start_phase_law': str(PHASE_RMS_DEG),
        'extraction': EXTRACTION,
        **provenance(TABLE),
        'machine': platform.platform(),
        'units': (
            'transmission 1, before the septum; w_mev MeV; spread 1; '
            'phase_rms_deg deg; wall_s s, one core'
        ),
    }
    outdir.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist(out).replace_schema_metadata(meta),
        outdir / 'case_r_tunes.parquet',
    )
    pq.write_table(
        pa.Table.from_pylist(by_particle).replace_schema_metadata(meta),
        outdir / 'case_r_particles.parquet',
    )
    report(out, by_particle)


def particle_rows(model, m, lions, row):
    """Fate of each macroparticle: the reduced model run with the same
    start, particle by particle through its own record.
    """
    fate, w_end, turns = model.fates(
        m.machine, m.ion, m.knobs, m.errors, m.draws
    )
    return [
        {
            'tune': row['tune'],
            'particle': p['particle'],
            'phi_start_deg': p['phi_start_deg'],
            'z0_m': p['z0_m'],
            'lions_fate': p['fate'],
            'lions_w': p['w_final_mev'],
            'lions_turns': p['turns'],
            'reduced_fate': str(fate[p['particle']]),
            'reduced_w': float(w_end[p['particle']]),
            'reduced_turns': int(turns[p['particle']]),
        }
        for p in lions
    ]


def lines_report(rows, title):
    """Median absolute error per line. A tune where LIONS extracts
    nothing has no energy or spread and is left out, counted. A model
    that extracts nothing where LIONS does counts as infinite error.
    """
    a = lambda k: np.array([r[k] for r in rows], dtype=float)  # noqa: E731
    beam = a('lions_extracted') > 0
    measured = {
        'energy': lambda n: np.abs(a(f'{n}_w_mev') / a('lions_w_mev') - 1),
        'spread': lambda n: np.abs(a(f'{n}_spread') - a('lions_spread')),
    }
    lines = {
        'transmission': lambda n: np.abs(
            a(f'{n}_transmission') - a('lions_transmission')
        ),
        **{
            k: lambda n, g=g: np.where(a(f'{n}_extracted') > 0, g(n), np.inf)[
                beam
            ]
            for k, g in measured.items()
        },
        'phase_deg': lambda n: np.where(
            np.isfinite(a(f'{n}_phase_rms_deg')),
            a(f'{n}_phase_rms_deg'),
            np.inf,
        )[beam],
    }
    print(
        f'{title}: {len(rows)} tunes, {int((~beam).sum())} with no '
        'LIONS beam, left out of energy, spread and phase'
    )
    outcome = []
    for key, f in lines.items():
        red, t0 = f('reduced'), f('t0')
        passed = np.median(red) < THRESHOLD[key] + (
            1e-12 if key == 'transmission' else 0
        ) and np.median(red) < np.median(t0)
        outcome.append(passed)
        print(
            f'  {key}: threshold {THRESHOLD[key]:.4g}; '
            f'reduced median {np.median(red):.4g} '
            f'[{red.min():.4g}, {red.max():.4g}]; '
            f'T0 median {np.median(t0):.4g} '
            f'[{t0.min():.4g}, {t0.max():.4g}]: '
            f'{"PASS" if passed else "FAIL"}'
        )
    print(f'  {sum(outcome)} of {len(outcome)} lines pass')


def report(rows, particles):
    """Print the pass time, transmission and per-particle fates."""
    a = lambda k: np.array([r[k] for r in rows], dtype=float)  # noqa: E731
    lines_report(rows, 'all tunes')
    lines_report([r for r in rows if r['method'] == 'coach'], 'coach')
    wall = a('wall_s')
    print(
        f"pass time: median {np.median(wall) * 1e3:.1f} ms, "
        f"max {wall.max() * 1e3:.1f} ms; "
        f"T0 median {np.median(a('wall_t0_s')) * 1e3:.1f} ms: "
        f"{'PASS' if wall.max() < THRESHOLD['pass_s'] else 'FAIL'}"
    )
    for q, label in (
        ('lions_transmission', 'LIONS'),
        ('reduced_transmission', 'reduced'),
        ('t0_transmission', 'T0'),
    ):
        print(
            f'{label} transmission median {np.median(a(q)):.3f} '
            f'[{a(q).min():.3f}, {a(q).max():.3f}]; '
            f'energy {np.nanmedian(a(q.split("_")[0] + "_w_mev")):.3f}; '
            f'spread {np.nanmedian(a(q.split("_")[0] + "_spread")):.4f}'
        )
    print(
        'turn limit, extracted summed over tunes: 1x',
        int(a('reduced_extracted').sum()),
        '2x',
        int(a('reduced_extracted_2x').sum()),
        '4x',
        int(a('reduced_extracted_4x').sum()),
        '; tunes changing 2x to 4x',
        int((a('reduced_extracted_2x') != a('reduced_extracted_4x')).sum()),
    )
    fate = np.array(
        [
            (
                p['lions_fate'] == case_v.EXTRACTED,
                p['reduced_fate'] == case_v.EXTRACTED,
            )
            for p in particles
        ]
    )
    print(
        'per particle: extracted in both',
        int((fate[:, 0] & fate[:, 1]).sum()),
        'LIONS only',
        int((fate[:, 0] & ~fate[:, 1]).sum()),
        'reduced only',
        int((~fate[:, 0] & fate[:, 1]).sum()),
        'in neither',
        int((~fate[:, 0] & ~fate[:, 1]).sum()),
        'of',
        len(fate),
    )
    phi = np.array([p['phi_start_deg'] for p in particles])
    for lo, hi in (
        (-180, -150),
        (-150, -120),
        (-120, -90),
        (-90, -60),
        (-60, 0),
        (0, 60),
        (60, 90),
        (90, 120),
        (120, 150),
        (150, 180),
    ):
        on = (phi >= lo) & (phi < hi)
        if on.any():
            print(
                f'  start phase {lo} to {hi}: n {int(on.sum())}, '
                f'extracted LIONS {fate[on, 0].mean():.2f}, '
                f'reduced {fate[on, 1].mean():.2f}'
            )
    both = fate[:, 0] & fate[:, 1]
    crest = both & (np.abs(phi) < 30)
    dt = np.array([p['reduced_turns'] - p['lions_turns'] for p in particles])
    print(
        'turns, reduced minus LIONS, extracted in both, '
        '|phi start| < 30: median',
        float(np.median(dt[crest])),
        'range',
        int(dt[crest].min()),
        int(dt[crest].max()),
    )
    print(
        'fates',
        {
            f: int(sum(p['reduced_fate'] == f for p in particles))
            for f in (
                case_v.EXTRACTED,
                case_v.Z_APERTURE,
                case_v.OUT_OF_MAP,
                case_v.TURN_LIMIT,
            )
        },
        'LIONS',
        {
            f: int(sum(p['lions_fate'] == f for p in particles))
            for f in (
                case_v.EXTRACTED,
                case_v.Z_APERTURE,
                case_v.OUT_OF_MAP,
                case_v.TURN_LIMIT,
            )
        },
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('step', choices=('tables', 'check', 'compare'))
    parser.add_argument(
        '--mapdir', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_p'
    )
    parser.add_argument(
        '--verdict', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_v'
    )
    parser.add_argument(
        '--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_r'
    )
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if args.step == 'tables':
        tables(args.mapdir.resolve(), args.workers)
    elif args.step == 'check':
        check(args.mapdir.resolve())
    else:
        compare(args.verdict.resolve(), args.outroot.resolve())
