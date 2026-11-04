"""Case O: does radial motion explain the energy spread of the LIONS
verdict? docs/training/PHYSICAL.md. Synthetic map: says nothing about a
real machine.

usage: uv run --extra game --extra rl --group bench --group dev python
benchmarks/case_o.py {record,analyse} [--outroot DIR] [--mapdir DIR]
[--tunes FILE] [--reference FILE] [--workers N] [--start-mev W] record:
the tunes of case V through LIONS, one row per particle and crossing of
theta = 0. analyse: fits, predictions and the registered outcome.
"""

import argparse
import multiprocessing
import os
import pickle
import shutil
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

import case_v  # noqa: E402 (after sys.path)

WINDOW, MIN_CROSSINGS = (
    60,
    30,
)  # crossings fitted before the last one, and the fewest accepted
EARLY_SKIP = 5  # crossings left out before the early window
TUNES_Q = np.linspace(
    0.005,
    0.12,
    231,
    # section frequency of the radial oscillation, turns^-1: Qr - 1 is
    # 0.014 to 0.034 on this map
)
NEAR_CREST_DEG = 30.0
SUPPORTED, REFUTED = (
    0.8,
    0.5,
)  # fraction of the within-tune energy variance explained
PHASE_BINS = (-180, -120, -90, -60, -30, 0, 30, 60, 90, 180)


class Recorded(case_v.Pass):
    """Pass that keeps every section crossing of every particle."""

    def __init__(self, *args, **kwargs):
        self.history, self.z_prev = [], None
        super().__init__(*args, **kwargs)

    def step(self, d):
        """After each step, keep the vertical positions."""
        carry_on = super().step(d)
        self.z_prev = d.bunch.z.copy()
        return carry_on

    def crossing(self, c, d, xp, yp, gp, tp):
        """Record r, energy, phase and z at the theta = 0 crossing."""
        b = d.bunch
        f = -yp[c] / (b.y[c] - yp[c])
        when = tp + f * (d.run.t - tp)
        zp = b.z if self.z_prev is None else self.z_prev
        phi = case_v.t0_phase(d.rf.omega * when, self.mp.crest_deg)
        self.history.append(
            (
                c.copy(),
                self.turns[c] + 1,
                when,
                xp[c] + f * (b.x[c] - xp[c]),
                (gp[c] + f * (b.gamma[c] - gp[c]) - 1) * case_v.REST_MEV,
                np.degrees(phi),
                zp[c] + f * (b.z[c] - zp[c]),
            )
        )
        super().crossing(c, d, xp, yp, gp, tp)


def one(job):
    """Track one tune through LIONS and pickle its crossings."""
    row, mapdir, parts, start_mev = job
    path = parts / f"{row['tune']:03d}.pkl"
    if path.exists():
        return row['tune']
    mp = case_v.Map(mapdir, start_mev)
    work = mp.workdir(parts / f"work_{row['tune']:03d}")
    ref = case_v.reference_orbit(work, mp)
    m = case_v.of_tune(row)
    p = Recorded(work, mp, ref, m.knobs, m.errors, m.draws)
    h = [np.concatenate(column) for column in zip(*p.history)]
    b = p.beam()
    tune = {
        'tune': row['tune'],
        'method': row['method'],
        'condition': row['condition'],
        'field_scale': float(p.s),
        'limit_turns': p.limit,
        'w_mean_mev': float(b.w_mean),
        'w_spread': float(b.w_spread),
        'n_extracted': int((p.fate == case_v.EXTRACTED).sum()),
    }
    particles = [
        {
            'tune': row['tune'],
            'particle': i,
            'phi_start_deg': float(np.degrees(p.phi[i])),
            'z0_m': float(p.z0[i]),
            'fate': str(p.fate[i]),
            'turns': int(p.turns[i]),
            'w_extraction_mev': float(p.w_last[i])
            if p.fate[i] == case_v.EXTRACTED
            else float('nan'),
        }
        for i in range(case_v.N)
    ]
    with open(path, 'wb') as f:
        pickle.dump((tune, particles, h), f)
    shutil.rmtree(work)
    print(
        row['tune'],
        row['method'],
        row['condition'],
        'extracted',
        tune['n_extracted'],
        'spread',
        tune['w_spread'],
        'crossings',
        len(h[0]),
        flush=True,
    )
    return row['tune']


UNITS = (
    'r_m m, w_mev MeV, phase_deg deg (T0 phase, 0 at the crest), '
    'z_m m, t_s s, turn crossing count, x_m a_m m, q 1/turn, w_* MeV'
)


def record(outdir, mp, tunes, workers, start_mev):
    """Run all tunes in a pool and write the history parquet tables."""
    rows = pq.read_table(tunes).to_pylist()
    parts = outdir / 'parts'
    parts.mkdir(parents=True, exist_ok=True)
    with multiprocessing.get_context('spawn').Pool(workers) as pool:
        list(
            pool.imap_unordered(
                one, [(r, mp.dir, parts, start_mev) for r in rows]
            )
        )
    done = [
        pickle.loads((parts / f"{r['tune']:03d}.pkl").read_bytes())
        for r in rows
    ]
    tables = [
        pa.table(
            {
                'tune': np.full(len(h[0]), t['tune'], dtype=np.int16),
                'particle': h[0].astype(np.int16),
                'turn': h[1].astype(np.int32),
                't_s': h[2],
                'r_m': h[3],
                'w_mev': h[4],
                'phase_deg': h[5],
                'z_m': h[6],
            }
        )
        for t, _, h in done
    ]
    history = pa.concat_tables(tables)
    case_v.write([], outdir / 'empty.parquet', 'x', mp, UNITS)
    meta = pq.read_schema(outdir / 'empty.parquet').metadata
    (outdir / 'empty.parquet').unlink()
    pq.write_table(
        history.replace_schema_metadata(
            {**meta, b'case': b'synthetic F18, section crossings per particle'}
        ),
        outdir / 'case_o_history.parquet',
        compression='zstd',
    )
    case_v.write(
        [d[0] for d in done],
        outdir / 'case_o_tunes.parquet',
        'radial motion per tune',
        mp,
        UNITS,
    )
    case_v.write(
        [p for d in done for p in d[1]],
        outdir / 'case_o_particles.parquet',
        'radial motion per particle',
        mp,
        UNITS,
    )
    shutil.rmtree(parts)
    print('crossings', history.num_rows)


def momentum(w):
    """Return p / (m c) for kinetic energy w, MeV."""
    return np.sqrt((1 + w / case_v.REST_MEV) ** 2 - 1)


def closed_orbit(mp):
    """Return the closed-orbit radius, m, as a spline of momentum."""
    if case_v.case_p.A != 1:
        raise ValueError(
            'case_o: w_mev in the orbit table is per nucleon, momentum() '
            'takes total energy: A = 1 only'
        )
    rec = sorted(
        (r['w_mev'], r['r_co_m'])
        for r in pq.read_table(mp.dir / 'case_p_orbits.parquet').to_pylist()
        if r['code'] == 'boris_analytic'
    )
    w, r = np.array(rec).T
    p = momentum(w)
    held = CubicSpline(p[::2], r[::2])(p[1::2]) - r[1::2]
    print(
        'closed-orbit table: max error at held-out nodes above 1 MeV',
        float(np.max(np.abs(held[w[1::2] > 1]))),
        'm',
    )
    return CubicSpline(p, r)


def fit(n, x):
    """x = c + a cos(2 pi q n) + b sin(2 pi q n), q on a grid. Returns
    q, c, a, b and the rms residual.
    """
    best = None
    for q in TUNES_Q:
        basis = np.column_stack(
            [
                np.ones_like(n, dtype=float),
                np.cos(2 * np.pi * q * n),
                np.sin(2 * np.pi * q * n),
            ]
        )
        coef, *_ = np.linalg.lstsq(basis, x, rcond=None)
        rss = float(np.sum((x - basis @ coef) ** 2))
        if best is None or rss < best[0]:
            best = (rss, q, *coef)
    return best[1], best[2], best[3], best[4], np.sqrt(best[0] / len(x))


def analyse(outdir, mp, reference):
    """Fit the radial oscillation of each extracted particle."""
    orbit = closed_orbit(mp)
    grid = np.linspace(
        momentum(mp.w0), momentum(case_v.case_p.W_EXT) * 1.15, 40000
    )
    r_grid = orbit(grid)

    def w_of(r, s):
        return (
            np.sqrt(1 + (s * np.interp(r, r_grid, grid)) ** 2) - 1
        ) * case_v.REST_MEV

    tunes = {
        r['tune']: r
        for r in pq.read_table(outdir / 'case_o_tunes.parquet').to_pylist()
    }
    particles = pq.read_table(outdir / 'case_o_particles.parquet').to_pylist()
    h = pq.read_table(
        outdir / 'case_o_history.parquet',
        columns=['tune', 'particle', 'turn', 'r_m', 'w_mev'],
    )
    tune_id, particle_id, turn_n, radius, energy = (
        h[c].to_numpy() for c in h.column_names
    )
    order = np.lexsort((turn_n, particle_id, tune_id))
    tune_id, particle_id, turn_n, radius, energy = (
        tune_id[order],
        particle_id[order],
        turn_n[order],
        radius[order],
        energy[order],
    )
    edges = (
        np.flatnonzero((np.diff(tune_id) != 0) | (np.diff(particle_id) != 0))
        + 1
    )
    starts = np.r_[0, edges]
    ends = np.r_[edges, len(order)]
    groups = {
        (int(tune_id[a]), int(particle_id[a])): (a, b)
        for a, b in zip(starts, ends)
    }
    was = {r['tune']: r for r in pq.read_table(reference).to_pylist()}
    print(
        'against the verdict table: max |w_spread diff|',
        max(
            abs(t['w_spread'] - was[k]['lions_w_spread'])
            for k, t in tunes.items()
        ),
        'max |w_mean diff| MeV',
        max(
            abs(t['w_mean_mev'] - was[k]['lions_w_mean'])
            for k, t in tunes.items()
        ),
    )
    fits = []
    for q in particles:
        if q['fate'] != case_v.EXTRACTED:
            continue
        first, last = groups[(q['tune'], q['particle'])]
        s = tunes[q['tune']]['field_scale']
        n, r, w = (
            turn_n[first:last].astype(float),
            radius[first:last],
            energy[first:last],
        )
        x = r - orbit(momentum(w) / s)
        out = {
            **q,
            'n_crossings': len(n),
            'r_last_m': r[-1],
            'r_prev_m': r[-2] if len(r) > 1 else float('nan'),
            'x_last_m': x[-1],
            'gain_last_mev': w[-1] - w[-2] if len(w) > 1 else float('nan'),
        }
        if len(n) > MIN_CROSSINGS:
            k = max(len(n) - 1 - WINDOW, 0)
            freq, c, a, b, res = fit(n[k:-1], x[k:-1])
            xhat = (
                c
                + a * np.cos(2 * np.pi * freq * n[-1])
                + b * np.sin(2 * np.pi * freq * n[-1])
            )
            amp = float(np.hypot(a, b))
            out.update(
                q_late=freq,
                c_m=c,
                a_m=amp,
                residual_m=res,
                x_pred_m=xhat,
                w_pred_mev=float(w_of(mp.r_ext - xhat, s)),
                w_offset_mev=float(w_of(mp.r_ext - c, s)),
                w_noosc_mev=float(
                    w_of(mp.r_ext, s) + np.mean(np.diff(w[k:])) / 2
                ),
                w_hi_mev=float(w_of(mp.r_ext + amp, s)),
                w_lo_mev=float(w_of(mp.r_ext - amp, s)),
                p_late=float(momentum(w[k:-1].mean())),
            )
            if len(n) > EARLY_SKIP + 2 * WINDOW:
                e = slice(EARLY_SKIP, EARLY_SKIP + WINDOW)
                fe, _, ae, be, rese = fit(n[e], x[e])
                out.update(
                    q_early=fe,
                    a_early_m=float(np.hypot(ae, be)),
                    residual_early_m=rese,
                    p_early=float(momentum(w[e].mean())),
                    w_early_mev=float(w[e].mean()),
                )
        fits.append(out)
    case_v.write(
        fits,
        outdir / 'case_o_fits.parquet',
        'radial oscillation fits per extracted particle',
        mp,
        UNITS,
    )
    report(fits, tunes, mp)


def explained(rows, key):
    """Pooled within-tune variance of the extraction energy removed by a
    predictor, offset per tune free, slope 1.
    """
    actual, pred = [], []
    by = {}
    for r in rows:
        by.setdefault(r['tune'], []).append(r)
    per_tune = []
    for s in by.values():
        if len(s) < 3:
            continue
        a = np.array([v['w_extraction_mev'] for v in s])
        p = np.array([v[key] for v in s])
        a, p = a - a.mean(), p - p.mean()
        actual.append(a)
        pred.append(p)
        per_tune.append(1 - np.sum((a - p) ** 2) / np.sum(a**2))
    a, p = np.concatenate(actual), np.concatenate(pred)
    return 1 - float(np.sum((a - p) ** 2) / np.sum(a**2)), float(
        np.median(per_tune)
    )


def spread(rows, key):
    """Return the median over tunes of std(key) / 18.0."""
    by = {}
    for r in rows:
        by.setdefault(r['tune'], []).append(r[key])
    return float(
        np.median([np.std(v) / 18.0 for v in by.values() if len(v) > 2])
    )


def report(fits, tunes, mp):
    """Print the fit statistics and the registered outcome."""
    usable = [
        f
        for f in fits
        if 'w_pred_mev' in f and np.isfinite(f.get('w_pred_mev', np.nan))
    ]
    print(
        'extracted',
        len(fits),
        'with a fit',
        len(usable),
        'median crossings',
        float(np.median([f['n_crossings'] for f in fits])),
    )
    print(
        'median |x_last - x_pred| mm',
        1e3
        * float(
            np.median([abs(f['x_last_m'] - f['x_pred_m']) for f in usable])
        ),
        'median fit residual mm',
        1e3 * float(np.median([f['residual_m'] for f in usable])),
    )
    for name, rows in (
        (
            'within 30 deg of the crest',
            [f for f in usable if abs(f['phi_start_deg']) <= NEAR_CREST_DEG],
        ),
        ('all extracted', usable),
    ):
        print(f'== {name}: {len(rows)} particles')
        for label, key in (
            ('oscillation fit', 'w_pred_mev'),
            ('constant offset, no oscillation', 'w_offset_mev'),
            ('no radial term, half a turn gain', 'w_noosc_mev'),
        ):
            pooled, median = explained(rows, key)
            print(
                f'  variance explained, {label}: pooled {pooled:.3f}, '
                f'median over tunes {median:.3f}; '
                f'predicted spread {spread(rows, key):.5f}'
            )
        print(
            "  actual spread, median over tunes "
            f"{spread(rows, 'w_extraction_mev'):.5f}"
        )
        pooled = explained(rows, 'w_pred_mev')[0]
        print(
            '  REGISTERED OUTCOME:',
            'supported'
            if pooled >= SUPPORTED
            else 'refuted'
            if pooled < REFUTED
            else 'partial',
            f'(pooled {pooled:.3f}; supported at {SUPPORTED}, '
            f'refuted below {REFUTED})',
        )
        inside = [
            f['w_lo_mev'] - 0.1 <= f['w_extraction_mev'] <= f['w_hi_mev'] + 0.1
            for f in rows
        ]
        print(
            '  extraction energy within w(r_ext -+ A), one turn gain '
            f'(0.1 MeV) of margin: {np.mean(inside):.3f}'
        )
        print(
            f"  A median mm {1e3 * np.median([f['a_m'] for f in rows]):.3f}, "
            "q median "
            f"{np.median([f['q_late'] for f in rows]):.4f}, "
            "c median mm "
            f"{1e3 * np.median([f['c_m'] for f in rows]):.3f}, "
            "c 5 to 95% mm "
            f"{1e3 * np.percentile([f['c_m'] for f in rows], 5):.3f} "
            f"{1e3 * np.percentile([f['c_m'] for f in rows], 95):.3f}"
        )
    print('== amplitude against start phase')
    for lo, hi in zip(PHASE_BINS[:-1], PHASE_BINS[1:]):
        s = [f for f in usable if lo <= f['phi_start_deg'] < hi]
        if s:
            early = [
                f['a_early_m'] * np.sqrt(f['p_early'] / f['p_late'])
                for f in s
                if 'a_early_m' in f
            ]
            print(
                f"  {lo:5d} to {hi:4d} deg: n {len(s):4d}; "
                "A late mm median "
                f"{1e3 * np.median([f['a_m'] for f in s]):.3f}, "
                "5 to 95% "
                f"{1e3 * np.percentile([f['a_m'] for f in s], 5):.3f} "
                f"{1e3 * np.percentile([f['a_m'] for f in s], 95):.3f}; "
                "early A carried by sqrt(p) mm "
                f"{1e3 * np.median(early) if early else float('nan'):.3f} "
                f"(n {len(early)}); energy spread "
                f"{np.std([f['w_extraction_mev'] for f in s]) / 18.0:.5f}"
            )
    near = [
        f
        for f in usable
        if abs(f['phi_start_deg']) <= NEAR_CREST_DEG and 'a_early_m' in f
    ]
    ratio = [
        (f['a_m'] ** 2 * f['p_late'])
        / max(f['a_early_m'] ** 2 * f['p_early'], 1e-30)
        for f in near
    ]
    print(
        'action ratio late / early, within 30 deg: median',
        float(np.median(ratio)),
        '5 to 95%',
        np.percentile(ratio, [5, 95]).tolist(),
        'n',
        len(ratio),
        'early window mean energy MeV median',
        float(np.median([f['w_early_mev'] for f in near])),
        'q early median',
        float(np.median([f['q_early'] for f in near])),
    )
    print('== turn separation and last radii')
    seps = {}
    for f in fits:
        if np.isfinite(f['r_prev_m']):
            seps.setdefault(f['tune'], []).append(f)
    med = [
        np.median([v['r_last_m'] - v['r_prev_m'] for v in s])
        for s in seps.values()
    ]
    print(
        'turn separation r_last - r_prev mm: median over tunes',
        1e3 * float(np.median(med)),
        'range',
        1e3 * float(np.min(med)),
        1e3 * float(np.max(med)),
    )
    all_sep = np.array(
        [
            f['r_last_m'] - f['r_prev_m']
            for f in fits
            if np.isfinite(f['r_prev_m'])
        ]
    )
    print(
        'all particles, mm: percentiles 1 5 50 95 99',
        (1e3 * np.percentile(all_sep, [1, 5, 50, 95, 99])).round(3).tolist(),
        'fraction with r_last <= r_prev',
        float(np.mean(all_sep <= 0)),
    )
    last = np.array([f['r_last_m'] - mp.r_ext for f in fits])
    prev = np.array(
        [f['r_prev_m'] - mp.r_ext for f in fits if np.isfinite(f['r_prev_m'])]
    )
    print(
        'r_last - r_ext mm: percentiles 1 5 50 95 99',
        (1e3 * np.percentile(last, [1, 5, 50, 95, 99])).round(3).tolist(),
    )
    print(
        'r_prev - r_ext mm: percentiles 1 5 50 95 99',
        (1e3 * np.percentile(prev, [1, 5, 50, 95, 99])).round(3).tolist(),
    )
    gains = np.array(
        [f['gain_last_mev'] for f in fits if np.isfinite(f['gain_last_mev'])]
    )
    print(
        'gain on the last turn MeV: percentiles 5 50 95',
        np.percentile(gains, [5, 50, 95]).round(4).tolist(),
        '; spread of tunes',
        float(np.median([t['w_spread'] for t in tunes.values()])),
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('step', choices=('record', 'analyse'))
    parser.add_argument(
        '--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_o'
    )
    parser.add_argument(
        '--mapdir', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_p'
    )
    parser.add_argument(
        '--tunes',
        type=Path,
        default=ROOT
        / 'benchmarks'
        / 'runs'
        / 'case_v'
        / 'case_v_tunes.parquet',
    )
    parser.add_argument(
        '--reference',
        type=Path,
        default=ROOT
        / 'benchmarks'
        / 'runs'
        / 'case_v'
        / 'case_v_verdict.parquet',
    )
    parser.add_argument('--start-mev', type=float)
    parser.add_argument(
        '--workers',
        type=int,
        default=int(os.environ.get('SLURM_CPUS_PER_TASK', 1)),
    )
    args = parser.parse_args()
    outdir = args.outroot.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    mp_ = case_v.Map(args.mapdir, args.start_mev)
    if args.step == 'record':
        record(outdir, mp_, args.tunes.resolve(), args.workers, args.start_mev)
    else:
        analyse(outdir, mp_, args.reference.resolve())
