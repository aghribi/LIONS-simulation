"""Case R: reduced model of stage B against LIONS and T0 on the F18 tunes, docs/training/PHYSICAL.md. Synthetic map: says nothing about a real machine.

usage: uv run --extra game --extra rl --group bench --group dev python benchmarks/case_r.py {tables,check,compare} [--mapdir DIR] [--verdict DIR] [--outroot DIR] [--workers N]
tables: LIONS Hermite closed orbits at the energies of case P, no RF, into the scenario table. check: the rules of the model, each against its reference.
compare: reduced model and T0 against LIONS on the 38 tunes, with pass times.
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

import case_p as P  # noqa: E402
import case_v as V  # noqa: E402
import synthetic as syn  # noqa: E402
from lions.game import machine as mc  # noqa: E402
from lions.game import reduced as R  # noqa: E402
from provenance import provenance  # noqa: E402

TABLE = ROOT / 'src' / 'lions' / 'game' / 'scenarios' / 'f18_synthetic_orbits.parquet'
R_MIN = 0.005       # m, end of the map grid, as the verdict
THRESHOLD = {'transmission': 1 / V.N, 'energy': 0.005, 'spread': 0.005, 'phase_deg': 5.06, 'pass_s': 0.15}


def orbit_job(job):
    workdir, w, guess = job
    return P.lions_check((workdir, w, guess))


def tables(mapdir, workers):
    mp = V.Map(mapdir)
    rec = [r for r in pq.read_table(mapdir / 'case_p_orbits.parquet').to_pylist() if r['code'] == 'boris_analytic']
    dirs = [mp.workdir(mapdir / f'table_{i}') for i in range(len(rec))]
    with concurrent.futures.ProcessPoolExecutor(workers) as pool:
        rows = list(pool.map(orbit_job, [(d, r['w_mev'], r['r_co_m']) for d, r in zip(dirs, rec)]))
    gains = [g for g in pq.read_table(mapdir / 'case_p_gain.parquet').to_pylist() if g['w_mev'] >= mp.w0]
    lnw, per_kv = np.log([g['w_mev'] for g in gains]), np.array([g['gain_amp_mev'] for g in gains]) / P.DEE_KV
    out = []
    for r, b in zip(rows, rec):
        w = r['w_mev']
        g = 1 + w / V.REST_MEV
        p = np.sqrt(g * g - 1)
        f = P.F0 * (1 + r['f_rel_error'])
        out.append({'w_mev': w, 'p': float(p), 'path_m': float(p / g * syn.C / f), 'f_hz': float(f), 'r_section_m': r['r_co_m'], 'qz': r['qz_frac'], 'qr': 1 + r['qr_frac'],
                    'gain_mev_per_kv': float(np.interp(np.log(w), lnw, per_kv)), 'gain_measured': bool(w >= mp.w0)})
        print({k: out[-1][k] for k in ('w_mev', 'f_hz', 'r_section_m', 'qz', 'qr')}, 'boris r', b['r_co_m'], flush=True)
    meta = {'case': 'synthetic F18, closed orbits of LIONS Hermite, no RF', 'label': 'SYNTHETIC', 'registration': 'docs/training/PHYSICAL.md', **provenance(*mp.inputs), 'machine': platform.platform(),
            'units': 'w_mev MeV, p beta gamma, path_m m, f_hz Hz, r_section_m m at theta = 0, gain_mev_per_kv MeV/kV peak per turn; below w_start the gain is held, not measured',
            'rest_mev': str(V.REST_MEV), 'w_start_mev': str(mp.w0), 'crest_deg': str(mp.crest_deg), 'r_ext_m': str(mp.r_ext), 'r_min_m': str(R_MIN),
            'b0_nominal_t': str(mc.nominal(P.MACHINE, P.ION, P.W_EXT)[0]), 'harmonic': str(P.HARMONIC)}
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(out).replace_schema_metadata(meta), TABLE)


class Scaled:
    def __init__(self, field, s):
        self.field, self.s = field, s

    def __call__(self, x, y, z):
        return tuple(self.s * c for c in self.field(x, y, z))


def check(mapdir):
    t = R.Tables(TABLE)
    data = pq.read_table(TABLE).to_pydict()
    # Interpolation of the path length: drop every other energy, predict it.
    lnp, lnl = np.log(data['p']), np.log(data['path_m'])
    from scipy.interpolate import CubicSpline
    miss = np.exp(CubicSpline(lnp[::2], lnl[::2])(lnp[1::2]) - lnl[1::2]) - 1
    print('path length, every other energy held out: max relative error', float(np.max(np.abs(miss[:-1]))), '; scaled grid against the table, max', float(np.max(np.abs(t.path(lnp) / np.array(data['path_m']) - 1))))
    # f_rev(p, s) = v(p) / L(p / s) against Boris closed orbits in a scaled field.
    field = P.load_field(mapdir)
    for s in (0.997, 1.002):
        for w in (1.0, 10.0):
            p = P.momentum(w)
            f = P.closed_orbit(Scaled(field, s), w, float(np.exp(np.interp(np.log(p / s), lnp, np.log(data['r_section_m'])))))[1]
            g = 1 + w / V.REST_MEV
            predicted = (p / g) * syn.C / t.path(np.log(p / s))
            print(f'scale {s}, w {w} MeV: f_rev Boris {f:.1f} Hz, rule {predicted:.1f} Hz, relative {predicted / f - 1:+.2e}', flush=True)
    # Gain linear in the dee voltage: LIONS, one turn on the crest.
    mp = V.Map(mapdir)
    work = mp.workdir(mapdir / 'check_r')
    gains = [g for g in pq.read_table(mapdir / 'case_p_gain.parquet').to_pylist() if g['w_mev'] >= mp.w0]
    orb = {r['w_mev']: r['r_co_m'] for r in pq.read_table(mapdir / 'case_p_orbits.parquet').to_pylist() if r['code'] == 'boris_analytic'}
    keep = P.DEE_KV
    for g in (gains[1], gains[3]):
        out = []
        for kv in (40.0, 60.0):
            P.DEE_KV = kv
            out.append(P.accelerate(work, g['w_mev'], orb[g['w_mev']], 0.0, g['crest_phi0_deg'], 2, one_turn=True))
        print(f"gain at {g['w_mev']:.3f} MeV: 40 kV {out[0]:.5f}, 60 kV {out[1]:.5f} MeV; ratio {out[1] / out[0]:.4f} against 1.5", flush=True)
    P.DEE_KV = keep
    # Betatron action against the off-axis tracks of case P.
    for tag, z0 in (('offaxis_1', 3.0), ('offaxis_2', 6.0)):
        path = mapdir / f'case_p_slip_{tag}.parquet'
        rows = pq.read_table(path).to_pylist()
        w, z = np.array([r['w_mev'] for r in rows]), np.abs([r['z_m'] for r in rows]) * 1e3
        for lo, hi in ((0.15, 0.5), (0.5, 2), (2, 6), (6, 12), (12, 18)):
            on = (w >= lo) & (w < hi)
            p = np.sqrt((1 + w[on].mean() / t.rest_mev)**2 - 1)
            p0 = np.sqrt((1 + t.w_start / t.rest_mev)**2 - 1)
            j = lambda q: np.interp(np.log(q), t.qz_lnp, t.qz) * q / t.path(np.log(q))  # noqa: E731
            print(f'{tag} {lo} to {hi} MeV: envelope of |z| {z[on].max():.3f} mm, action rule {z0 * np.sqrt(j(p0) / j(p)):.3f} mm', flush=True)


def errors_of(curves_lions, phase_model, probe):
    counts, phase = np.array(curves_lions[0]), np.array(curves_lions[1])
    on = (counts > 0) & (probe > 0)
    return float(np.sqrt(np.average((phase_model[on] - phase[on])**2, weights=counts[on])))


def compare(verdict, outdir):
    model = R.Reduced(TABLE)
    rows = pq.read_table(verdict / 'case_v_tunes.parquet').to_pylist()
    parts = pq.read_table(verdict / 'case_v_particles.parquet').to_pylist()
    curves = pq.read_table(verdict / 'case_v_curves.parquet').to_pylist()
    out, by_particle, times = [], [], []
    for row in rows:
        m = V.of_tune(row)
        lions = sorted((p for p in parts if p['tune'] == row['tune']), key=lambda p: p['particle'])
        extracted = np.array([p['fate'] == V.EXTRACTED for p in lions])
        w_ext = np.array([p['w_extraction_mev'] for p in lions])[extracted]
        c = sorted((x for x in curves if x['tune'] == row['tune']), key=lambda x: x['bin'])
        counts = np.array([x['lions_probe'] for x in c]) * V.N * V.N     # the tables hold counts / N^2
        phase_lions = np.array([x['lions_phase_deg'] for x in c])
        began = time.perf_counter()
        b = model.run(m.machine, m.ion, m.knobs, m.errors, m.draws)
        wall = time.perf_counter() - began
        times.append(wall)
        t0 = m.beam
        began = time.perf_counter()
        mc.run(m.machine, m.ion, m.knobs, m.errors, m.draws)
        wall_t0 = time.perf_counter() - began
        converged = {f: model.run(m.machine, m.ion, m.knobs, m.errors, m.draws, limit_factor=f) for f in (2, 4)}
        rec = {'tune': row['tune'], 'method': row['method'], 'condition': row['condition'], 'wall_s': wall, 'wall_t0_s': wall_t0, 'turns_reduced': b.turns, 'turns_t0': t0.turns,
               'lions_extracted': int(extracted.sum()), 'lions_transmission': float(extracted.mean()), 'lions_w_mev': float(w_ext.mean()), 'lions_spread': float(w_ext.std() / w_ext.mean())}
        for name, beam in (('reduced', b), ('t0', t0)):
            e = beam.energies
            rec |= {f'{name}_extracted': len(e), f'{name}_transmission': len(e) / V.N, f'{name}_w_mev': float(e.mean()) if len(e) else float('nan'),
                    f'{name}_spread': float(e.std() / e.mean()) if len(e) > 1 else 0.0, f'{name}_phase_rms_deg': errors_of((counts, phase_lions), beam.phase, beam.probe)}
        rec |= {f'reduced_extracted_{f}x': len(x.energies) for f, x in converged.items()}
        out.append(rec)
        print(row['tune'], row['method'], f"extracted LIONS {rec['lions_extracted']} reduced {rec['reduced_extracted']} T0 {rec['t0_extracted']}; spread {rec['lions_spread']:.4f} {rec['reduced_spread']:.4f} {rec['t0_spread']:.4f}; phase rms {rec['reduced_phase_rms_deg']:.1f} {rec['t0_phase_rms_deg']:.1f}; {wall * 1e3:.0f} ms", flush=True)
        by_particle += particle_rows(model, m, lions, row)
    meta = {'case': 'synthetic F18, reduced model against LIONS and T0', 'label': 'SYNTHETIC', 'registration': 'docs/training/PHYSICAL.md', **provenance(TABLE), 'machine': platform.platform(),
            'units': 'transmission 1, before the septum; w_mev MeV; spread 1; phase_rms_deg deg; wall_s s, one core'}
    outdir.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(out).replace_schema_metadata(meta), outdir / 'case_r_tunes.parquet')
    pq.write_table(pa.Table.from_pylist(by_particle).replace_schema_metadata(meta), outdir / 'case_r_particles.parquet')
    report(out, by_particle)


def particle_rows(model, m, lions, row):
    """Fate of each macroparticle: the reduced model run with the same start, particle by particle through its own record."""
    fate, w_end, turns = model.fates(m.machine, m.ion, m.knobs, m.errors, m.draws)
    return [{'tune': row['tune'], 'particle': p['particle'], 'phi_start_deg': p['phi_start_deg'], 'z0_m': p['z0_m'], 'lions_fate': p['fate'], 'lions_w': p['w_final_mev'], 'lions_turns': p['turns'],
             'reduced_fate': str(fate[p['particle']]), 'reduced_w': float(w_end[p['particle']]), 'reduced_turns': int(turns[p['particle']])} for p in lions]


def report(rows, particles):
    a = lambda k: np.array([r[k] for r in rows], dtype=float)  # noqa: E731
    lines = {'transmission': lambda n: np.abs(a(f'{n}_transmission') - a('lions_transmission')), 'energy': lambda n: np.abs(a(f'{n}_w_mev') / a('lions_w_mev') - 1),
             'spread': lambda n: np.abs(a(f'{n}_spread') - a('lions_spread')), 'phase_deg': lambda n: a(f'{n}_phase_rms_deg')}
    for key, f in lines.items():
        red, t0 = f('reduced'), f('t0')
        verdict = 'PASS' if np.median(red) < THRESHOLD[key] + (1e-12 if key == 'transmission' else 0) and np.median(red) < np.median(t0) else 'FAIL'
        print(f'{key}: threshold {THRESHOLD[key]:.4g}; reduced median {np.median(red):.4g} [{red.min():.4g}, {red.max():.4g}]; T0 median {np.median(t0):.4g} [{t0.min():.4g}, {t0.max():.4g}]: {verdict}')
    wall = a('wall_s')
    print(f"pass time: median {np.median(wall) * 1e3:.1f} ms, max {wall.max() * 1e3:.1f} ms; T0 median {np.median(a('wall_t0_s')) * 1e3:.1f} ms: {'PASS' if wall.max() < THRESHOLD['pass_s'] else 'FAIL'}")
    for q, label in (('lions_transmission', 'LIONS'), ('reduced_transmission', 'reduced'), ('t0_transmission', 'T0')):
        print(f'{label} transmission median {np.median(a(q)):.3f} [{a(q).min():.3f}, {a(q).max():.3f}]; energy {np.median(a(q.split("_")[0] + "_w_mev")):.3f}; spread {np.median(a(q.split("_")[0] + "_spread")):.4f}')
    print('turn limit, extracted summed over tunes: 1x', int(a('reduced_extracted').sum()), '2x', int(a('reduced_extracted_2x').sum()), '4x', int(a('reduced_extracted_4x').sum()), '; tunes changing 2x to 4x', int((a('reduced_extracted_2x') != a('reduced_extracted_4x')).sum()))
    fate = np.array([(p['lions_fate'] == V.EXTRACTED, p['reduced_fate'] == V.EXTRACTED) for p in particles])
    print('per particle: extracted in both', int((fate[:, 0] & fate[:, 1]).sum()), 'LIONS only', int((fate[:, 0] & ~fate[:, 1]).sum()), 'reduced only', int((~fate[:, 0] & fate[:, 1]).sum()), 'in neither', int((~fate[:, 0] & ~fate[:, 1]).sum()), 'of', len(fate))
    phi = np.array([p['phi_start_deg'] for p in particles])
    for lo, hi in ((-180, -150), (-150, -120), (-120, -90), (-90, -60), (-60, 0), (0, 60), (60, 90), (90, 120), (120, 150), (150, 180)):
        on = (phi >= lo) & (phi < hi)
        if on.any():
            print(f'  start phase {lo} to {hi}: n {int(on.sum())}, extracted LIONS {fate[on, 0].mean():.2f}, reduced {fate[on, 1].mean():.2f}')
    both = fate[:, 0] & fate[:, 1]
    crest = both & (np.abs(phi) < 30)
    dt = np.array([p['reduced_turns'] - p['lions_turns'] for p in particles])
    print('turns, reduced minus LIONS, extracted in both, |phi start| < 30: median', float(np.median(dt[crest])), 'range', int(dt[crest].min()), int(dt[crest].max()))
    print('fates', {f: int(sum(p['reduced_fate'] == f for p in particles)) for f in (V.EXTRACTED, V.Z_APERTURE, V.OUT_OF_MAP, V.TURN_LIMIT)}, 'LIONS', {f: int(sum(p['lions_fate'] == f for p in particles)) for f in (V.EXTRACTED, V.Z_APERTURE, V.OUT_OF_MAP, V.TURN_LIMIT)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('step', choices=('tables', 'check', 'compare'))
    parser.add_argument('--mapdir', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_p')
    parser.add_argument('--verdict', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_v')
    parser.add_argument('--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_r')
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if args.step == 'tables':
        tables(args.mapdir.resolve(), args.workers)
    elif args.step == 'check':
        check(args.mapdir.resolve())
    else:
        compare(args.verdict.resolve(), args.outroot.resolve())
