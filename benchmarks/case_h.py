"""Case H: field-model and integrator options for LIONS (DESIGN.md section 23), phase 1 in the harness.

usage: uv run --group bench python benchmarks/case_h.py {validate,matrix,cost} [--outroot DIR]
validate: F1 node reproduction and C1 continuity; RK4 order. matrix: closed orbits, tunes, det(M) for every map,
field model, integrator and step count. cost: wall time per particle-step and memory of each option.
"""
import argparse
import contextlib
import io
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
sys.path.insert(0, str(ROOT / 'tests'))

import case_g1  # noqa: E402
import case_g2  # noqa: E402
import fieldmodels  # noqa: E402
import lions  # noqa: E402
import lions.cli  # noqa: E402
import orbits  # noqa: E402
import synthetic as syn  # noqa: E402
from provenance import provenance  # noqa: E402

RUNS = ROOT / 'benchmarks' / 'runs'
G1_MAPS = {'3-sector': (RUNS / 'case_g1' / '3sector_master', RUNS / 'case_g1' / 'step1'),
           '10-sector': (RUNS / 'case_g1' / '10sector_master', RUNS / 'case_g1' / '10sector' / 'step1')}
PSI_STEP3 = RUNS / 'case_g2' / 'step3' / 'case_g2_step3.parquet'
SPACINGS = {'3-sector': (1, 3), '10-sector': (1, 3), 'psi': (1,)}   # master columns per map column
STEPS = (720, 1440, 2880)
ORDER_STEPS = (90, 180, 360)
INTEGRATORS = {'boris': orbits.BorisTracker, 'rk4': orbits.RK4Tracker}
COST_PARTICLES = (1000, 100000)


def g1_map(name):
    """grids, attrs, cases [(label, t, p, q_ref)], mass_kg, charge_c for a G1 master."""
    master, step1 = G1_MAPS[name]
    grids, attrs = case_g1.load_master(master / 'master.h5')
    ref = pq.read_table(step1 / 'case_g1_step1.parquet').to_pylist()
    q_rft = {(r['a'], r['metric']): r['value'] for r in ref if r['samples'] == max(case_g1.SAMPLES)}
    species = case_g1.lions_species(float(attrs['mass_ev']), float(attrs['charge_e']))
    cases = []
    for t in pq.read_table(master / 'truth.parquet').to_pylist():
        p = np.sqrt(t['gamma']**2 - 1)
        q_ref = [case_g2.fold(t['qr']) + q_rft[(t['a'], 'qr_error')], case_g2.fold(t['qz']) + q_rft[(t['a'], 'qz_error')]]
        cases.append((f"a={t['a']}", t, p, q_ref))
    mass_kg = case_g1.lc.AMU * (species[0] - case_g1.lc.ELECTRON_MASS_U)
    return grids, attrs, cases, mass_kg, species[2] * case_g1.lc.E_CHARGE


def psi_map():
    """As g1_map for the PSI Ring; references from RF-Track on the quintic spline (G2 step 3)."""
    grids, attrs, _, mass_kg, _, fixpo = case_g2.setup()
    rft = {(round(r['energy_mev'], 3), r['metric']): r['value']
           for r in pq.read_table(PSI_STEP3).to_pylist() if r['code'] == 'rftrack'}
    cases = []
    for energy, r_ic, pr_ic, qr_ref, qz_ref in [*fixpo[::10], fixpo[-1]]:
        e3 = round(energy, 3)
        p, t, q_fixpo = case_g2.reference(energy, r_ic, pr_ic, qr_ref, qz_ref)
        t = {'r_theta0_m': t['r_theta0_m'] * (1 + rft[(e3, 'r_eo_rel_error')]),
             'pr_over_p_theta0': t['pr_over_p_theta0'] + rft[(e3, 'pr_eo_error_over_p')],
             'f_rev_hz': t['f_rev_hz'] / (1 + rft[(e3, 'turn_time_rel_error')])}
        q_ref = [q_fixpo[0] + rft[(e3, 'qr_error')], q_fixpo[1] + rft[(e3, 'qz_error')]]
        cases.append((f'E={energy:.3f}', t, p, q_ref))
    return grids, attrs, cases, mass_kg, case_g1.lc.E_CHARGE


def load_map(name):
    return psi_map() if name == 'psi' else g1_map(name)


def field_models(grids, attrs, every):
    return {'F0_bilinear': case_g1.bilinear_field(grids, attrs, every),
            'F1_hermite': fieldmodels.hermite_field(*fieldmodels.node_arrays(grids, attrs, every)),
            'F2_cubic': case_g1.smooth_field(grids, attrs, every, 3),
            'F3_quintic': case_g1.smooth_field(grids, attrs, every, 5)}


def write(records, outdir, name, case, *inputs):
    outdir.mkdir(parents=True, exist_ok=True)
    meta = {'case': case, **provenance(*inputs), 'machine': platform.platform(), 'python': platform.python_version()}
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / f'case_h_{name}.parquet')


def validate(outdir):
    records = []
    # F1 on the PSI map (steepest gradients): node reproduction and C1 continuity across cell edges
    grids, attrs, *_ = psi_map()
    r, th, bz, dbr, dbt = fieldmodels.node_arrays(grids, attrs, 1)
    surface = fieldmodels.hermite_surface(r, th, bz, dbr, dbt)
    rr, tt = np.meshgrid(r[:-1], th[:-1], indexing='ij')
    got = surface(rr.ravel(), tt.ravel())
    for k, (name, exact) in enumerate((('bz', bz), ('dbz_dr', dbr), ('dbz_dtheta', dbt))):
        err = np.max(np.abs(got[k] - exact[:-1, :-1].ravel())) / np.max(np.abs(exact))
        records.append({'check': 'F1_nodes', 'quantity': name, 'value': float(err)})
    rng = np.random.default_rng(0)
    i = rng.integers(1, len(r) - 2, 2000)
    j = rng.integers(1, len(th) - 2, 2000)
    eps_r, eps_t = 1e-9 * (r[1] - r[0]), 1e-9 * (th[1] - th[0])
    t_mid = th[j] + rng.random(2000) * (th[1] - th[0])
    r_mid = r[i] + rng.random(2000) * (r[1] - r[0])
    for edge, (lo, hi) in (('r_edge', (surface(r[i] - eps_r, t_mid), surface(r[i] + eps_r, t_mid))),
                           ('theta_edge', (surface(r_mid, th[j] - eps_t), surface(r_mid, th[j] + eps_t)))):
        for k, (name, exact) in enumerate((('bz', bz), ('dbz_dr', dbr), ('dbz_dtheta', dbt))):
            jump = np.max(np.abs(hi[k] - lo[k])) / np.max(np.abs(exact))
            records.append({'check': f'F1_{edge}_jump', 'quantity': name, 'value': float(jump)})
    # RK4 order on the 3-sector quintic field, outer orbit (reference exact to 1e-13)
    grids, attrs, cases, mass_kg, charge_c = g1_map('3-sector')
    field = case_g1.smooth_field(grids, attrs, 1, 5)
    label, t, p, q_ref = cases[2]
    for n in ORDER_STEPS:
        values = case_g1.orbit_metrics(orbits.RK4Tracker(field, mass_kg, charge_c, p, t['f_rev_hz'], n).track,
                                       t, p, q_ref)
        records.append({'check': f'rk4_order_{label}', 'quantity': f'r_eo_rel_error_{n}',
                        'value': float(values['r_eo_rel_error'][0])})
    for rec in records:
        print(rec, file=sys.stderr)
    write(records, outdir, 'validate', 'Case H validation: F1 Hermite, RK4 order')


def matrix(outdir):
    for name in ('3-sector', '10-sector', 'psi'):
        grids, attrs, cases, mass_kg, charge_c = load_map(name)
        records = []
        for every in SPACINGS[name]:
            spacing = case_g1.EVERY[every]
            for model, field in field_models(grids, attrs, every).items():
                for integrator, cls in INTEGRATORS.items():
                    for n in STEPS:
                        for label, t, p, q_ref in cases:
                            base = {'map': name, 'spacing': spacing, 'field': model, 'integrator': integrator,
                                    'steps_per_turn': n, 'case': label}
                            t0 = time.time()
                            try:
                                values = case_g1.orbit_metrics(cls(field, mass_kg, charge_c, p, t['f_rev_hz'], n).track,
                                                               t, p, q_ref)
                            except Exception as exc:
                                records.append({**base, 'metric': 'failed', 'value': None,
                                                'unit': f'{type(exc).__name__}: {exc}'[:200]})
                                print(f'{base} failed: {exc}', file=sys.stderr)
                                continue
                            for k, (v, u) in values.items():
                                records.append({**base, 'metric': k, 'value': float(v), 'unit': u})
                            print(f'{name} {spacing} {model} {integrator} {n} {label} {time.time() - t0:.0f} s ' +
                                  ' '.join(f'{k}={v[0]:+.1e}' for k, v in values.items()), file=sys.stderr)
        write(records, outdir, f'matrix_{name}', f'Case H matrix, {name}')


def cost(outdir):
    """Wall time per particle-step on the PSI map at 1/3 deg (300 MeV protons spread over the map), and LIONS."""
    grids, attrs, cases, mass_kg, charge_c = psi_map()
    p = cases[11][2]
    rng = np.random.default_rng(1)
    records = []
    for n_part in COST_PARTICLES:
        r = rng.uniform(2.2, 4.2, n_part)
        t = rng.uniform(0, 2 * np.pi, n_part)
        pos = np.column_stack([r * np.cos(t), r * np.sin(t), rng.uniform(-5e-3, 5e-3, n_part)])
        mom = p * np.column_stack([-np.sin(t), np.cos(t), np.zeros(n_part)])
        reps = 20 if n_part <= 1000 else 2
        for model, field in field_models(grids, attrs, 1).items():
            for integrator, cls in INTEGRATORS.items():
                tr = cls(field, mass_kg, charge_c, p, case_g2.F_RF / case_g2.HARMONIC, 1440)

                def step():
                    # one step of each scheme: Boris, one field call and rotation; RK4, four rates
                    if integrator == 'boris':
                        tr.kick(pos + tr.c_dt * mom, mom, 1.0)
                    else:
                        for _ in range(4):
                            tr.rate(pos, mom)

                step()
                samples = []
                for _ in range(5):
                    t0 = time.perf_counter()
                    for _ in range(reps):
                        step()
                    samples.append((time.perf_counter() - t0) / reps / n_part)
                records.append({'code': f'{model}_{integrator}', 'particles': n_part, 'metric': 's_per_particle_step',
                                'value': float(np.median(samples))})
                print(f'{model} {integrator} {n_part}: {np.median(samples):.2e} s/particle-step', file=sys.stderr)
        records.append({'code': 'lions', 'particles': n_part, 'metric': 's_per_particle_step',
                        'value': lions_cost(grids, attrs, pos, mom, outdir)})
        print(f"lions {n_part}: {records[-1]['value']:.2e} s/particle-step", file=sys.stderr)
    for model, nbytes in memory(grids, attrs).items():
        records.append({'code': model, 'particles': None, 'metric': 'bytes', 'value': float(nbytes)})
    write(records, outdir, 'cost', 'Case H cost on the PSI map, 1/3 deg, M3 CPU')


def lions_cost(grids, attrs, pos, mom, outdir, steps=20):
    """LIONS itself: plain tracking loop over `steps` steps, per particle-step."""
    work = outdir / 'lions_cost'
    work.mkdir(parents=True, exist_ok=True)
    with contextlib.chdir(work):
        case_g1.write_lions_map(grids, attrs, 1)
        syn.write_deck('deck.dat', fhf=case_g2.F_RF / case_g2.HARMONIC, r0=3.0, w0=300.0, turns=1, npaspt=steps,
                       nbrion=len(pos))
        d = lions.State()
        d.run.swap_bytes = 'n'
        with contextlib.redirect_stdout(io.StringIO()):
            lions.load(d, 'deck.dat')
            b = d.bunch
            b.x[:], b.y[:], b.z[:] = pos.T
            b.px[:], b.py[:], b.pz[:] = mom.T
            b.gamma[:] = np.sqrt(1 + (mom**2).sum(axis=1))
            t0 = time.perf_counter()
            lions.cli.track(d)
            elapsed = time.perf_counter() - t0
    return float(elapsed / (steps + 1) / len(pos))


def memory(grids, attrs):
    """Bytes held by each field model at 1/3 deg: node arrays (F0, F1), spline knots and coefficients (F2, F3)."""
    _, _, bz, _, _ = fieldmodels.node_arrays(grids, attrs, 1)
    out = {'F0_bilinear': 3 * bz.nbytes, 'F1_hermite': 4 * bz.nbytes}
    r = attrs['ri_m'] + attrs['rstep_m'] * np.arange(attrs['nrad'])
    nth = grids['bz'].shape[1]
    jj = np.arange(-case_g1.PAD, nth + case_g1.PAD)   # the nodes smooth_field splines (every = 1)
    for model, k in (('F2_cubic', 3), ('F3_quintic', 5)):
        total = 0
        for g in ('bz', 'dbz_dr', 'dbz_dtheta'):
            s = case_g1.RectBivariateSpline(r, attrs['theta_step_rad'] * jj, grids[g][:, jj % nth], kx=k, ky=k)
            total += s.get_coeffs().nbytes + sum(t.nbytes for t in s.get_knots())
        out[model] = total
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('task', choices=('validate', 'matrix', 'cost'))
    parser.add_argument('--outroot', type=Path, default=RUNS / 'case_h')
    args = parser.parse_args()
    {'validate': validate, 'matrix': matrix, 'cost': cost}[args.task](args.outroot.resolve())
