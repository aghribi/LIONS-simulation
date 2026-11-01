"""Case G1: codes against FROM_ORBIT manufactured solutions (DESIGN.md section 21).

usage: uv run --group bench python benchmarks/case_g1.py MASTER_DIR [--step 1|2|3] [--outroot DIR]
MASTER_DIR holds master.h5 and truth.parquet from benchmarks/g1_maps.py.
Step 1: precision of the references, RF-Track on a quintic spline of the exact map.
Step 2 (--step 2): LIONS on the bilinear map at several azimuthal spacings.
Step 3 (--step 3): the LIONS Boris step on bilinear (validation against LIONS), bicubic and quintic fields.
"""
import argparse
import contextlib
import platform
import sys
from pathlib import Path

import h5py
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.interpolate import RectBivariateSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'benchmarks'))

import lions.constants as lc  # noqa: E402
import RF_Track as rft  # noqa: E402
import synthetic as syn  # noqa: E402
import orbits  # noqa: E402
from provenance import provenance  # noqa: E402

SAMPLES = [1440, 2880]
STEPS = [360, 720, 1440]
EVERY = {1: '1/3 deg', 3: '1 deg', 12: '4 deg'}   # master columns per map column (1/3 deg master)
TOL = 1e-12
H = 1e-6            # m, finite-difference step of the Jacobian
RESIDUAL_MAX = 1e-9   # m, closed-orbit Newton residual above which the run counts as failed
PAD = 5             # theta nodes wrapped on each side so the periodic spline has no edge


def load_master(path):
    with h5py.File(path) as f:
        grids = {k: f[k][()] for k in ('bz', 'dbz_dr', 'dbz_dtheta')}
        attrs = dict(f.attrs)
    if not np.isclose(attrs['theta_step_rad'], np.radians(1 / 3), rtol=1e-12):
        raise ValueError(f"{path}: master step {np.degrees(attrs['theta_step_rad'])} deg, the EVERY labels need 1/3 deg")
    return grids, attrs


def spline_field(grids, attrs, every=1):
    """b(x, y, z) in the LIONS frame and sign (-Bz of FROM_ORBIT), first order in z, from quintic splines.

    every > 1 builds the splines on a subset of the nodes (interpolation-error estimate).
    """
    n_sec = int(attrs['sector_count'])
    r = attrs['ri_m'] + attrs['rstep_m'] * np.arange(attrs['nrad'])
    dth = attrs['theta_step_rad']
    nth = grids['bz'].shape[1]
    th = dth * np.arange(-PAD * every, nth + PAD * every)
    idx = np.arange(-PAD * every, nth + PAD * every) % nth
    splines = [RectBivariateSpline(r[::every], th[::every], -g[::every][:, idx][:, ::every], kx=5, ky=5)
               for g in (grids['bz'], grids['dbz_dr'], grids['dbz_dtheta'])]
    sector = 2 * np.pi / n_sec

    def b(x, y, z):
        rr = np.hypot(x, y)
        t = np.arctan2(y, x)
        ts = t % sector
        bz, dr, dt = (float(s.ev(rr, ts)) for s in splines)
        br, bt = z * dr, z * dt / rr
        return br * np.cos(t) - bt * np.sin(t), br * np.sin(t) + bt * np.cos(t), bz

    return b, splines


def interpolation_error(grids, attrs):
    """Spline on every second node, at the skipped nodes, against the master (relative, max over the grid)."""
    _, coarse = spline_field(grids, attrs, every=2)
    r = attrs['ri_m'] + attrs['rstep_m'] * np.arange(attrs['nrad'])
    th = attrs['theta_step_rad'] * np.arange(grids['bz'].shape[1])
    ri, tj = r[1:-1:2], th[1::2]
    approx = coarse[0](ri, tj)
    exact = -grids['bz'][1:-1:2][:, 1::2]
    return float(np.max(np.abs(approx / exact - 1)))


def lions_map_arrays(grids, attrs, every):
    """LIONS map nodes: master columns every `every`, full turn with the endpoint, LIONS sign; (r, theta, arrays)."""
    n_sec, nth = int(attrs['sector_count']), grids['bz'].shape[1]
    if nth % every:
        raise ValueError(f'{nth} master columns per sector are not a multiple of {every}')
    ntet = n_sec * nth // every + 1
    cols = (np.arange(ntet) * every) % nth
    th = attrs['theta_step_rad'] * every * np.arange(ntet)
    r = attrs['ri_m'] + attrs['rstep_m'] * np.arange(attrs['nrad'])
    bz, dr, dt = (-grids[k][:, cols] for k in ('bz', 'dbz_dr', 'dbz_dtheta'))
    dbzx = np.cos(th) * dr - np.sin(th) / r[:, None] * dt
    dbzy = np.sin(th) * dr + np.cos(th) / r[:, None] * dt
    return r, th, {'bz': bz, 'dbzx': dbzx, 'dbzy': dbzy}


def write_lions_map(grids, attrs, every):
    """bmap.dat (and empty potentials) in the working directory; returns the number of map columns."""
    r, th, m = lions_map_arrays(grids, attrs, every)
    syn.write_pot_central('pot1.bin')
    syn.write_pot_outer('pot2.dat')
    syn.write_bmap('bmap.dat', bz=lambda *_: m['bz'], dbzx=lambda *_: m['dbzx'], dbzy=lambda *_: m['dbzy'],
                   rmin=r[0], rmax=r[-1], nr=len(r), ntet=len(th))
    return len(th)


def bilinear_field(grids, attrs, every):
    """The LIONS field model on the same nodes: bilinear Bz, z x bilinear dBz/dx and dBz/dy (d2Bz = 0)."""
    r, th, m = lions_map_arrays(grids, attrs, every)
    dr, dth = r[1] - r[0], th[1] - th[0]

    def field(x, y, z):
        rr = np.hypot(x, y)
        t = np.arctan2(y, x) % (2 * np.pi)
        i = ((rr - r[0]) / dr).astype(int)
        j = (t / dth).astype(int)
        a, b = rr - r[0] - i * dr, t - j * dth
        w = [(dr - a) * (dth - b), a * (dth - b), a * b, (dr - a) * b]
        cell = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
        interp = {k: sum(wk * g[c] for wk, c in zip(w, cell)) / (dr * dth) for k, g in m.items()}
        return z * interp['dbzx'], z * interp['dbzy'], interp['bz']

    return field


def smooth_field(grids, attrs, every, k):
    """Splines of degree k (periodic in theta) of Bz, dBz/dr, dBz/dtheta on the same nodes; vectorised b(x, y, z)."""
    n_sec, nth = int(attrs['sector_count']), grids['bz'].shape[1]
    r = attrs['ri_m'] + attrs['rstep_m'] * np.arange(attrs['nrad'])
    jj = np.arange(-PAD * every, nth + PAD * every, every)
    th = attrs['theta_step_rad'] * jj
    splines = [RectBivariateSpline(r, th, -g[:, jj % nth], kx=k, ky=k)
               for g in (grids['bz'], grids['dbz_dr'], grids['dbz_dtheta'])]
    sector = 2 * np.pi / n_sec

    def field(x, y, z):
        rr = np.hypot(x, y)
        t = np.arctan2(y, x)
        ts = t % sector
        bz, dr, dt = (s.ev(rr, ts) for s in splines)
        br, bt = z * dr, z * dt / rr
        return br * np.cos(t) - bt * np.sin(t), br * np.sin(t) + bt * np.cos(t), bz

    return field


def lions_species(mass_ev, charge):
    """(atomic mass, mass number, charge) whose LIONS ion mass AMU (A - m_e q) is mass_ev."""
    atomic_mass = mass_ev * lc.E_CHARGE / (lc.AMU * lc.C_LIGHT**2) + lc.ELECTRON_MASS_U * charge
    return atomic_mass, round(atomic_mass), charge


def step2(master_dir, outdir):
    grids, attrs = load_master(master_dir / 'master.h5')
    truth = pq.read_table(master_dir / 'truth.parquet').to_pylist()
    ref = pq.read_table(outdir.parent / 'step1' / 'case_g1_step1.parquet').to_pylist()
    q_rft = {(r['a'], r['metric']): r['value'] for r in ref if r['samples'] == max(SAMPLES)}
    species = lions_species(float(attrs['mass_ev']), float(attrs['charge_e']))
    records = []
    for every, label in EVERY.items():
        work = outdir / f'map_every{every}'
        work.mkdir(parents=True, exist_ok=True)
        with contextlib.chdir(work):
            ntet = write_lions_map(grids, attrs, every)
        for t in truth:
            p = np.sqrt(t['gamma']**2 - 1)
            r0, pr0 = t['r_theta0_m'], p * t['pr_over_p_theta0']
            q_true = [np.arccos(np.cos(2 * np.pi * q)) / (2 * np.pi) for q in (t['qr'], t['qz'])]
            q_ref = [q_true[0] + q_rft[(t['a'], 'qr_error')], q_true[1] + q_rft[(t['a'], 'qz_error')]]
            w = t['ek_mev'] / species[1]
            for n in STEPS:
                lt = orbits.LionsTracker(work, w, t['f_rev_hz'], n, species)
                try:
                    s, res = orbits.closed_orbit(lt.track, (r0, pr0))
                    if res > RESIDUAL_MAX:
                        raise RuntimeError(f'closed orbit not converged: residual {res:.1e} m')
                    m = orbits.jacobian(lt.track, s, H, H * p / r0)
                    qr, qz = orbits.linear_tunes(m)
                    tt, pos, mom = lt.track(np.array([s]), 2.2)
                    cross = orbits.crossings(tt, pos[:, 0], mom[:, 0])
                except Exception as exc:   # a map too coarse for the orbit is a result, recorded as such
                    records.append({'map': label, 'map_columns': ntet, 'a': t['a'], 'steps_per_turn': n,
                                    'metric': 'failed', 'value': None, 'unit': f'{type(exc).__name__}: {exc}'[:200]})
                    print(f"{label} a={t['a']} n={n} failed: {exc}", file=sys.stderr)
                    continue
                values = {'r_eo_rel_error': (s[0] / r0 - 1, '1'), 'pr_eo_error_over_p': ((s[1] - pr0) / p, '1'),
                          'turn_time_rel_error': ((cross[1, 0] - cross[0, 0]) * t['f_rev_hz'] - 1, '1'),
                          'qr_error': (qr - q_ref[0], '1'), 'qz_error': (qz - q_ref[1], '1'),
                          'det_error': (orbits.symplectic_error(m)[1], '1'), 'residual': (res, 'm')}
                for k, (v, u) in values.items():
                    records.append({'map': label, 'map_columns': ntet, 'a': t['a'], 'steps_per_turn': n,
                                    'metric': k, 'value': float(v), 'unit': u})
                print(f"{label} a={t['a']} n={n} " + ' '.join(f'{k}={v[0]:+.2e}' for k, v in values.items()),
                      file=sys.stderr)
    meta = {'case': 'G1 step 2: LIONS on the bilinear FROM_ORBIT map', 'master': str(master_dir),
            **{f'master_{k}': str(v) for k, v in attrs.items()}, 'species': repr(species), 'fd_step_m': str(H),
            'tune_reference': f'RF-Track step 1, {max(SAMPLES)} samples per turn',
            **provenance(master_dir / 'master.h5', master_dir / 'truth.parquet',
                         outdir.parent / 'step1' / 'case_g1_step1.parquet'),
            'machine': platform.platform(), 'python': platform.python_version()}
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / 'case_g1_step2.parquet')


def orbit_metrics(track, t, p, q_ref):
    """Closed orbit, tunes, isochronism of one configuration against the references (as in step 2)."""
    r0, pr0 = t['r_theta0_m'], p * t['pr_over_p_theta0']
    s, res = orbits.closed_orbit(track, (r0, pr0))
    if res > RESIDUAL_MAX:
        raise RuntimeError(f'closed orbit not converged: residual {res:.1e} m')
    m = orbits.jacobian(track, s, H, H * p / r0)
    qr, qz = orbits.linear_tunes(m)
    tt, pos, mom = track(np.array([s]), 2.2)
    cross = orbits.crossings(tt, pos[:, 0], mom[:, 0])
    return {'r_eo_rel_error': (s[0] / r0 - 1, '1'), 'pr_eo_error_over_p': ((s[1] - pr0) / p, '1'),
            'turn_time_rel_error': ((cross[1, 0] - cross[0, 0]) * t['f_rev_hz'] - 1, '1'),
            'qr_error': (qr - q_ref[0], '1'), 'qz_error': (qz - q_ref[1], '1'),
            'det_error': (orbits.symplectic_error(m)[1], '1'), 'residual': (res, 'm')}


def step3(master_dir, outdir):
    grids, attrs = load_master(master_dir / 'master.h5')
    truth = pq.read_table(master_dir / 'truth.parquet').to_pylist()
    ref = pq.read_table(outdir.parent / 'step1' / 'case_g1_step1.parquet').to_pylist()
    q_rft = {(r['a'], r['metric']): r['value'] for r in ref if r['samples'] == max(SAMPLES)}
    mass_ev, charge = float(attrs['mass_ev']), float(attrs['charge_e'])
    species = lions_species(mass_ev, charge)
    mass_kg = lc.AMU * (species[0] - lc.ELECTRON_MASS_U * species[2])   # the LIONS ion mass
    charge_c = species[2] * lc.E_CHARGE
    records = []

    # 3a: the numpy Boris step on the LIONS field model against LIONS itself, 1 deg map, 360 steps
    t = truth[1]
    p = np.sqrt(t['gamma']**2 - 1)
    start = np.array([[t['r_theta0_m'], p * t['pr_over_p_theta0'], 0.0, 0.0]])
    work = outdir / 'validation'
    work.mkdir(parents=True, exist_ok=True)
    with contextlib.chdir(work):
        write_lions_map(grids, attrs, 3)
    lions_turn = orbits.one_turn(orbits.LionsTracker(work, t['ek_mev'] / species[1], t['f_rev_hz'], 360, species).track,
                                 start)
    boris = orbits.BorisTracker(bilinear_field(grids, attrs, 3), mass_kg, charge_c, p, t['f_rev_hz'], 360)
    port_turn = orbits.one_turn(boris.track, start)
    diff = float(np.max(np.abs(port_turn - lions_turn)[0, :2] / np.array([t['r_theta0_m'], p])))
    records.append({'interpolation': 'bilinear', 'map': '1 deg', 'a': t['a'], 'steps_per_turn': 360,
                    'metric': 'validation_vs_lions', 'value': diff, 'unit': '1'})
    print(f'3a: numpy Boris on the LIONS field vs LIONS, one turn: {diff:.2e}', file=sys.stderr)
    if diff > 1e-12:
        raise RuntimeError(f'numpy Boris differs from LIONS by {diff:.2e}: step 3b not run')

    # 3b: same step, same nodes, smooth interpolations
    for interp, k in (('bicubic', 3), ('quintic', 5)):
        for every in (1, 3):
            field = smooth_field(grids, attrs, every, k)
            for t in truth:
                p = np.sqrt(t['gamma']**2 - 1)
                q_true = [np.arccos(np.cos(2 * np.pi * q)) / (2 * np.pi) for q in (t['qr'], t['qz'])]
                q_ref = [q_true[0] + q_rft[(t['a'], 'qr_error')], q_true[1] + q_rft[(t['a'], 'qz_error')]]
                for n in STEPS:
                    tracker = orbits.BorisTracker(field, mass_kg, charge_c, p, t['f_rev_hz'], n)
                    try:
                        values = orbit_metrics(tracker.track, t, p, q_ref)
                    except Exception as exc:
                        records.append({'interpolation': interp, 'map': EVERY[every], 'a': t['a'],
                                        'steps_per_turn': n, 'metric': 'failed', 'value': None,
                                        'unit': f'{type(exc).__name__}: {exc}'[:200]})
                        continue
                    for key, (v, u) in values.items():
                        records.append({'interpolation': interp, 'map': EVERY[every], 'a': t['a'],
                                        'steps_per_turn': n, 'metric': key, 'value': float(v), 'unit': u})
                    print(f"{interp} {EVERY[every]} a={t['a']} n={n} " +
                          ' '.join(f'{key}={v[0]:+.2e}' for key, v in values.items()), file=sys.stderr)
    meta = {'case': 'G1 step 3: LIONS Boris step (numpy) on bilinear, bicubic, quintic fields',
            'master': str(master_dir), 'species': repr(species), 'fd_step_m': str(H),
            'tune_reference': f'RF-Track step 1, {max(SAMPLES)} samples per turn',
            **provenance(master_dir / 'master.h5', master_dir / 'truth.parquet',
                         outdir.parent / 'step1' / 'case_g1_step1.parquet'),
            'machine': platform.platform(), 'python': platform.python_version()}
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / 'case_g1_step3.parquet')


def main(master_dir, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    grids, attrs = load_master(master_dir / 'master.h5')
    truth = pq.read_table(master_dir / 'truth.parquet').to_pylist()
    mass_ev, charge = float(attrs['mass_ev']), float(attrs['charge_e'])
    if charge != 1:
        raise ValueError(f'RFTrackTracker tracks charge 1 e, the example has {charge}')
    mass_kg = mass_ev * syn.QE / syn.C**2      # RF-Track gets mass_ev / 1e6 back exactly
    interp_err = interpolation_error(grids, attrs)
    print(f'spline on every second node: max relative error {interp_err:.2e}', file=sys.stderr)
    b, _ = spline_field(grids, attrs)
    half_size = attrs['ri_m'] + attrs['rstep_m'] * (attrs['nrad'] - 1)
    records = [{'a': None, 'samples': None, 'metric': 'spline_half_grid_max_rel_error', 'value': interp_err,
                'unit': '1'}]
    fold = orbits.linear_tunes   # trace convention: fractional tune folded into [0, 0.5]
    for t in truth:
        p = np.sqrt(t['gamma']**2 - 1)
        r0, pr0 = t['r_theta0_m'], p * t['pr_over_p_theta0']
        q_true = [np.arccos(np.cos(2 * np.pi * q)) / (2 * np.pi) for q in (t['qr'], t['qz'])]
        for n in SAMPLES:
            rt = orbits.RFTrackTracker(b, half_size, mass_kg, p, 1 / t['f_rev_hz'], n, TOL)
            s, res = orbits.closed_orbit(rt.track, (r0, pr0))
            m = orbits.jacobian(rt.track, s, H, H * p / r0)
            qr, qz = fold(m)
            tt, pos, mom = rt.track(np.array([s]), 2.2)
            cross = orbits.crossings(tt, pos[:, 0], mom[:, 0])
            values = {'r_eo_rel_error': (s[0] / r0 - 1, '1'), 'pr_eo_error_over_p': ((s[1] - pr0) / p, '1'),
                      'turn_time_rel_error': ((cross[1, 0] - cross[0, 0]) * t['f_rev_hz'] - 1, '1'),
                      'qr_error': (qr - q_true[0], '1'), 'qz_error': (qz - q_true[1], '1'),
                      'det_error': (orbits.symplectic_error(m)[1], '1'), 'residual': (res, 'm')}
            for k, (v, u) in values.items():
                records.append({'a': t['a'], 'samples': n, 'metric': k, 'value': float(v), 'unit': u})
            print(f"a={t['a']} n={n} " + ' '.join(f'{k}={v[0]:+.2e}' for k, v in values.items()), file=sys.stderr)
    meta = {'case': 'G1 step 1: RF-Track on a quintic spline of the exact FROM_ORBIT map',
            'master': str(master_dir), **{f'master_{k}': str(v) for k, v in attrs.items()},
            'rf_track': str(rft.RF_TRACK_VERSION), 'tolerance': str(TOL), 'fd_step_m': str(H),
            **provenance(master_dir / 'master.h5', master_dir / 'truth.parquet'),
            'machine': platform.platform(), 'python': platform.python_version()}
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / 'case_g1_step1.parquet')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('master_dir', type=Path)
    parser.add_argument('--step', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_g1')
    args = parser.parse_args()
    out = args.outroot.resolve() / f'step{args.step}'
    {1: main, 2: step2, 3: step3}[args.step](args.master_dir.resolve(), out)
