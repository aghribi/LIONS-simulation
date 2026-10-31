"""Case G1: codes against FROM_ORBIT manufactured solutions (DESIGN.md section 21).

usage: uv run --group bench python benchmarks/case_g1.py MASTER_DIR [outdir]
MASTER_DIR holds master.h5 and truth.parquet from benchmarks/g1_maps.py.
Step 1 (this script so far): precision of the references, RF-Track on a quintic spline of the exact map.
"""
import argparse
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

import RF_Track as rft  # noqa: E402
import synthetic as syn  # noqa: E402
import orbits  # noqa: E402
from provenance import provenance  # noqa: E402

SAMPLES = [1440, 2880]
TOL = 1e-12
H = 1e-6            # m, finite-difference step of the Jacobian
PAD = 5             # theta nodes wrapped on each side so the periodic spline has no edge


def load_master(path):
    with h5py.File(path) as f:
        grids = {k: f[k][()] for k in ('bz', 'dbz_dr', 'dbz_dtheta')}
        attrs = dict(f.attrs)
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
    parser.add_argument('outdir', nargs='?', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_g1' / 'step1')
    args = parser.parse_args()
    main(args.master_dir.resolve(), args.outdir.resolve())
