# /// script requires-python = "==3.12.*" dependencies =
# ["numpy==1.26.4", "scipy==1.13.1", "ndsplines==0.2.0",
# "sympy==1.13.3", "matplotlib==3.9.2",
#                 "optimparallel==0.1.2", "h5py==3.12.1",
#                 "pyarrow==18.1.0"]
# ///
"""Case G1 inputs: exact field map and references from a FROM_ORBIT
example (DESIGN.md section 21).

usage: MPLBACKEND=Agg uv run --script benchmarks/g1_maps.py EXAMPLE
OUTDIR [--ri --rstep --nrad --ntheta] EXAMPLE: 3-sector
(example/3-sector/3sector-ideal.py) or 10-sector
(example/10-sector/10sector-2h.py).

FROM_ORBIT (GPL-3.0) stays in its clone
(benchmarks/external/from_orbit): the published example file is run as
is, cut before its field-map and CYCLOPS section, so no FROM_ORBIT code
enters this repository. The map is the exact field of the orbit family
(FROM_ORBIT def_B: root finding for a(r, theta), then B = brho(a)
curv(a, theta) / r_inf), not the output of its write_field_map, whose
smoothing spline of a(r, theta) departs from the exact field by up to
1e-4. Derivatives by central differences of the exact field. Writes
OUTDIR/master.h5 (Bz, dBz/dr, dBz/dtheta on one sector, FROM_ORBIT sign
and units) and OUTDIR/truth.parquet (closed orbits, revolution
frequency, tunes).
"""

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import os
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.constants import c

ROOT = Path(__file__).resolve().parents[1]
CLONE = ROOT / 'benchmarks' / 'external' / 'from_orbit'
EXAMPLES = {
    '3-sector': 'example/3-sector/3sector-ideal.py',
    '10-sector': 'example/10-sector/10sector-2h.py',
}
# the example's own marker: everything after it writes its map and runs
# CYCLOPS
CUT = '#create a field map'
H_R, H_TH = 1e-6, 1e-6  # m, rad: central-difference steps on the exact field


def run_example(name):
    """Run the published example up to its field-map section; return its
    namespace.
    """
    path = CLONE / EXAMPLES[name]
    src = path.read_text()
    if CUT not in src:
        raise RuntimeError(
            f'{path}: marker {CUT!r} not found; the example changed'
        )
    src = src.split(CUT)[0].replace('import cyclops_wrapper as cyclops\n', '')
    sys.path.insert(0, str(CLONE))
    ns = {'__name__': 'from_orbit_example'}
    exec(compile(src, str(path), 'exec'), ns)
    return ns, hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    """Build the exact field map and references; write OUTDIR files."""
    parser = argparse.ArgumentParser()
    parser.add_argument('example', choices=EXAMPLES)
    parser.add_argument('outdir', type=Path)
    parser.add_argument(
        '--ri', type=float, default=0.40, help='first radius [m]'
    )
    parser.add_argument(
        '--rstep', type=float, default=0.01, help='radial step [m]'
    )
    parser.add_argument('--nrad', type=int, default=161)
    parser.add_argument(
        '--ntheta', type=int, default=360, help='azimuthal points per sector'
    )
    parser.add_argument(
        '--orbits',
        type=int,
        nargs='+',
        default=[1, 2, 3],
        help='indices into the example orbit grid '
        '(0 is the centre; the last bounds the field)',
    )
    args = parser.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(
        out
    )  # the example writes tunes.txt and plots in the working directory

    t0 = time.time()
    ns, example_sha = run_example(args.example)
    orb = ns['orb']
    n_sec, rinf, mass_ev, charge = (
        orb['N'],
        orb['r-inf'],
        orb['mass'],
        orb['charge-state'],
    )
    ns['def_B'](orb)
    print(f'orbits set up in {time.time() - t0:.0f} s', file=sys.stderr)

    t0 = time.time()
    r = args.ri + args.rstep * np.arange(args.nrad)
    th = 2 * np.pi / n_sec / args.ntheta * np.arange(args.ntheta)
    field = orb['B(r,th)']
    bz, dr, dth = (np.empty((args.nrad, args.ntheta)) for _ in range(3))
    log = io.StringIO()
    with contextlib.redirect_stdout(
        log
    ):  # FROM_ORBIT prints a line for each poorly converged root
        for i, ri in enumerate(r):
            for j, tj in enumerate(th):
                bz[i, j] = field(ri, tj)
                dr[i, j] = (field(ri + H_R, tj) - field(ri - H_R, tj)) / (
                    2 * H_R
                )
                dth[i, j] = (field(ri, tj + H_TH) - field(ri, tj - H_TH)) / (
                    2 * H_TH
                )
    issues = log.getvalue().count('possible issue')
    print(
        f'exact map {args.nrad} x {args.ntheta} in '
        f'{time.time() - t0:.0f} s, {issues} root warnings',
        file=sys.stderr,
    )
    if issues:
        raise RuntimeError(
            f'{issues} root-finding warnings: the grid leaves the orbit family'
        )

    commit = subprocess.run(
        ['git', '-C', str(CLONE), 'rev-parse', 'HEAD'],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    meta = {
        'from_orbit_commit': commit,
        'example': EXAMPLES[args.example],
        'example_sha256': example_sha,
        'versions': ';'.join(
            f'{p}={importlib.metadata.version(p)}'
            for p in ('numpy', 'scipy', 'ndsplines', 'sympy')
        ),
        'sector_count': str(n_sec),
        'r_inf_m': repr(rinf),
        'mass_ev': repr(mass_ev),
        'charge_e': repr(charge),
        'sign': 'FROM_ORBIT convention: Bz > 0; LIONS uses -Bz for the '
        'same orbits travelled counterclockwise',
    }
    with h5py.File('master.h5', 'w') as f:
        for name, data, unit in (
            ('bz', bz, 'T'),
            ('dbz_dr', dr, 'T/m'),
            ('dbz_dtheta', dth, 'T/rad'),
        ):
            f.create_dataset(name, data=data).attrs['units'] = unit
        f.attrs.update(
            {
                'axes': 'r, theta',
                'ri_m': args.ri,
                'rstep_m': args.rstep,
                'nrad': args.nrad,
                'ntheta_per_sector': args.ntheta,
                'theta_step_rad': float(th[1] - th[0]),
                'theta_endpoint': 'excluded (periodic)',
                'fd_step_m': H_R,
                'fd_step_rad': H_TH,
                **meta,
            }
        )

    # references per prescribed orbit
    nur, nuz = ns['tunes'](orb['grid-a'], orb)
    rows = []
    for k in args.orbits:
        a = orb['grid-a'][k]
        beta = orb['beta(a)'](a)
        gamma = 1 / np.sqrt(1 - beta**2)
        r0 = orb['r'](a, 0.0) * rinf
        drdth0 = orb['drdth'](a, 0.0) * rinf
        ths = np.linspace(0, 2 * np.pi / n_sec, 2001)
        length = n_sec * np.trapz(orb['dsdth'](a, ths), ths) * rinf
        r_orbit = orb['r'](a, ths) * rinf
        if r_orbit.min() < r[0] + 0.03 or r_orbit.max() > r[-1] - 0.03:
            raise RuntimeError(
                f'orbit a = {a} reaches {r_orbit.min():.3f} to '
                f'{r_orbit.max():.3f} m, map is '
                f'{r[0]:.3f} to {r[-1]:.3f} m '
                '(LIONS needs 2 cm inside the map)'
            )
        rows.append(
            {
                'a': a,
                'beta': beta,
                'gamma': gamma,
                'ek_mev': (gamma - 1) * mass_ev / 1e6,
                'r_theta0_m': r0,
                'pr_over_p_theta0': drdth0 / np.hypot(r0, drdth0),
                'r_min_m': r_orbit.min(),
                'r_max_m': r_orbit.max(),
                'f_rev_hz': beta * c / length,
                'f_rev_isochronous_hz': c / (2 * np.pi * rinf),
                'qr': float(nur[k]),
                'qz': float(nuz[k]),
            }
        )
    pq.write_table(
        pa.Table.from_pylist(rows).replace_schema_metadata(meta),
        'truth.parquet',
    )
    for row in rows:
        print(
            {
                k: round(v, 6) if isinstance(v, float) else v
                for k, v in row.items()
            }
        )


if __name__ == '__main__':
    main()
