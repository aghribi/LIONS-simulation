"""Case D2: long-term canonical angular momentum in the B0 field
(DESIGN.md section 16).

usage: uv run --group bench python benchmarks/case_d.py [outdir] Writes
case_d.parquet. OPAL runs in Docker; outdir must be visible to Docker.
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

import RF_Track as rft  # noqa: E402, N813 (after sys.path; external name)
import synthetic as syn  # noqa: E402 (after sys.path)
from provenance import provenance  # noqa: E402 (after sys.path)
import orbits  # noqa: E402 (after sys.path)
from case_a import opal_version  # noqa: E402 (after sys.path)
from case_b import exact_b0, lions_state  # noqa: E402 (after sys.path)
from fields_b import B0, K, R0, b_b0, write_field  # noqa: E402 (after sys.path)

W = 1.0  # MeV
STEPS = 360
TURNS = 10_000
OFFSET = 1e-3  # m, radial offset from the closed orbit
TOL = 1e-12


def invariant(pos, mom, q_over_mc):
    """Canonical angular momentum in beta gamma m (DESIGN.md section
    16), plane z = 0.
    """
    x, y = pos[:, 0], pos[:, 1]
    r = np.hypot(x, y)
    return (
        x * mom[:, 1]
        - y * mom[:, 0]
        + q_over_mc * B0 * ((1 + K * R0) * r**2 / 2 - K * r**3 / 3)
    )


def metrics(t, pos, mom, q_over_mc):
    """Relative drift of the invariant and of |p|, and turns tracked."""
    lz = invariant(pos, mom, q_over_mc)
    p = np.linalg.norm(mom, axis=1)
    slope = np.polyfit(t - t[0], lz - lz[0], 1)[0]
    return {
        'l_secular_rel': (slope * (t[-1] - t[0]) / abs(lz[0]), '1'),
        'l_max_dev_rel': (np.max(np.abs(lz - lz[0])) / abs(lz[0]), '1'),
        'p_rel_change': (p[-1] / p[0] - 1, '1'),
        'turns_tracked': ((t[-1] - t[0]) / (t[1] - t[0]) / STEPS, '1'),
    }


def main(outdir):
    """Track 10 000 turns in each code; write case_d.parquet."""
    outdir.mkdir(parents=True, exist_ok=True)
    versions = {
        'lions': importlib.metadata.version('lions'),
        'opal': opal_version(),
        'rf_track': str(rft.RF_TRACK_VERSION),
    }
    records = []

    def add(code, integrator, values):
        for metric, (value, unit) in values.items():
            records.append(
                {
                    'code': code,
                    'version': versions[code],
                    'integrator': integrator,
                    'metric': metric,
                    'value': float(value),
                    'unit': unit,
                }
            )

    with contextlib.chdir(outdir):
        write_field('B0', 361)
        d = lions_state(W, STEPS)
        p, mass, charge = d.bunch.p[0], d.bunch.mass[0], d.bunch.charge[0]
        q_over_mc = charge / (mass * syn.C)
        f_rev = d.rf.frequency / d.rf.harmonic
        state = np.array(
            [exact_b0(p, mass, charge)['r_eo'][0] + OFFSET, 0.0, 0.0, 0.0]
        )

        t, pos, mom = orbits.LionsTracker(
            '.', W, syn.revolution_frequency(W, B0), STEPS
        ).track(state, TURNS)
        add('lions', 'boris', metrics(t, pos[:, 0], mom[:, 0], q_over_mc))
        print('lions done', file=sys.stderr)
        for integ in ('LF2', 'RK4'):
            ot = orbits.OpalTracker(
                outdir / 'opal', d.field, mass, p, f_rev, STEPS, integ
            )
            add('opal', integ, metrics(*ot.track_one(state, TURNS), q_over_mc))
            print(f'opal {integ} done', file=sys.stderr)
        rt = orbits.RFTrackTracker(
            b_b0, d.field.rmax, mass, p, 1 / f_rev, STEPS, TOL
        )
        add(
            'rf_track',
            'rkf45',
            metrics(*rt.track_one(state, TURNS), q_over_mc),
        )

    table = pa.Table.from_pylist(records).replace_schema_metadata(
        {
            'case': f'D2: B0 field Bz = {B0} (1 - {K} (r - {R0})) T, '
            f'{W} MeV proton, r_eo + {OFFSET} m, z = 0, '
            f'{TURNS} turns at {STEPS} steps per turn',
            'rf_track_tolerance': str(TOL),
            **provenance(),
            'machine': platform.platform(),
            'python': platform.python_version(),
            **{f'version_{k}': v for k, v in versions.items()},
        }
    )
    pq.write_table(table, outdir / 'case_d.parquet')
    return table


if __name__ == '__main__':
    out = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else ROOT / 'benchmarks' / 'runs' / 'case_d'
    )
    for row in main(out.resolve()).to_pylist():
        print(row)
