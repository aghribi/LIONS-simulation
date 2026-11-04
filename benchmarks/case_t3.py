"""Mission tuning, third amendment of docs/training/RL_RUNS.md: score
and stop on the readings, box from the tunes of the coach and PPO,
BoTorch baselines, stay-in-spec.

usage: uv run --extra game --extra rl --extra bo --group bench python
benchmarks/case_t3.py STEP ... tunes MISSION --policies DIR: coach and
PPO with the old stop rule, one Parquet per mission. box: containment
and the widened box, from the tunes. classical MISSION METHOD BOX: 20
conditions, one Parquet. BOX is original or widened. rl MISSION SEED:
100k steps, then 20 conditions. timing: one pass per mission. collect:
every part into one table per step.
"""

import argparse
import json
import multiprocessing
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

from lions.game import tuners as t  # noqa: E402 (after sys.path)
from lions.game.reduced import install_from_env  # noqa: E402 (after sys.path)
from lions.game.env import mission_by_id  # noqa: E402 (after sys.path)
from lions.game.missions import BUDGET, STEP_S  # noqa: E402 (after sys.path)
from provenance import provenance  # noqa: E402 (after sys.path)

MODEL = install_from_env()
MISSIONS = ('EYE', 'DEEP', 'F18', 'NEON', 'ISOL', 'QUAL')
RL_SEEDS = range(5)
HOLD = 50
CLASSICAL = {
    'coach': t.coach,
    'nelder_mead': t.nelder_mead,
    'random': t.random_search,
    'bayesian': t.bayesian,
    'nsga2': t.nsga2,
    'cma_es': t.cma_es,
    'turbo': t.turbo,
    'constrained_bo': t.constrained_bo,
    'mobo': t.mobo,
}
NO_BOX = ('coach',)  # the box does not bind these methods


def write(records, path, step):
    """Write records to parquet with provenance metadata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        'case': f'mission tuning, third amendment, {step}',
        'registration': 'docs/training/RL_RUNS.md',
        'units': 'knobs as in scenarios/missions.yaml, elapsed_s in s',
        'model': MODEL or 'T0',
        **provenance(),
        'machine': platform.platform(),
        'python': platform.python_version(),
    }
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta), path
    )


def floats(terms):
    """Return the terms as floats, or None when empty."""
    return {k: float(v) for k, v in terms.items()} if terms else None


def row(mission, method, box, training_seed, condition, m, run, elapsed):
    """The run has stopped: the first measured read in spec, or the
    budget. Then the knobs stay and the machine drifts.
    """
    in_spec = run.reached is not None
    true_at_stop = bool(m.ready())
    knobs_at_stop, field_at_stop = json.dumps(m.knobs), float(m.errors.field)
    held = []
    if in_spec:
        for _ in range(HOLD):
            m.wait(STEP_S)
            m.refresh()
            held.append(bool(m.ready()))
    return {
        'mission': mission,
        'method': method,
        'box': box,
        'training_seed': training_seed,
        'condition': condition,
        'in_spec_measured': in_spec,
        'evals_to_spec': float(run.reached) if in_spec else float('nan'),
        'true_at_stop': true_at_stop and in_spec,
        'hold_fraction': float(np.mean(held)) if held else float('nan'),
        'held_all': bool(in_spec and all(held)),
        'score': float(run.best),
        'best_terms': json.dumps(
            {k: floats(v) for k, v in run.best_terms.items()}
        )
        if run.best_terms
        else '',
        'model_failures': run.model_failures,
        'final_knobs': json.dumps({k: float(m.knobs[k]) for k in run.keys}),
        'knobs_at_stop': knobs_at_stop,
        'field_error_at_stop': field_at_stop,
        'budget': BUDGET,
        'elapsed_s': elapsed,
    }


def load_box(outdir, mission):
    """{knob: (below, above)} for the knobs the rule widens."""
    rows = pq.read_table(outdir / 'case_t3_box.parquet').to_pylist()
    return {
        r['knob']: (r['below'], r['above'])
        for r in rows
        if r['mission'] == mission and r['widened']
    }


def classical_job(job):
    """Tune one (mission, method, box, seed) and return its row."""
    mission, method, box_name, box, seed = job
    import torch

    torch.set_num_threads(1)
    m = mission_by_id(mission)
    m.reset(seed)
    t0 = time.perf_counter()
    run = t.tune(m, CLASSICAL[method], seed, box=box)
    elapsed = time.perf_counter() - t0
    out = row(mission, method, box_name, -1, seed, m, run, elapsed)
    print(
        f'{mission} {method} {box_name} {seed}: {run.reached} '
        f'evals {run.evals} {elapsed:.0f} s',
        file=sys.stderr,
        flush=True,
    )
    return out


def tune_old_rule(mission, method, seed, **options):
    """Tune with the stop rule m.ready, return machine, start, run."""
    m = mission_by_id(mission)
    m.reset(seed)
    start = dict(m.knobs)
    run = t.tune(m, method, seed, stop=m.ready, **options)
    return m, start, run


def tunes(mission, policies, outdir):
    """In-spec tunes of the coach and of every policy, by the old stop
    rule, as offsets from their own start.
    """
    rows = []
    jobs = [
        ('coach', -1, t.coach, {}),
        *(
            (
                'ppo',
                s,
                t.policy,
                {'model': t.load_policy(policies / f'{mission}_ppo_{s}.zip')},
            )
            for s in RL_SEEDS
        ),
    ]
    for name, training_seed, method, options in jobs:
        for condition in t.HELD_OUT:
            m, start, run = tune_old_rule(
                mission, method, condition, **options
            )
            original = _original_box(mission, condition)
            for spec in m.specs:
                k = spec['key']
                rows.append(
                    {
                        'mission': mission,
                        'method': name,
                        'training_seed': training_seed,
                        'condition': condition,
                        'in_spec': run.reached is not None,
                        'evals': float(run.reached or 'nan'),
                        'knob': k,
                        'start': start[k],
                        'tune': m.knobs[k],
                        'lo': original[k][0],
                        'hi': original[k][1],
                    }
                )
    write(rows, outdir / 'tunes' / f'case_t3_tunes_{mission}.parquet', 'tunes')


def _original_box(mission, condition):
    m = mission_by_id(mission)
    m.reset(condition)
    run = t.Run(m)
    return {
        k: (float(lo), float(hi))
        for k, lo, hi in zip(run.keys, run.lo, run.hi)
    }


def box(outdir):
    """The rule of the amendment, applied per mission and knob on the
    in-spec tunes.
    """
    records = pq.read_table(outdir / 'tunes').to_pylist()
    out = []
    for mission in MISSIONS:
        specs = {s['key']: s for s in mission_by_id(mission).specs}
        for k, spec in specs.items():
            sel = [
                r
                for r in records
                if r['mission'] == mission and r['knob'] == k and r['in_spec']
            ]
            width = t.BOX_ERRORS * spec['start_error']
            d = [r['tune'] - r['start'] for r in sel]
            below, above = (
                max(width, -min(d, default=0.0)),
                max(width, max(d, default=0.0)),
            )
            outside = sum(not r['lo'] <= r['tune'] <= r['hi'] for r in sel)
            out.append(
                {
                    'mission': mission,
                    'knob': k,
                    'tunes': len(sel),
                    'tunes_outside': outside,
                    'width': width,
                    'below': below,
                    'above': above,
                    'widened': below > width or above > width,
                    'min_offset': min(d, default=float('nan')),
                    'max_offset': max(d, default=float('nan')),
                }
            )
            print(out[-1], flush=True)
    write(out, outdir / 'case_t3_box.parquet', 'box')


def classical(mission, method, box_name, outdir, workers):
    """Run a classical method on all held-out conditions."""
    chosen = load_box(outdir, mission) if box_name == 'widened' else None
    if box_name == 'widened' and (not chosen or method in NO_BOX):
        print(
            f'{mission} {method}: the box does not change, '
            'widened rows repeat the original ones',
            flush=True,
        )
        return
    jobs = [(mission, method, box_name, chosen, s) for s in t.HELD_OUT]
    with multiprocessing.get_context('spawn').Pool(workers) as pool:
        records = pool.map(classical_job, jobs, chunksize=1)
    write(
        records,
        outdir
        / 'classical'
        / f'case_t3_{mission}_{method}_{box_name}.parquet',
        'classical',
    )


def rl(mission, seed, outdir):
    """Train one PPO policy and tune the held-out conditions."""
    t0 = time.perf_counter()
    model = t.train_policy(
        mission, 'ppo', seed, outdir / 'policies' / f'{mission}_ppo_{seed}'
    )
    train_s = time.perf_counter() - t0
    records = []
    for condition in t.HELD_OUT:
        m = mission_by_id(mission)
        m.reset(condition)
        t1 = time.perf_counter()
        run = t.tune(m, t.policy, condition, model=model)
        records.append(
            {
                **row(
                    mission,
                    'ppo',
                    'none',
                    seed,
                    condition,
                    m,
                    run,
                    time.perf_counter() - t1,
                ),
                'train_s': train_s,
            }
        )
    print(
        f'{mission} ppo {seed}: '
        f'{sum(r["in_spec_measured"] for r in records)}/20, '
        f'train {train_s:.0f} s',
        flush=True,
    )
    write(
        records, outdir / 'rl' / f'case_t3_{mission}_ppo_{seed}.parquet', 'rl'
    )


def timing():
    """Print the time per pass of each mission at its start."""
    for mission in MISSIONS:
        m = mission_by_id(mission)
        m.reset(t.HELD_OUT[0])
        t0 = time.perf_counter()
        for _ in range(20):
            m.refresh()
        print(
            mission,
            f'{(time.perf_counter() - t0) / 20 * 1e3:.1f} ms '
            'per pass at the start',
            flush=True,
        )


def collect(outdir):
    """Merge the per-run parquet files of classical and rl."""
    for sub in ('classical', 'rl'):
        if (outdir / sub).exists():
            table = pq.read_table(outdir / sub)
            pq.write_table(
                table.replace_schema_metadata(
                    pq.read_schema(
                        next((outdir / sub).glob('*.parquet'))
                    ).metadata
                ),
                outdir / f'case_t3_{sub}.parquet',
            )
            print(sub, table.num_rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'step',
        choices=('tunes', 'box', 'classical', 'rl', 'timing', 'collect'),
    )
    parser.add_argument('args', nargs='*')
    parser.add_argument(
        '--policies',
        type=Path,
        default=ROOT / 'benchmarks' / 'runs' / 'case_t_drift',
    )
    parser.add_argument(
        '--outroot',
        type=Path,
        default=ROOT / 'benchmarks' / 'runs' / 'case_t3',
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=int(os.environ.get('SLURM_CPUS_PER_TASK', 1)),
    )
    options = parser.parse_args()
    outdir, a = options.outroot.resolve(), options.args
    if options.step == 'tunes':
        tunes(a[0], options.policies, outdir)
    elif options.step == 'box':
        box(outdir)
    elif options.step == 'classical':
        classical(a[0], a[1], a[2], outdir, options.workers)
    elif options.step == 'rl':
        rl(a[0], int(a[1]), outdir)
    elif options.step == 'timing':
        timing()
    else:
        collect(outdir)
