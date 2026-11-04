"""Mission tuning: coach, classical baselines and RL at equal budget
(docs/training/RL_RUNS.md). The methods are lions.game.tuners.

usage: uv run --extra game --extra rl --group bench python
benchmarks/case_t.py --step classical|rl [--missions F18 EYE] [--algos
ppo] [--methods cma_es] [--outroot DIR]
"""

import argparse
import multiprocessing
import platform
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

from lions.game import tuners as t  # noqa: E402
from lions.game.env import mission_by_id  # noqa: E402
from lions.game.missions import BUDGET  # noqa: E402
from provenance import provenance  # noqa: E402

RL_SEEDS = range(5)
ALL = ('EYE', 'DEEP', 'F18', 'NEON', 'ISOL', 'QUAL')
CLASSICAL = {
    'coach': t.coach,
    'nelder_mead': t.nelder_mead,
    'random': t.random_search,
    'bayesian': t.bayesian,
    'nsga2': t.nsga2,
    'cma_es': t.cma_es,
}
WORKERS = 6


def row(mission, method, training_seed, condition, run):
    return {
        'mission': mission,
        'method': method,
        'training_seed': training_seed,
        'condition': condition,
        'in_spec': run.reached is not None,
        'evals_to_spec': float(run.reached) if run.reached else float('nan'),
        'score': float(run.best),
        'budget': BUDGET,
    }


def classical(job):
    mission, method, seed = job
    m = mission_by_id(mission)
    m.reset(seed)
    run = t.tune(m, CLASSICAL[method], seed)
    print(
        f'{mission} {method} {seed}: {run.reached}',
        file=sys.stderr,
        flush=True,
    )
    return row(mission, method, -1, seed, run)


def rl(job):
    mission, algo, seed, outdir = job
    model = t.train_policy(
        mission, algo, seed, outdir / f'{mission}_{algo}_{seed}'
    )
    rows = []
    for condition in t.HELD_OUT:
        m = mission_by_id(mission)
        m.reset(condition)
        rows.append(
            row(
                mission,
                algo,
                seed,
                condition,
                t.tune(m, t.policy, condition, model=model),
            )
        )
    print(
        f'{mission} {algo} {seed}: {sum(r["in_spec"] for r in rows)}/{len(rows)}',
        file=sys.stderr,
        flush=True,
    )
    return rows


def write(records, path, step):
    meta = {
        'case': f'mission tuning, {step}',
        'registration': 'docs/training/RL_RUNS.md',
        **provenance(),
        'machine': platform.platform(),
        'python': platform.python_version(),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta), path
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--step', choices=('classical', 'rl'), required=True)
    parser.add_argument('--missions', nargs='+', default=ALL)
    parser.add_argument('--algos', nargs='+', default=('ppo',))
    parser.add_argument('--methods', nargs='+', default=tuple(CLASSICAL))
    parser.add_argument(
        '--outroot',
        type=Path,
        default=ROOT / 'benchmarks' / 'runs' / 'case_t_drift',
    )
    args = parser.parse_args()
    outdir = args.outroot.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = (
        ''
        if set(args.methods) == set(CLASSICAL)
        else '_' + '_'.join(args.methods)
    )
    with multiprocessing.get_context('spawn').Pool(WORKERS) as pool:
        if args.step == 'classical':
            jobs = [
                (m, method, s)
                for m in args.missions
                for method in args.methods
                for s in t.HELD_OUT
            ]
            records = pool.map(classical, jobs, chunksize=1)
        else:
            jobs = [
                (m, a, s, outdir)
                for m in args.missions
                for a in args.algos
                for s in RL_SEEDS
            ]
            records = [
                r for part in pool.imap_unordered(rl, jobs) for r in part
            ]
    write(records, outdir / f'case_t_{args.step}{suffix}.parquet', args.step)
