"""Post-hoc diagnostic of case T, not the registered metric. After the
first read in spec: hold the knobs 50 passes with drift, read without
noise, replay the actions blind.

usage: uv run --extra game --extra rl --group bench python
benchmarks/case_t_check.py [--outroot DIR]
"""

import argparse
import multiprocessing
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

import case_t  # noqa: E402 (after sys.path)
from lions.game import tuners as t  # noqa: E402 (after sys.path)
from lions.game.env import MissionEnv, mission_by_id  # noqa: E402 (after sys.path)
from lions.game.missions import STEP_S  # noqa: E402 (after sys.path)

HOLD = 50
FRESH = list(range(20000, 20020))


def verdict(m):
    """Hold the knobs HOLD passes, then read the truth without noise."""
    held = 0
    for _ in range(HOLD):
        m.wait(STEP_S)
        m.refresh()
        held += bool(m.ready())
    noise, m.p['noise'] = m.p['noise'], 0.0
    m.refresh()
    truth = bool(m.ready())
    m.p['noise'] = noise
    return {'hold': held / HOLD, 'truth': truth}


def classical(job):
    """Tune with a classical method, then run the verdict."""
    mission, method, seed = job
    m = mission_by_id(mission)
    m.reset(seed)
    run = t.tune(m, case_t.CLASSICAL[method], seed)
    done = (
        verdict(m) if run.reached else {'hold': float('nan'), 'truth': False}
    )
    return {
        'mission': mission,
        'method': method,
        'training_seed': -1,
        'condition': seed,
        'set': 'held_out',
        'in_spec': run.reached is not None,
        'evals': float(run.reached or 'nan'),
        'open_loop': False,
        **done,
    }


def rl(job):
    """Tune with a trained policy, run the verdict, replay blind."""
    mission, seed, conditions, outdir = job
    model = t.load_policy(outdir / f'{mission}_ppo_{seed}.zip')
    rows, replay = [], None
    for c in conditions:
        m = mission_by_id(mission)
        m.reset(c)
        actions = []
        original = model.predict

        def predict(obs, deterministic=True):
            action, state = original(obs, deterministic=deterministic)
            actions.append(action.copy())
            return action, state

        run = t.tune(
            m,
            t.policy,
            c,
            model=type('Recorded', (), {'predict': staticmethod(predict)}),
        )
        done = (
            verdict(m)
            if run.reached
            else {'hold': float('nan'), 'truth': False}
        )
        if run.reached and replay is None:
            replay = actions
        rows.append(
            {
                'mission': mission,
                'method': 'ppo',
                'training_seed': seed,
                'condition': c,
                'set': 'fresh' if c >= 20000 else 'held_out',
                'in_spec': run.reached is not None,
                'evals': float(run.reached or 'nan'),
                **done,
            }
        )
    for r, c in zip(rows, conditions):
        r['open_loop'] = False
        if replay is not None:
            m = mission_by_id(mission)
            m.reset(c)
            env = MissionEnv(mission)
            env.m = m
            run = t.Run(m)
            try:
                for action in replay:
                    env.act(action)
                    run.passed()
            except t.Run.Done:
                r['open_loop'] = run.reached is not None
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--outroot',
        type=Path,
        default=ROOT / 'benchmarks' / 'runs' / 'case_t_drift',
    )
    args = parser.parse_args()
    outdir = args.outroot.resolve()
    with multiprocessing.get_context('spawn').Pool(case_t.WORKERS) as pool:
        rows = pool.map(
            classical,
            [
                (m, a, s)
                for m in case_t.ALL
                for a in case_t.CLASSICAL
                for s in t.HELD_OUT
            ],
            chunksize=1,
        )
        jobs = [
            (m, s, t.HELD_OUT + FRESH, outdir)
            for m in case_t.ALL
            for s in case_t.RL_SEEDS
        ]
        rows += [r for part in pool.map(rl, jobs, chunksize=1) for r in part]
    pq.write_table(pa.Table.from_pylist(rows), outdir / 'check.parquet')
