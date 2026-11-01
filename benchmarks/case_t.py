"""Mission tuning: coach, classical baselines and RL at equal budget (docs/training/RL_RUNS.md).

usage: uv run --extra game --extra rl --group bench python benchmarks/case_t.py --step classical|rl [--missions F18 EYE] [--outroot DIR]
"""
import argparse
import multiprocessing
import platform
import sys
import warnings
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.optimize import minimize
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

from lions.game.coach import Coach  # noqa: E402
from lions.game.env import STEP_S, MissionEnv, mission_by_id  # noqa: E402
from provenance import provenance  # noqa: E402

BUDGET = 200
HELD_OUT = list(range(10000, 10020))
TRAIN_SEEDS = list(range(1000))
RL_SEEDS = range(5)
ALL = ('EYE', 'DEEP', 'F18', 'NEON', 'ISOL', 'QUAL')
RL_STEPS = {'ppo': 100_000, 'sac': 30_000}
BO_INITIAL, BO_CANDIDATES = 20, 2000
WORKERS = 6
BOX_ERRORS = 3


class Counter:
    """The budget: one call is one pass of the machine. Stops the optimiser when in spec or out of budget."""

    class Done(Exception):
        pass

    def __init__(self, m):
        self.m, self.evals, self.reached, self.keys = m, 0, None, [s['key'] for s in m.specs]
        centre = np.array([m.knobs[k] for k in self.keys])
        width = np.array([BOX_ERRORS * s['start_error'] for s in m.specs])
        self.lo = np.maximum(centre - width, [s['min'] for s in m.specs])
        self.hi = np.minimum(centre + width, [s['max'] for s in m.specs])
        self.start = (centre - self.lo) / (self.hi - self.lo)
        self.best = -1.0

    def __call__(self, u):
        u = np.clip(u, 0, 1)
        self.m.knobs.update(zip(self.keys, (self.lo + (self.hi - self.lo) * u).tolist()))
        self.m.refresh()
        self.evals += 1
        score = self.m.tune_score()
        self.best = max(self.best, score)
        if self.m.ready():
            self.reached = self.evals
            raise Counter.Done
        if self.evals >= BUDGET:
            raise Counter.Done
        return -score


def nelder_mead(f, rng):
    n = len(f.start)
    start = f.start
    while True:     # restarts from the start with a wider simplex until the budget is gone
        simplex = np.vstack([start, start + 0.1 * np.eye(n)])
        minimize(f, start, method='Nelder-Mead', options={'maxfev': BUDGET, 'initial_simplex': simplex, 'xatol': 1e-4, 'fatol': 1e-5})
        start = rng.uniform(0, 1, n)


def random_search(f, rng):
    while True:
        f(rng.uniform(0, 1, len(f.start)))


def bayesian(f, rng):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
    n = len(f.start)
    xs = [f.start, *rng.uniform(0, 1, (BO_INITIAL - 1, n))]
    ys = [f(x) for x in xs]
    while True:
        kernel = ConstantKernel() * Matern(nu=2.5, length_scale=np.full(n, 0.3)) + WhiteKernel(1e-4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            gp = GaussianProcessRegressor(kernel, normalize_y=True, n_restarts_optimizer=0).fit(np.array(xs), np.array(ys))
        candidates = rng.uniform(0, 1, (BO_CANDIDATES, n))
        mean, std = gp.predict(candidates, return_std=True)
        gain = min(ys) - mean
        z = gain / np.maximum(std, 1e-9)
        x = candidates[int(np.argmax(gain * norm.cdf(z) + std * norm.pdf(z)))]
        ys.append(f(x))
        xs.append(x)


def coach_run(m):
    """The coach through the environment's action: one knob, STEP_S of travel, its fine step at 0.1."""
    coach, evals = Coach(m), 0
    while evals < BUDGET:
        advice = coach.advise()
        if advice.knob:
            spec = coach.keys[advice.knob]
            step = advice.direction * spec['rate'] * STEP_S * (0.1 if advice.fine else 1)
            m.knobs[advice.knob] = m.clip(spec, m.knobs[advice.knob] + step)
        m.refresh()
        evals += 1
        if m.ready():
            return evals
    return None


def classical(job):
    mission, method, seed = job
    m = mission_by_id(mission)
    m.reset(seed)
    rng = np.random.default_rng(seed)
    if method == 'coach':
        reached, evals = coach_run(m), BUDGET
        best = m.tune_score()
    else:
        f = Counter(m)
        try:
            {'nelder_mead': nelder_mead, 'random': random_search, 'bayesian': bayesian}[method](f, rng)
        except Counter.Done:
            pass
        reached, evals, best = f.reached, f.evals, f.best
    print(f'{mission} {method} {seed}: {reached}', file=sys.stderr, flush=True)
    return {'mission': mission, 'method': method, 'training_seed': -1, 'condition': seed, 'in_spec': reached is not None,
            'evals_to_spec': float(reached) if reached else float('nan'), 'score': float(best), 'budget': BUDGET}


def rl(job):
    import torch
    from stable_baselines3 import PPO, SAC
    from stable_baselines3.common.monitor import Monitor
    torch.set_num_threads(1)
    mission, algo, seed, outdir = job
    cls = {'ppo': PPO, 'sac': SAC}[algo]
    model = cls('MlpPolicy', Monitor(MissionEnv(mission, TRAIN_SEEDS)), seed=seed, device='cpu', verbose=0)
    model.learn(RL_STEPS[algo])
    model.save(outdir / f'{mission}_{algo}_{seed}')
    rows = []
    env = MissionEnv(mission)
    for condition in HELD_OUT:
        env.seeds = [condition]
        obs, _ = env.reset(seed=condition)
        reached, done = None, False
        while not done:
            action, _ = model.predict(obs, deterministic=True)
            obs, _, terminated, truncated, info = env.step(action)
            done = terminated or truncated
            if terminated:
                reached = env.steps
        rows.append({'mission': mission, 'method': algo, 'training_seed': seed, 'condition': condition, 'in_spec': reached is not None,
                     'evals_to_spec': float(reached) if reached else float('nan'), 'score': float(info['score']), 'budget': BUDGET})
    print(f'{mission} {algo} {seed}: {sum(r["in_spec"] for r in rows)}/{len(rows)}', file=sys.stderr, flush=True)
    return rows


def write(records, path, step):
    meta = {'case': f'mission tuning, {step}', 'registration': 'docs/training/RL_RUNS.md', **provenance(),
            'machine': platform.platform(), 'python': platform.python_version()}
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--step', choices=('classical', 'rl'), required=True)
    parser.add_argument('--missions', nargs='+', default=None)
    parser.add_argument('--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_t')
    args = parser.parse_args()
    outdir = args.outroot.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    with multiprocessing.get_context('spawn').Pool(WORKERS) as pool:
        if args.step == 'classical':
            jobs = [(m, method, s) for m in (args.missions or ALL) for method in ('coach', 'nelder_mead', 'random', 'bayesian') for s in HELD_OUT]
            records = pool.map(classical, jobs, chunksize=1)
        else:
            jobs = [(m, a, s, outdir) for m in (args.missions or ('F18', 'EYE')) for a in RL_STEPS for s in RL_SEEDS]
            records = [r for part in pool.imap_unordered(rl, jobs) for r in part]
    write(records, outdir / f'case_t_{args.step}.parquet', args.step)
