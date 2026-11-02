"""Tuning methods on a mission: one pass of the machine is one evaluation. The arcade and benchmarks/case_t.py run this same code.

Registered comparison in docs/training/RL_RUNS.md. Each method sees control-room readings or the score built from them, never the hidden errors.
"""
import copy
import json
import threading
import warnings
from importlib.util import find_spec
from pathlib import Path

import numpy as np

from lions.game.coach import Coach
from lions.game.missions import BUDGET, STEP_S

BOX_ERRORS = 3
BO_INITIAL, BO_CANDIDATES = 20, 2000
NSGA_POPULATION = 40
RL_STEPS = {'ppo': 100_000, 'sac': 30_000}
TRAIN_SEEDS = list(range(1000))
HELD_OUT = list(range(10000, 10020))
POLICIES = Path.home() / '.lions_arcade' / 'policies'
NO_BEAM = 1e12    # worse than any term with beam out


class Run:
    """The budget and the search box. A method calls it once per pass; it raises Done when in spec or out of budget."""

    class Done(Exception):
        pass

    def __init__(self, m, budget=BUDGET, on_pass=None):
        self.m, self.budget, self.on_pass = m, budget, on_pass
        self.evals, self.reached, self.score, self.best, self.best_knobs = 0, None, -1.0, -1.0, dict(m.knobs)
        self.keys = [s['key'] for s in m.specs]
        centre = np.array([m.knobs[k] for k in self.keys])
        width = np.array([BOX_ERRORS * s['start_error'] for s in m.specs])
        self.lo = np.maximum(centre - width, [s['min'] for s in m.specs])
        self.hi = np.minimum(centre + width, [s['max'] for s in m.specs])
        self.start = (centre - self.lo) / (self.hi - self.lo)

    def passed(self):
        """The knobs are set: time passes, the machine runs once."""
        m = self.m
        m.wait(STEP_S)
        m.refresh()
        self.evals += 1
        self.score = m.tune_score()
        if self.score > self.best:
            self.best, self.best_knobs = self.score, dict(m.knobs)
        if self.on_pass:
            self.on_pass(self)
        if m.ready():
            self.reached = self.evals
            raise Run.Done
        if self.evals >= self.budget:
            raise Run.Done

    def set(self, u):
        u = np.clip(u, 0, 1)
        self.m.knobs.update(zip(self.keys, (self.lo + (self.hi - self.lo) * u).tolist()))
        self.passed()

    def __call__(self, u):
        self.set(u)
        return -self.score

    def terms(self, u):
        """One objective per term of the client's spec. Without beam every term is worse than with beam, closer on the radial probe."""
        self.set(u)
        m = self.m
        return list(m.spec_terms().values()) if m.has_beam() else [NO_BEAM * (2 - m.reach())] * len(m.spec_terms())


def coach(run, rng):
    """The rules of lions.game.coach: one knob per pass, its fine step at 0.1."""
    rules = Coach(run.m)
    while True:
        advice = rules.advise()
        if advice.knob:
            spec = rules.keys[advice.knob]
            step = advice.direction * spec['rate'] * STEP_S * (0.1 if advice.fine else 1)
            run.m.knobs[advice.knob] = run.m.clip(spec, run.m.knobs[advice.knob] + step)
        run.passed()


def nelder_mead(run, rng):
    from scipy.optimize import minimize
    n = len(run.start)
    start = run.start
    while True:     # restarts from a random point until the budget is gone
        simplex = np.vstack([start, start + 0.1 * np.eye(n)])
        minimize(run, start, method='Nelder-Mead', options={'maxfev': run.budget, 'initial_simplex': simplex, 'xatol': 1e-4, 'fatol': 1e-5})
        start = rng.uniform(0, 1, n)


def random_search(run, rng):
    while True:
        run(rng.uniform(0, 1, len(run.start)))


def bayesian(run, rng):
    from scipy.stats import norm
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
    n = len(run.start)
    xs = [run.start, *rng.uniform(0, 1, (BO_INITIAL - 1, n))]
    ys = [run(x) for x in xs]
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
        ys.append(run(x))
        xs.append(x)


def nsga2(run, rng):
    """Multi-objective genetic algorithm: the spec terms are kept apart, no weights between them."""
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.problem import Problem

    class Spec(Problem):
        def _evaluate(self, x, out, *args, **kwargs):
            out['F'] = np.array([run.terms(u) for u in x])

    n = len(run.start)
    first = np.vstack([run.start, rng.uniform(0, 1, (NSGA_POPULATION - 1, n))])
    algorithm = NSGA2(pop_size=NSGA_POPULATION, sampling=first)
    algorithm.setup(Spec(n_var=n, n_obj=len(run.m.spec_terms()), xl=0.0, xu=1.0), seed=int(rng.integers(2**31)))
    while True:
        algorithm.next()


def policy(run, rng, model):
    """A trained policy, deterministic: one action per pass, every knob at once."""
    from lions.game.env import MissionEnv
    env = MissionEnv(run.m.mission_id)
    env.m = run.m
    while True:
        action, _ = model.predict(env.observe(), deterministic=True)
        env.act(action)
        run.passed()


def tune(m, method, seed, budget=BUDGET, on_pass=None, **options):
    """Runs a method on mission m until in spec or out of budget. Same seed, same search."""
    run = Run(m, budget, on_pass)
    try:
        method(run, np.random.default_rng(seed), **options)
    except Run.Done:
        pass
    return run


def policy_path(mission_id, algo='ppo', seed=0):
    return POLICIES / f'{mission_id}_{algo}_{seed}.zip'


def train_policy(mission_id, algo, seed, path, steps=None, progress=None):
    """SB3 defaults, MlpPolicy, training conditions apart from the held-out ones. Saves to path."""
    import torch
    from stable_baselines3 import PPO, SAC
    from stable_baselines3.common.callbacks import BaseCallback
    from stable_baselines3.common.monitor import Monitor
    from lions.game.env import MissionEnv
    torch.set_num_threads(1)

    class Report(BaseCallback):
        def _on_step(self):
            if progress and self.num_timesteps % 1000 == 0:
                progress(self.num_timesteps)
            return True

    model = {'ppo': PPO, 'sac': SAC}[algo]('MlpPolicy', Monitor(MissionEnv(mission_id, TRAIN_SEEDS)), seed=seed, device='cpu', verbose=0)
    model.learn(steps or RL_STEPS[algo], callback=Report())
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    model.save(path)
    return model


def load_policy(path):
    from stable_baselines3 import PPO, SAC
    return (SAC if '_sac_' in Path(path).name else PPO).load(path, device='cpu')


def policy_report(mission_id):
    """Held-out result saved beside a policy by lions.game.train, or None."""
    path = policy_path(mission_id).with_suffix('.json')
    return json.loads(path.read_text()) if path.exists() else None


class Method:
    def __init__(self, name, fn, needs, what):
        self.name, self.fn, self.needs, self.what = name, fn, needs, what

    def missing(self):
        return [module for module in self.needs if find_spec(module) is None]


METHODS = [
    Method('COACH', coach, (), 'Written rules, one knob per pass, from the probes and readings. Not the CIME operator procedure.'),
    Method('NELDER-MEAD', nelder_mead, ('scipy',), 'Simplex on the score, no model, no gradient. Restarts from a random point when it stalls.'),
    Method('BAYESIAN OPT.', bayesian, ('scipy', 'sklearn'), 'Gaussian process on the score, expected improvement. Random points first, then the most promising.'),
    Method('NSGA-II (MOGA)', nsga2, ('pymoo',), 'Genetic algorithm on each term of the spec, kept apart: no weights between them.'),
    Method('PPO (RL)', policy, ('stable_baselines3',), 'A policy trained on other conditions, 100k passes. Reads the control room, moves every knob at once.'),
    Method('RANDOM', random_search, (), 'Random settings in the search box. What any method must beat.'),
]


class Tuning:
    """A method on a private copy of the mission, in a thread owned by the caller. The screen replays the passes as they come."""

    def __init__(self, stage, method, seed):
        self.method, self.seed, self.stopping = method, seed, False
        self.passes, self.run, self.error = [], None, None
        self.m = copy.deepcopy(stage)
        self.thread = threading.Thread(target=self.work, daemon=True)
        self.thread.start()

    def work(self):
        options = {}
        try:
            if self.method.fn is policy:
                options['model'] = load_policy(policy_path(self.m.mission_id))
            self.run = tune(self.m, self.method.fn, self.seed, on_pass=self.record, **options)
        except Exception as error:      # shown on the screen: a tuner must not take the game down
            self.error = f'{type(error).__name__}: {error}'

    def record(self, run):
        if self.stopping:
            raise Run.Done
        self.passes.append({'knobs': dict(run.m.knobs), 'field': run.m.errors.field, 'score': run.score, 'best': run.best})

    def finished(self):
        return not self.thread.is_alive()

    def stop(self):
        self.stopping = True
        self.thread.join(timeout=5)
