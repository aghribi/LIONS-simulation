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
CMA_SIGMA = 0.3    # of the search box
RL_STEPS = {'ppo': 100_000, 'sac': 30_000}
TRAIN_SEEDS = list(range(1000))
HELD_OUT = list(range(10000, 10020))
POLICIES = Path.home() / '.lions_arcade' / 'policies'
NO_BEAM = 1e12    # worse than any term with beam out
TERM_CLAMP = 1e3  # per-term models: log(term) within +- log of this
ACQ_RESTARTS, ACQ_RAW_SAMPLES = 4, 256


class Run:
    """The budget and the search box. A method calls it once per pass; it raises Done when in spec or out of budget."""

    class Done(Exception):
        pass

    def __init__(self, m, budget=BUDGET, on_pass=None, box=None, stop=None):
        """box: {knob: (below, above)}, offsets from the start replacing BOX_ERRORS. stop: the end test, the measured one by default."""
        self.m, self.budget, self.on_pass = m, budget, on_pass
        self.stop = stop or m.measured_ready
        self.evals, self.reached, self.score, self.best, self.best_knobs = 0, None, -1.0, -1.0, dict(m.knobs)
        self.best_terms, self.model_failures = None, 0
        self.keys = [s['key'] for s in m.specs]
        centre = np.array([m.knobs[k] for k in self.keys])
        width = np.array([BOX_ERRORS * s['start_error'] for s in m.specs])
        below, above = (np.array([box[k][i] if box and k in box else w for k, w in zip(self.keys, width)]) for i in (0, 1))
        self.lo = np.maximum(centre - below, [s['min'] for s in m.specs])
        self.hi = np.minimum(centre + above, [s['max'] for s in m.specs])
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
            self.best_terms = {'measured': m.measured_terms(), 'true': m.spec_terms()} if m.has_beam() else None
        if self.on_pass:
            self.on_pass(self)
        if self.stop():
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
        return list(m.measured_terms().values()) if m.has_beam() else [NO_BEAM * (2 - m.reach())] * len(m.spec_terms())

    def log_terms(self, u):
        """The same terms for a model: log, clamped. Without beam worse than any beam, closer on the radial probe."""
        self.set(u)
        m = self.m
        if not m.has_beam():
            return [np.log(TERM_CLAMP) + 1 - m.reach()] * len(m.spec_terms())
        return np.log(np.clip(list(m.measured_terms().values()), 1 / TERM_CLAMP, TERM_CLAMP)).tolist()


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


def cma_es(run, rng):
    """Covariance matrix adaptation on the score, default population. Restarts from a random point when it stops."""
    import cma
    start = run.start
    while True:
        es = cma.CMAEvolutionStrategy(start, CMA_SIGMA, {'bounds': [0, 1], 'seed': int(rng.integers(1, 2**31)), 'verbose': -9, 'verb_log': 0, 'verb_disp': 0})
        while not es.stop():
            solutions = es.ask()
            es.tell(solutions, [run(np.asarray(x)) for x in solutions])
        start = rng.uniform(0, 1, len(start))


def _fitted_gp(x, y):
    """One independent GP per column of y: Matern 5/2 with automatic relevance and inferred noise, BoTorch's dimension-scaled prior."""
    import torch
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.outcome import Standardize
    from botorch.models.utils.gpytorch_modules import get_covar_module_with_dim_scaled_prior
    from gpytorch.mlls import ExactMarginalLogLikelihood
    outputs = y.shape[-1]
    kernel = get_covar_module_with_dim_scaled_prior(x.shape[-1], torch.Size([outputs]) if outputs > 1 else torch.Size(), use_rbf_kernel=False)
    gp = SingleTaskGP(x, y, covar_module=kernel, outcome_transform=Standardize(outputs))
    fit_gpytorch_mll(ExactMarginalLogLikelihood(gp.likelihood, gp))
    return gp


def _per_term_bo(run, rng, acquisition):
    """20 points as Bayesian optimisation, then one model per spec term. acquisition(gp, x): an acquisition function on the log terms."""
    import torch
    from botorch.exceptions.errors import ModelFittingError
    from botorch.optim import optimize_acqf
    from linear_operator.utils.errors import NotPSDError
    torch.manual_seed(int(rng.integers(2**31)))
    n = len(run.start)
    xs = [run.start, *rng.uniform(0, 1, (BO_INITIAL - 1, n))]
    ys = [run.log_terms(x) for x in xs]
    bounds = torch.stack([torch.zeros(n, dtype=torch.double), torch.ones(n, dtype=torch.double)])
    while True:
        x, y = torch.tensor(np.array(xs), dtype=torch.double), torch.tensor(np.array(ys), dtype=torch.double)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                candidate, _ = optimize_acqf(acquisition(_fitted_gp(x, y), x), bounds, q=1, num_restarts=ACQ_RESTARTS, raw_samples=ACQ_RAW_SAMPLES)
            next_x = candidate[0].numpy()
        except (ModelFittingError, NotPSDError):
            run.model_failures += 1
            next_x = rng.uniform(0, 1, n)
        ys.append(run.log_terms(next_x))
        xs.append(next_x)


def constrained_bo(run, rng):
    """One GP per spec term, the point most likely to meet them all."""
    from botorch.acquisition.analytic import LogProbabilityOfFeasibility
    _per_term_bo(run, rng, lambda gp, x: LogProbabilityOfFeasibility(gp, {i: (None, 0.0) for i in range(gp.num_outputs)}))


def mobo(run, rng):
    """qLogNEHVI on the spec terms, reference point at the spec limits: objectives -log(term), in spec above 0."""
    import torch
    from botorch.acquisition.multi_objective.logei import qLogNoisyExpectedHypervolumeImprovement
    from botorch.acquisition.multi_objective.objective import WeightedMCMultiOutputObjective
    from botorch.acquisition.multi_objective.utils import get_default_partitioning_alpha

    def acquisition(gp, x):
        n = gp.num_outputs
        return qLogNoisyExpectedHypervolumeImprovement(gp, torch.zeros(n, dtype=torch.double), x, prune_baseline=True, alpha=get_default_partitioning_alpha(n),
                                                       objective=WeightedMCMultiOutputObjective(-torch.ones(n, dtype=torch.double)))

    _per_term_bo(run, rng, acquisition)


def turbo(run, rng):
    """TuRBO-1 on the score: Thompson sampling in a trust region around the best point, restarted on collapse."""
    import torch
    from botorch.exceptions.errors import ModelFittingError
    from botorch.generation.sampling import MaxPosteriorSampling
    from linear_operator.utils.errors import NotPSDError
    torch.manual_seed(int(rng.integers(2**31)))
    n = len(run.start)
    candidates, perturb = min(5000, 200 * n), min(1.0, 20 / n)
    first = run.start
    while True:     # one trust region
        length, success, failure, fail_limit = 0.8, 0, 0, max(4, n)
        xs = [first, *rng.uniform(0, 1, (2 * n - 1, n))]
        ys = [-run(x) for x in xs]
        while length >= 0.5**7:
            x, y = torch.tensor(np.array(xs), dtype=torch.double), torch.tensor(np.array(ys), dtype=torch.double).unsqueeze(-1)
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    gp = _fitted_gp(x, y)
                    scales = gp.covar_module.lengthscale.detach().squeeze()
                    scales = scales / scales.mean() / scales.prod().pow(1.0 / n)
                    centre = x[y.argmax()]
                    lo, hi = (centre - scales * length / 2).clamp(0, 1), (centre + scales * length / 2).clamp(0, 1)
                    pert = lo + (hi - lo) * torch.tensor(rng.uniform(0, 1, (candidates, n)), dtype=torch.double)
                    mask = torch.tensor(rng.uniform(0, 1, (candidates, n)) <= perturb)
                    empty = ~mask.any(dim=1)
                    mask[empty, torch.tensor(rng.integers(0, n, int(empty.sum())), dtype=torch.long)] = True
                    pool = torch.where(mask, pert, centre.expand(candidates, n))
                    with torch.no_grad():
                        next_x = MaxPosteriorSampling(model=gp, replacement=False)(pool, num_samples=1)[0].numpy()
            except (ModelFittingError, NotPSDError):
                run.model_failures += 1
                next_x = rng.uniform(0, 1, n)
            best = max(ys)
            value = -run(next_x)
            xs.append(next_x)
            ys.append(value)
            success, failure = (success + 1, 0) if value > best + 1e-3 * abs(best) else (0, failure + 1)
            if success == 3:
                length, success = min(2 * length, 1.6), 0
            elif failure == fail_limit:
                length, failure = length / 2, 0
        first = rng.uniform(0, 1, n)


def policy(run, rng, model):
    """A trained policy, deterministic: one action per pass, every knob at once."""
    from lions.game.env import MissionEnv
    env = MissionEnv(run.m.mission_id)
    env.m = run.m
    while True:
        action, _ = model.predict(env.observe(), deterministic=True)
        env.act(action)
        run.passed()


def tune(m, method, seed, budget=BUDGET, on_pass=None, box=None, stop=None, **options):
    """Runs a method on mission m until in spec or out of budget. Same seed, same search."""
    run = Run(m, budget, on_pass, box, stop)
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
    Method('CMA-ES', cma_es, ('cma',), 'Evolution strategy on the score: learns which knobs move together, no model of the machine.'),
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
