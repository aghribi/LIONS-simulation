# Mission tuning: baselines and RL (registered before running, approved 2026-11-01)

**Question.** From control-room readings, do PPO and SAC reach a client's spec on held-out conditions more often than the scripted coach,
Nelder-Mead, Bayesian optimisation and random search, at the same budget per episode?

**Task.** `lions.game.env.MissionEnv`: tune until `Mission.ready()`. Delivery is not part of it. Noise, drift and hidden errors as in
`scenarios/missions.yaml`, all illustrative.

**Metric.** Fraction of 20 held-out conditions (seeds 10000 to 10019) in spec within 200 machine evaluations. Also evaluations to spec
(median over those reached) and the final `tune_score`.

**Methods**, 200 evaluations per condition:
- Coach: rules of `lions.game.coach`, one knob per action, its fine step at 0.1 of full.
- Nelder-Mead, Bayesian optimisation (Gaussian process, Matern kernel, expected improvement on 2000 random candidates, 20 initial points),
  random search. Objective: `tune_score`, stop when in spec. Search box: the start +- 3 `start_error` per knob, clipped to the knob range.
  They set knobs directly, with no rate limit: an advantage over the coach and the policies, stated here.
- PPO, 100k training steps; SAC, 30k training steps. SB3 defaults, `MlpPolicy`, 5 training seeds. Training conditions: random seeds in 0 to
  999. Evaluation: deterministic policy, held-out conditions. Training cost is reported apart from the 200-evaluation budget.

**Scope.** Classical baselines and coach on all six missions. RL on F18 and EYE first, then the others if these are worth it.

**Reported as:** per method and mission, success fraction. For RL, median and range over training seeds. Training curves are not results.
An agent that does not beat Bayesian optimisation and the coach on held-out conditions is reported as such.

**Gaming checks, in `tests/test_game.py`, before any training:** no beam scores at most 0.3 and reads no position; every knob parked at a limit
scores at most 0.5; too much current is out of spec.

**What the simulator omits.** Everything in `lions.game.machine`'s docstring: no field map, no misalignments of the dees, no hysteresis, no
finite actuator speed beyond the knob rate, closed forms in place of tracking. A policy trained here is not a candidate for any machine.

## Results

Source: `benchmarks/runs/case_t/case_t_{classical,rl}.parquet`, from `benchmarks/case_t.py`. Classical at e6ad608, RL at fd5a5e9. Between them
`machine.py`, `physics.py` and `env.py` are unchanged and the `missions.py` diff touches attributes only.

In spec out of 20 held-out conditions, with median evaluations to spec among those reached:

| Mission | Coach | Nelder-Mead | Bayesian | Random |
|---|---|---|---|---|
| EYE | 19 (37) | 15 (84) | 1 (186) | 0 |
| DEEP | 8 (43.5) | 18 (107) | 3 (137) | 0 |
| F18 | 8 (40) | 17 (60) | 13 (154) | 3 (86) |
| NEON | 8 (97.5) | 5 (126) | 0 | 0 |
| ISOL | 20 (48.5) | 19 (84) | 10 (106.5) | 0 |
| QUAL | 3 (48) | 0 | 0 | 0 |

RL, in spec out of 20, median and range over 5 training seeds:

| Mission | PPO, 100k steps | SAC, 30k steps |
|---|---|---|
| F18 | 20 (20 to 20) | 0 (0 to 1) |
| EYE | 0 (0 to 5) | 0 (0 to 0) |

## Deviation from the registration

Drift is in the registration and in no run. `Mission.advance` applies it; `MissionEnv.step`, `Counter` and `coach_run` call `refresh` directly.
The accumulated field drift over 100 s is 1e-4 (1 sigma), against a hidden field error of 5e-4 to 1e-3. The results above are without drift.

## Post-hoc checks, not registered

Run after PPO on F18 reached spec in 3 to 4 evaluations. `benchmarks/case_t_check.py`, output `check_{classical,rl}.parquet`.

- Transient reads: rejected as an explanation. Holding the knobs 50 passes, every in-spec run stayed in spec at every pass, and a noise-free read agreed
  (170 classical, 7 PPO on EYE, 200 PPO on F18).
- Memorised sequence: rejected. Replaying one successful action sequence blind reaches spec only on its own condition, plus one more of 20 on the
  fresh set, at every training seed. The policy uses its readings.
- Reproduction: PPO on F18, 20/20 at each of 5 seeds on fresh conditions 20000 to 20019, median 5 evaluations, maximum 8.
- Knob count: PPO moves all 8 knobs per evaluation, the coach one. Evaluations to spec are not comparable as operator effort.
- Training cost: 100k machine evaluations per PPO seed, outside the 200 budget.

## Amendment, registered 2026-11-02 before running (new analysis, decided after the results above)

Output in `benchmarks/runs/case_t_drift`, code pinned at one commit in a worktree. Metric, conditions and budget as above.

- Drift is on: the field error takes one `Mission.wait(0.5 s)` per pass, for every method, in training too. `lions.game.tuners` holds the methods for the game and this benchmark.
- New baseline: NSGA-II (pymoo 0.6.2 defaults, population 40, first population the start plus random points). Objectives: the terms of `spec_terms`; without beam every term is 1e12 (2 - reach). Stop when in spec. Search box as the other optimisers. It sets knobs directly, as they do.
- PPO on all six missions: 100k steps, 5 training seeds. SAC is not rerun.
- `score` is the best `tune_score` of the run for every method. The coach had reported its last value.
- Comparison with the runs above mixes drift with a changed random stream: readings draw from the same generator as drift.
- `benchmarks/case_t_check.py` reruns the post-hoc checks, with drift in the hold. Not a registered metric.

### Second amendment, registered 2026-11-02 before running

- New baseline: CMA-ES (`cma` 4.5.0 defaults, sigma 0.3 of the search box, starting at the start point, restart from a random point when it stops). Objective `tune_score`, stop when in spec. Same box, same 20 conditions, 200 evaluations. It sets knobs directly, as the other optimisers.
- Run apart from the first amendment, into `case_t_classical_cma_es.parquet`, from a later commit. `machine.py`, `missions.py`, `env.py` and `physics.py` must be unchanged since eea287d.

## Results with drift

Source: `benchmarks/runs/case_t_drift/`, from `benchmarks/case_t.py` at eea287d, CMA-ES at 8a6e9e2. In spec out of 20 held-out conditions, with median evaluations to spec among those reached:

| Mission | Coach | Nelder-Mead | Bayesian | CMA-ES | NSGA-II | Random |
|---|---|---|---|---|---|---|
| EYE | 20 (50) | 13 (102) | 1 (184) | 1 (48) | 0 | 0 |
| DEEP | 13 (81) | 10 (103.5) | 0 | 0 | 0 | 0 |
| F18 | 11 (46) | 18 (59.5) | 11 (143) | 5 (124) | 0 | 2 (126) |
| NEON | 12 (109.5) | 3 (56) | 0 | 0 | 0 | 0 |
| ISOL | 20 (43) | 16 (65.5) | 7 (73) | 3 (82) | 1 (161) | 0 |
| QUAL | 4 (60.5) | 0 | 0 | 0 | 0 | 0 |

PPO, 100k steps, in spec out of 20 per training seed 0 to 4, and the median:

| Mission | Seeds | Median | Evaluations to spec, median over seeds |
|---|---|---|---|
| F18 | 20, 20, 20, 20, 20 | 20 | 4 |
| ISOL | 20, 20, 20, 20, 20 | 20 | 4 |
| EYE | 1, 0, 0, 0, 1 | 0 | 59.5 |
| DEEP | 8, 0, 0, 0, 0 | 0 | 12.5 |
| NEON | 0, 0, 0, 0, 0 | 0 | none |
| QUAL | 0, 0, 0, 0, 5 | 0 | 41 |

PPO does not beat the coach on EYE, DEEP, NEON or QUAL. Why it reaches all or nothing is unknown.

## Post-hoc checks with drift, not registered

`benchmarks/case_t_check.py` at eea287d, output `check.parquet`, written unlabelled. In-spec runs, knobs held 50 passes with drift, then one more pass:

| Method, mission | In spec | Stayed in spec at all 50 passes | In spec at the extra pass |
|---|---|---|---|
| PPO F18, held-out | 100/100 | 81 | 91 |
| PPO F18, fresh | 100/100 | 85 | 94 |
| PPO ISOL, held-out | 100/100 | 74 | 87 |
| PPO ISOL, fresh | 100/100 | 76 | 90 |
| Nelder-Mead F18 | 18/20 | 8 | 13 |
| Coach F18 | 11/20 | 6 | 8 |
| Nelder-Mead ISOL | 16/20 | 7 | 13 |
| Coach ISOL | 20/20 | 17 | 17 |
| Coach EYE | 20/20 | 17 | 19 |

- A first read in spec overstates the result: a run stops at its first such read, and drift can then take it out. The first runs, without drift, gave 100% on the same hold.
- Memorised sequence: replaying one sequence blind reaches spec on 5/100 for F18 and 7/100 for ISOL on held-out conditions. The policy uses its readings.
- Fresh conditions 20000 to 20019: PPO F18 100/100, median 5 evaluations, maximum 8; PPO ISOL 100/100, median 4, maximum 6.
- Reproducibility: Bayesian optimisation on F18 gives 11/20 in the registered run and 10/20 here with the same seeds. Unexplained.

## Deviation found after the runs: the score and the spec test read no noise

`Mission.tune_score` and `Mission.ready` use the true beam, not `Mission.readings`. Checked at 5% noise on F18: 30 passes give 30 different transmission readings and one score.

- The noise of the registration reaches only `readings`. Optimisers and the stop rule saw an exact objective in every run above, an advantage not stated in the registration.
- PPO observes `readings`, with noise.
- The "noise-free read" of `case_t_check.py` changes nothing, so it measures one more pass of drift. Its agreement in the first post-hoc check, and the rejection of transient reads there, say nothing about noise.
