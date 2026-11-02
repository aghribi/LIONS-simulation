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
