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
