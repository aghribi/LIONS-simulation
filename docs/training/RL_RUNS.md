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

## Third amendment, registered 2026-11-02 before running (new analysis)

**Question.** With the score and stop rule reading what the control room reads, and a search box that holds the tunes other methods found, which methods reach the client's spec, and does the spec hold afterwards? Are the QUAL zeros of the optimisers an artefact of the box, and which spec terms still block QUAL?

Output in `benchmarks/runs/case_t3`, from `benchmarks/case_t3.py`, on CRIANN. Same 20 held-out conditions, 200 passes, 0.5 s of drift per pass. Earlier results and files stay as they are.

**Noise.**
- `Mission.tune_score`, `Run` and the end of a `MissionEnv` episode read `Mission.measured_terms`: the spec terms built from `Mission.readings`, one draw per pass. `Mission.spec_terms` and `Mission.ready` stay on the true beam and score the result.
- Quantities derive from the readings: current, energy, x, y, spot, lost power. Target power from current and energy. ISOL density from that power and the spot. QUAL flux from current and spot, sample energy from the energy reading, uniformity from spot and position. Depth of the Bragg peak moves by the range difference between the energy reading and the true energy.
- No reading has noise on spread, the radial probe or the phase probe, as before. No spec threshold changes.
- 0.5% on the energy of a 230 MeV beam moves the range by about 0.29 cm, against half widths of 0.4 and 0.8 cm (EYE, DEEP), and the sample energy of QUAL by about 5%, its tolerance (estimates from the range law). Measured in spec can then be chance. "True in spec at that point" below separates it.
- The coach still reads the true beam in part (`w_mean`, phase profile, `spot()`), and has no box. An advantage, stated here.
- PPO observes readings as before. Its reward and end of episode use the measured score and test.

**Box check, then rule, written before looking.**
- Tunes: coach and PPO with the old stop rule (true beam), the five policies of `benchmarks/runs/case_t_drift`, on the 20 held-out conditions. In-spec counts are compared with `case_t_drift` and differences reported.
- Containment, per mission and knob: is each tune inside its own condition's original box (start +- 3 `start_error`, clipped)?
- Rule: for each knob of a mission, with d = tune - that condition's start over all such tunes, the box of every condition becomes [start + min(-3 `start_error`, min d), start + max(3 `start_error`, max d)], clipped to the knob range. Other knobs are unchanged.
- The coach and PPO define the widened box, so the optimisers use their knowledge. Every optimiser is reported on both boxes. Where a mission's box does not change, its widened rows repeat the original ones, flagged.

**Methods, 200 passes per condition.** Coach, Nelder-Mead, random, Bayesian optimisation, CMA-ES, NSGA-II as in the amendments above. New, BoTorch 0.18, GP models with Matern 5/2 and automatic relevance, inferred noise, refit at every pass, double precision:
- TuRBO-1 on `tune_score`: 2d initial points, the first at the start; trust region 0.8, range 0.5^7 to 1.6, 3 successes to grow, max(4, d) failures to shrink; Thompson sampling on min(5000, 200 d) candidates, each coordinate perturbed with probability min(1, 20 / d). Restart from new random points on collapse.
- Per-term constrained BO: one GP per spec term on y = log(term), maximising the log probability that all y <= 0 (`LogProbabilityOfFeasibility`).
- MOBO: `qLogNoisyExpectedHypervolumeImprovement` on the same models, objectives -y, reference point 0 (the spec limits), `prune_baseline=True`, partitioning `alpha` from `get_default_partitioning_alpha`.

Both start like Bayesian optimisation: 20 points, the first at the start, then the acquisition optimised by `optimize_acqf` with 4 restarts and 256 raw samples. Terms are clamped to [1e-3, 1e3] before the log. Without beam every y = log(1e3) + 1 - reach, worse than any beam. A failed GP fit gives one random point for that pass and is counted. Stop on measured all terms <= 1.

PPO: 100k steps, 5 seeds, SB3 defaults, all six missions, trained on the noisy score. Registered defaults only, no tuning after seeing results. If the cost of the BoTorch methods rules out this configuration, the run stops and the change is documented as a deviation.

**Metrics, per mission and method.**
- In spec at the first measured read, out of 20, with median evaluations.
- True in spec at that point (`Mission.ready`).
- Stayed in spec: knobs held 50 more passes with drift, fraction of passes truly in spec, and whether all 50 held.
- Best spec terms, measured and true, for the blocking question.
- PPO: median and range over seeds. An agent that does not beat the coach and the classical methods is reported as such.

### Results of the third amendment

Source: `benchmarks/runs/case_t3/`, all at e1481c1 on CRIANN, Python 3.13, CPU torch at the locked version. 2480 runs, 242 core-h, no failed job, no failed GP fit. Elapsed per run, median [max]: Nelder-Mead 6 s [18], TuRBO 70 s [263], Bayesian 81 s [130], constrained BO 329 s [3253], MOBO 739 s [7864]. Peak memory of the sizing runs 1.7 GB; `sacct` returned no memory for these jobs.

**Box check.** In-spec counts of the coach and PPO match `case_t_drift` on all six missions. Per condition, PPO differs: 2 QUAL conditions change in or out of spec with the total unchanged, and 5 others reach spec at a different evaluation (EYE, DEEP, QUAL). Knobs the rule widens (tunes outside the original box):

| Mission | Tunes | Widened knobs |
|---|---|---|
| EYE | 22 | `rf_mhz` (1), `steer_x` (2) |
| DEEP | 21 | `steer_x` (4) |
| F18 | 111 | `steer_x` (1) |
| NEON | 12 | none |
| ISOL | 120 | `dee_kv` (1), `quad` (4) |
| QUAL | 9 | `quad` (9), `phase_deg` (5), `source_ua` (5), `steer_x` (5), `rf_mhz` (1) |

All 9 QUAL tunes lie outside the original `quad` box.

**Results.** Out of 20 held-out conditions: in spec at the first measured read (median evaluations) / truly in spec at that point / held all 50 passes. NEON's box does not change, so its widened rows repeat the original ones. The coach has no box.

Original box:

| Mission | Coach | Nelder-Mead | Random | Bayesian | CMA-ES | NSGA-II | TuRBO | Constrained BO | MOBO |
|---|---|---|---|---|---|---|---|---|---|
| EYE | 20 (50) / 19 / 19 | 5 (80) / 3 / 3 | 0 | 2 (104.5) / 2 / 1 | 0 | 0 | 15 (82) / 8 / 6 | 2 (85.5) / 2 / 2 | 7 (104) / 7 / 5 |
| DEEP | 13 (81) / 12 / 11 | 11 (67) / 8 / 7 | 0 | 1 (83) / 1 / 0 | 0 | 0 | 13 (84) / 13 / 10 | 1 (79) / 0 / 0 | 5 (118) / 5 / 5 |
| F18 | 11 (46) / 9 / 6 | 16 (63) / 16 / 12 | 2 (126) / 2 / 2 | 10 (108) / 10 / 7 | 6 (146.5) / 6 / 5 | 3 (130) / 3 / 3 | 20 (58) / 18 / 11 | 20 (43.5) / 20 / 17 | 20 (63) / 19 / 20 |
| NEON | 12 (110) / 12 / 12 | 0 | 0 | 0 | 0 | 0 | 1 (42) / 0 / 0 | 4 (114) / 4 / 2 | 8 (121) / 5 / 4 |
| ISOL | 20 (43) / 20 / 17 | 15 (72) / 12 / 2 | 0 | 14 (110.5) / 14 / 10 | 5 (164) / 5 / 4 | 1 (161) / 1 / 1 | 20 (65.5) / 20 / 10 | 18 (56.5) / 16 / 9 | 19 (91) / 19 / 12 |
| QUAL | 4 (62) / 4 / 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Widened box:

| Mission | Nelder-Mead | Random | Bayesian | CMA-ES | NSGA-II | TuRBO | Constrained BO | MOBO |
|---|---|---|---|---|---|---|---|---|
| EYE | 6 (73) / 4 / 3 | 0 | 1 (155) / 1 / 1 | 1 (73) / 1 / 1 | 0 | 9 (115) / 5 / 4 | 6 (150.5) / 4 / 2 | 9 (97) / 8 / 7 |
| DEEP | 12 (78) / 11 / 8 | 0 | 0 | 0 | 0 | 14 (89.5) / 13 / 10 | 3 (73) / 2 / 2 | 6 (114) / 6 / 6 |
| F18 | 16 (63) / 16 / 13 | 2 (126) / 2 / 2 | 9 (178) / 8 / 6 | 3 (111) / 3 / 2 | 3 (130) / 3 / 3 | 20 (52) / 19 / 13 | 20 (55.5) / 19 / 12 | 20 (51.5) / 19 / 17 |
| ISOL | 16 (60.5) / 15 / 4 | 0 | 12 (122) / 12 / 8 | 9 (140) / 9 / 8 | 2 (98.5) / 1 / 0 | 19 (63) / 19 / 13 | 20 (51.5) / 19 / 9 | 18 (79) / 17 / 13 |
| QUAL | 0 | 0 | 1 (150) / 0 / 0 | 0 | 0 | 13 (97) / 3 / 2 | 8 (119) / 6 / 3 | 12 (143.5) / 7 / 3 |

PPO, 100k steps, 5 seeds, per seed, median and range:

| Mission | Measured in spec | Truly in spec at that point | Held all 50 | Evaluations, median |
|---|---|---|---|---|
| F18 | 20, 20, 20, 20, 20: 20 (20 to 20) | 19, 19, 20, 20, 20: 20 (19 to 20) | 16, 14, 16, 17, 16: 16 (14 to 17) | 4 |
| ISOL | 20, 20, 20, 20, 20: 20 (20 to 20) | 20, 18, 20, 20, 20: 20 (18 to 20) | 15, 14, 13, 17, 18: 15 (13 to 18) | 4 |
| EYE | 3, 0, 0, 0, 0: 0 (0 to 3) | 2, 0, 0, 0, 0: 0 (0 to 2) | 2, 0, 0, 0, 0: 0 (0 to 2) | 19 |
| DEEP | 1, 0, 0, 0, 1: 0 (0 to 1) | 1, 0, 0, 0, 1: 0 (0 to 1) | 1, 0, 0, 0, 0: 0 (0 to 1) | 10 |
| NEON | 0 at every seed | 0 | 0 | none |
| QUAL | 0, 0, 0, 0, 8: 0 (0 to 8) | 0, 0, 0, 0, 5: 0 (0 to 5) | 0, 0, 0, 0, 1: 0 (0 to 1) | 47.5 |

PPO does not beat the coach on EYE, DEEP, NEON or QUAL.

**Reading.**
- A true-in-spec run can hold more often than it starts: the "held" count can exceed "truly in spec at that point" when drift brings a run back (F18, MOBO, original box).
- Measured in spec overstates most on QUAL, as registered: TuRBO on the widened box reads 13 and is truly in spec on 3.
- Results are single runs per condition and method. Differences of 1 or 2 out of 20 between methods or boxes are within what the box check shows for a platform change.

**QUAL.** The zeros of the eight optimisers come from the box. On the original box the true uniformity term is out at the best point of all 20 runs, for every optimiser. On the widened box three methods reach spec at some read: TuRBO 13 measured / 3 true, constrained BO 8 / 6, MOBO 12 / 7, against the coach's 4 / 4. Nelder-Mead, random, CMA-ES and NSGA-II stay at 0, Bayesian at 1 / 0. Terms still out at the best point on the widened box, true beam, out of 20:

| Method | Energy | Uniformity | Flux high | Flux low |
|---|---|---|---|---|
| TuRBO | 14 | 5 | 0 | 0 |
| Constrained BO | 10 | 9 | 1 | 1 |
| MOBO | 11 | 7 | 0 | 0 |
| Nelder-Mead | 20 | 20 | 0 | 1 |
| CMA-ES | 14 | 20 | 1 | 1 |

Energy on the sample blocks most where uniformity is solved: 3 measured against 14 true for TuRBO, 2 against 11 for MOBO. The noise on the 230 MeV energy reading hides it. Lost power is never out.

**Game.** `Mission.ready` is unchanged, so the coach's advice, the Experiment hold and delivery behave as before. The arcade draws no spec light from the noisy test, so no new flicker on screen. One change: a tuning method's IN SPEC screen now reports the measured read, so the beam a player takes after it can be out of spec on the true beam, without any warning. New methods are not in the arcade's list.
