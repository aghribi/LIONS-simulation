# LIONS training game: specification (draft for review)

Status: draft, nothing built. Decisions marked **D1** to **D8** are open for A. Ghribi and O. Kamalou.
"Lecture §n" is chapter n of the lecture notes in `docs/lecture`. Other "§n" are sections of this document.

## 1. Purpose

A game-style trainer for sector cyclotrons, built on the LIONS physics, with a GUI. One tool, four uses:

| Audience | Wants | Track |
|---|---|---|
| Students | how a cyclotron works, from the Lorentz force up | Discovery (levels 0 to 2) |
| Engineers | what each knob does to the beam, field maps, tolerances | Machine (levels 1 to 4) |
| Operators | tune CIME faster and more safely, recover from drifts | Control room (levels 3 to 5) |
| Physicists, ML researchers | compare tuning strategies, train and test RL policies | Lab (all levels, scripted) |

The same scenario runs for a human at the GUI, for a classical optimiser, and for an RL policy, on the same budget and the
same observations. That makes the game both a trainer and a benchmark.

## 2. Principles

- **Physics from LIONS.** Every level runs the LIONS physics or a model validated against it. A level whose fast model
  disagrees with LIONS beyond its stated tolerance does not ship.
- **Two models, one truth.** A fast model for play and RL training (PyTorch port, Lecture §10); LIONS for validation and
  for the scored end of each run. Disagreement is shown, not hidden.
- **Control-room observations.** In the operator track the player sees what an operator sees: radial probe currents,
  phase probes, extracted current, beam on viewers. The full phase space is a student-track view, labelled as such.
- **Sim-to-real is a setting.** Noise, drift, hysteresis, finite actuator speed and misalignments are switches with
  stated magnitudes (§6). Every score records which were on.
- **Scores are results.** Budget, seed count and baseline are recorded with every score (§7).

## 3. Engine

| Tier | Model | Use | Speed (measured) |
|---|---|---|---|
| T0 | closed forms (uniform field, weak focusing, smooth approximation; Lecture §2, §5) | level 0, instant feedback | analytic |
| T1 | LIONS, single particles | levels 1 and 2, validation | about 1e-7 s per particle-step (Lecture §9) |
| T2 | PyTorch port, batched, CPU or GPU | levels 3 to 5, optimisers, RL | 1.4e9 particle-steps/s on H100, 9e7 on M3 CPU at 1e5 particles (Lecture §10) |
| T3 | LIONS, full bunch | scored end of a run, held-out evaluation | as T1 |

Gaps to close before levels 3 to 5 (not built):
- The PyTorch port has the bilinear field and the Boris push only: no RF dees, no central region, no extraction elements.
  It must reach parity with LIONS on these, validated as in Lecture §7 (field parity, then tracking parity).
- The CIME field map is needed for the CIME scenarios (**D1**: availability and licence of `cime700rec.1804.pbob1.lio.dat`
  and the CIME decks for training use).
- Field model: `interpolation=hermite` for tune and Jacobian work (Lecture §9); the bilinear default for Fortran
  reproduction.

## 4. Levels

Each level states the learning goal, the controls, what is shown, the success test and the physics source.

### Level 0: how a cyclotron works (Discovery)
Short interactive scenes, T0 and T1, no score pressure:
1. **Lorentz force.** Drop a proton in a uniform B: circle, radius rho = p/(qB), frequency qB/(gamma m). Change B, energy,
   charge, mass; watch radius and period. Sign of B and direction of rotation.
2. **The dees.** Two D-shaped electrodes with an RF voltage. Find the RF frequency that accelerates; see the spiral.
   Wrong frequency: the particle falls out of phase.
3. **Relativity breaks the clock.** Accelerate to high energy in a uniform field: gamma grows, the revolution slows,
   the beam slips in phase. Fix: make B grow with radius (isochronism, B = gamma B0).
4. **Then the beam falls down.** A field growing with radius defocuses vertically (Qz^2 = -k). Watch particles leave the
   median plane.
5. **Hills and valleys.** Add sectors (flutter) and spiral: vertical focusing returns (Thomas). Read Qz from the
   smooth-approximation formula and from tracking side by side.
6. **The real machine.** Central region, extraction by a deflector, the CIME layout.

Success: a short check per scene (predict the radius before running, find the accelerating frequency, choose the
flutter that restores Qz above a threshold).

### Level 1: single-particle dynamics (Discovery, Machine)
Closed orbit, betatron oscillations, tunes. Controls: field level, gradient, flutter, spiral angle, starting offset.
Shown: Poincaré section at theta = 0, turn-by-turn signal, its spectrum. Success: read Qr and Qz within a tolerance from
the section or the spectrum. Physics: Lecture §2, §6.

### Level 2: find the closed orbit (Discovery, Machine)
The player adjusts r and pr at the section until the orbit closes; then the Newton method does it in three steps.
Shows why a near-integer Qr makes this hard (sensitivity 1/(2 sin pi Q), Lecture §6). Success: residual below 1e-4 m in
a given number of turns tracked.

### Level 3: accelerate a bunch (Machine, Control room)
Controls: B level (b_factor), dB/B, dee voltages and phases, RF frequency. Observables: radial probe current versus
radius, phase probe, beam losses. Success: transmission to a target radius above a threshold, energy spread below a
threshold, within a knob-change budget.

### Level 4: extract (Machine, Control room)
Add the electrostatic deflector (position, voltage) and the extraction grids. Success: extracted current and energy on
target, losses on the septum below a threshold.

### Level 5: tune CIME (Control room, Lab)
The full CIME deck from a perturbed start, sim-to-real switches on. Success: the operator target (transmission,
extracted energy, beam quality) reached within the time and knob budget. Variants: recover from a drift; change ion
species; restart after a field cycle (hysteresis).

## 5. GUI

**D2**: build on the existing PyQt5 GUI (`lions.gui`: parameter form, live plots, field viewer; `uv run lions-gui`) or a
separate web front end. Recommendation: extend the PyQt5 GUI for the first version. The engine, form and plots already
exist, and the scenario API (§8) keeps a later web front end possible.

Screens:
- **Level map** with tracks and progress.
- **Machine view**: the cyclotron seen from above, field map, orbits, dees, deflector; animated particles.
- **Control panel**: knobs with real ranges and step sizes, actuator speed shown when on.
- **Diagnostics**: probe currents, phase, extracted current, viewer images. A "god view" (phase space, tunes, closed
  orbit) available in the Discovery and Machine tracks only.
- **Coach**: hints tied to the lecture sections; after a run, a replay with the knob trace and what changed.
- **Arena**: the same scenario played by a human, an optimiser and a policy, side by side or as a replay.

## 6. Sim-to-real switches

| Effect | Model | Default magnitude |
|---|---|---|
| Diagnostic noise | Gaussian on probe readings | **D3** (needs CIME probe data) |
| Diagnostic saturation | clip at a full-scale value | from probe specifications (**D3**) |
| Field drift | slow random walk on b_factor | **D3** |
| Hysteresis | field depends on knob history | **D3** |
| Actuator speed | knob setpoint reached at a finite rate | **D3** |
| Misalignments | dee and deflector offsets drawn per episode | **D3** |
| Initial beam | sigma matrix drawn around nominal | from the CIME deck |

Each magnitude needs a source (machine data or operator knowledge). Until it has one, the switch is labelled
"illustrative" in the GUI and in every score.

## 7. Measures and baselines

Every scored run records: scenario, seed, switches, budget (knob changes, LIONS evaluations, wall time), the final
metrics, and the code version.

Metrics: success (target met), time and number of evaluations to target, transmission, extracted energy and spread, losses,
and the number of unsafe actions (§8).

Baselines, at equal evaluation budget and on the same observations:
- **Operator procedure**, scripted from the CIME tuning practice (**D4**: who writes it down; Omar and the operators).
- **Nelder-Mead**.
- **Bayesian optimisation**.
- **Random search**.

Results are reported over at least 5 seeds with the spread, never as a single run.

## 8. RL and policies

### Environment
A Gymnasium environment per level (3 to 5), on the same scenario API as the GUI:
- **Observation**: the control-room diagnostics of §5, with the active noise. Never the hidden state (closed orbit, true
  phase space) in the operator and lab tracks.
- **Action**: bounded increments of the level's knobs, at most one step size per action, as an operator would turn them.
  Safe ranges enforced by the environment.
- **Episode**: a perturbed start drawn from the scenario distribution; ends at target, at budget, or on an unsafe state.
- **Reward**: progress toward the target metrics, minus a cost per action. **Gaming checks**, written before training:
  what the reward gives when the beam is lost (must be the worst value, not a neutral one), when a probe saturates (must
  not read as a high current), and when the agent parks at the edge of a safe range.

### Training and evaluation
- Train on the fast model (T2) with domain randomisation over the switches of §6.
- Evaluate on **held-out conditions**: seeds, misalignments, drift profiles and map variants never seen in training,
  plus the scored end on LIONS (T3). A policy tested only on its training distribution is not reported as generalising.
- At least 5 training seeds, with the spread. Training curves are diagnostics; the held-out evaluation is the result.
- Against the baselines of §7 at equal budget. An agent that does not beat Bayesian optimisation and the operator
  procedure on held-out conditions is reported as such.
- Libraries: **D5** (Stable-Baselines3 for PPO and SAC first; PufferLib if throughput limits).

### Policies in the game
- **Watch**: replay a policy on a scenario, with its observations and actions.
- **Duel**: human against policy on the same seed and budget; leaderboard for humans, optimisers and policies.
- **Advisor**: the policy suggests the next knob change; the player accepts or not. Logged.
- **Inspect**: action traces, sensitivity of the policy to each observation.

### Human data
Every human run is logged (observations, actions, outcomes) as Parquet with units and metadata. It feeds imitation
learning and offline RL, and a comparison of human and agent strategies. **D6**: consent and anonymisation for operator
and student logs.

### Toward the real machine
Out of scope for this specification. A policy is a candidate for CIME only after held-out evaluation, a named sim-to-real
gap analysis, and read-only shadow mode on machine data. Any setpoint write needs explicit approval, bounded steps and
logged previous values.

## 9. Scenario API

One scenario definition drives the GUI, the baselines and the RL environment:
- a LIONS deck (or the level's T0 model), the knobs with ranges and step sizes, the observables, the target, the budget,
  the switches and their magnitudes, and the seed;
- `reset(seed) -> observation`, `step(action) -> observation, metrics, done`, `score() -> record`;
- format: **D7** (YAML on top of the existing `Deck` dataclass is the obvious option).

## 10. Validation before release

- Level 0 and 1 physics: closed forms against LIONS (case A and B0 tolerances, Lecture §8).
- Fast model against LIONS on every level scenario: closed orbit, tunes, energy gain, transmission, within stated
  tolerances, registered before the comparison.
- Baselines reproduce across seeds; an optimiser that never reaches the target on a scenario flags the scenario.
- Gaming checks of §8 run as tests.

## 11. Versions

| Version | Content |
|---|---|
| v0 | Level 0 scenes (T0, T1) in the PyQt5 GUI; scenario API; logging |
| v1 | Levels 1 and 2; Arena with Nelder-Mead, BO, random search |
| v2 | PyTorch port with dees and extraction, validated; levels 3 and 4; Gymnasium environments; first policies |
| v3 | CIME (level 5) with the CIME map and decks; sim-to-real switches with sourced magnitudes; held-out evaluation |

**D8**: priority between the educational track (v0, v1) and the RL track (v2) if they compete for time.

## 12. Open decisions

| Id | Decision |
|---|---|
| D1 | CIME map and decks: availability and licence for a training tool |
| D2 | GUI platform: extend PyQt5 (recommended) or web |
| D3 | Magnitudes of noise, saturation, drift, hysteresis, actuator speed, misalignments (CIME data) |
| D4 | Who writes down the operator tuning procedure |
| D5 | RL library |
| D6 | Consent and anonymisation of human logs |
| D7 | Scenario file format |
| D8 | Priority between educational and RL tracks |
