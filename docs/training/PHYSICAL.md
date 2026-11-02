# A physical game: LIONS verdict, reduced model, surrogate (registered 2026-11-02, before running)

No CIME map is in hand, so everything here runs on a synthetic map, labelled as such in the game. Results show self-consistency of the game's cyclotron physics. They say nothing about a real machine.

## One model slot

A pass is `run(machine, ion, knobs, errors, draws) -> Beam` (`lions.game.machine.run`, tier T0). The reduced model and the surrogate return the same `Beam` from the same arguments, so they plug in at the two call sites, `Mission.refresh` and `behind.py`. LIONS is too slow for this slot: it runs in the background on a settled tune and reports a verdict.

Downstream of extraction (deflector, steerers, quadrupole, degrader, target) stays closed forms in every model.

## Machine: F18

Radius 0.45 m, 4 sectors, flutter 0.1, no spiral, 2 dees of 30 deg, protons, harmonic 2, 18 MeV (`scenarios/missions.yaml`). LIONS interpolation: Hermite.

LIONS has no space charge (`docs/guide/01_overview.md`). The T0 space charge term has no LIONS counterpart, so `source_ua` only scales current in the LIONS verdict. This is reported beside each verdict.

## Stage A: LIONS judges the tune

**Question.** How often do T0 and LIONS agree on a tune T0 calls in spec?

**Gates on the synthetic machine, before any verdict.**
- Map: at the design setting, cumulative RF phase slip below 20 deg up to extraction, with no energy gain beyond the deck's. Qr and Qz real from injection to extraction.
- Dees: energy gain per turn within 5% of T0's `energy_gain` at the same voltage and phase.

**Gate results and amendment, 2026-11-02, before any verdict.** `benchmarks/case_p.py`, single deterministic tracks, centred particle.
- Map: passed. Largest phase slip 9.44 deg, Qr and Qz real at 30 energies from 0.03 to 18.54 MeV.
- Dees: gain within 5% at 0.091 MeV and above (0.989 to 1.020), failed at 0.03 MeV (0.510). A 30 keV proton cannot lose what T0's cosine formula asks, and the crest gain matches.
- Amendment: verdict runs start at 0.091 MeV. Below it lies the invented central region, so LIONS does not judge the first turns.
- Added gate, before the verdicts: off-axis particles, since the constant flutter at small radius is unphysical near the centre. From 0.091 MeV at the crest, z0 = 3 mm (T0's `z_injection`) and 6 mm. Passes if the proton reaches 18 MeV, |z| at the theta = 0 section stays below the 15 mm of `z_aperture`, and phase slip stays below 20 deg.

**Off-axis gate result, CRIANN job 3316753 at 684ffdc.** Passed. Same setting and phase, single deterministic tracks, z at the theta = 0 section only:

| z0 | Final energy | Largest z | Largest phase slip |
|---|---|---|---|
| 0 mm | 18.040 MeV | 0 mm | 5.06 deg |
| 3 mm | 18.042 MeV | 3.35 mm | 4.97 deg |
| 6 mm | 18.050 MeV | 6.92 mm | 4.96 deg |

The rerun on CRIANN reproduces the laptop's map and dee gate numbers as printed (9.44 deg, gain ratios 0.510 and 0.989 to 1.020). Particles off the design radius or momentum remain untested.

**Test set.** The tunes in spec from `benchmarks/runs/case_t_drift`, F18: 20 conditions for PPO seed 0, and those Nelder-Mead reached. The tunes must be regenerated from the pinned code, since the Parquet files hold outcomes only.

**What a verdict counts, fixed before any run (`benchmarks/case_v.py`).** Cyclotron only, 120 macroparticles, Hermite, 360 steps per turn. The tunes are compared with the tables of `case_t_drift`, condition by condition, and any difference stops the run.
- Machine setting: `rf_mhz`, `dee_kv`, harmonic from T0's rule, field scale = `b0 (1 + field error)` / T0's nominal `b0`, the field error being its value at the in-spec pass. The scale multiplies LIONS' per-particle magnetic coefficient, since `b_factor` in a deck is only printed. A check compares it with a rewritten map.
- Start: T0's phases (`machine.start_phases`, buncher and phase error included) and heights (`machine.start_heights`, same draws), `pz` 0, 0.091 MeV. T0 starts at 0.03 MeV.
- Phase: `lions_phi0` is the one conversion, T0 phase + crest, the crest being the median fitted crest at 0.091 MeV and above in `case_p_gain.parquet`. LIONS has one RF clock per bunch, so a phase is an arrival time: each particle starts before theta = 0 on the closed orbit, under half a turn away, with its electric kick off until it arrives. The closed orbit radius is scaled by 1 / field scale, exact for a uniform field. A check compares the sign and size of the phase motion with T0's update.
- Extraction: first theta = 0 crossing at radius >= the closed-orbit radius at that section whose mean radius equals the machine radius. T0's `rho >= radius` is a mean radius, so this is its counterpart, and a stronger field extracts at higher energy in both.
- Losses in LIONS: |z| above the 15 mm of `z_aperture`, leaving the map, and the turn limit of T0 (`machine.turn_limit`), still circulating counts as lost.
- `cos(phi) < 0.2` has no LIONS counterpart: the primary verdict applies no such rule, and particles behind the RF follow the dynamics. A second verdict adds T0's rule on top, lost at the first theta = 0 crossing below it. Both are reported.
- Space charge: none in LIONS. T0's term stays in the T0 column only.
- Septum, slit, deflector, steerers, quadrupole, spec terms: T0's closed forms (`machine.extract`, `Mission.spec_terms`) on both beams. The tunes keep T0's deflector and steerer settings.

**Metric.** In-spec agreement fraction; transmission, extraction energy and phase against radius for T0 and LIONS. Disagreement is reported as the result. There is no success threshold for this stage.

Every LIONS run writes one row to the dataset below.

**Results, 2026-11-02.** `benchmarks/runs/case_v/`, `benchmarks/case_v.py` at a39faa0 on CRIANN, one deterministic pass per tune, synthetic map. Tunes regenerated on the laptop at a653200: all 40 conditions end as in `case_t_drift`. On CRIANN, Nelder-Mead reaches the same 18 conditions but 13 of them at another evaluation, and the file with the map differs by platform, so tunes come from the laptop and maps from the site.

Tunes T0 calls in spec, and what LIONS says:

| Method | Tunes | LIONS in spec | LIONS with T0's `cos(phi) < 0.2` rule |
|---|---|---|---|
| PPO, seed 0 | 20 | 9 | 18 |
| Nelder-Mead | 18 | 5 | 17 |
| All | 38 | 14 | 35 |

Median [min, max] over the 38 tunes:

| Quantity | T0 | LIONS | LIONS with the rule |
|---|---|---|---|
| Transmission | 0.324 [0.212, 0.449] | 0.393 [0.263, 0.498] | 0.340 [0.216, 0.474] |
| Current, uA | 32.9 [30.2, 44.5] | 41.5 [33.8, 55.4] | 34.4 [30.3, 47.0] |
| Energy, MeV | 17.99 [17.93, 18.12] | 18.05 [17.90, 18.29] | 18.08 [17.95, 18.33] |
| Relative spread | 0.0015 [0.0012, 0.0019] | 0.022 [0.009, 0.036] | 0.018 [0.008, 0.037] |
| Lost power, kW (limit 0.6) | 0.548 [0.427, 0.600] | 0.620 [0.448, 0.944] | 0.516 [0.413, 0.624] |

- Every LIONS failure is the `lost` term: 24 of 24, and 3 of 3 with the rule. Energy, current, power and axis stay in spec.
- 19 of 38 tunes already sit above 0.9 of the lost-power limit in T0, as a tuner stops at its first read in spec.
- Survival against start phase: T0's window is symmetric. In LIONS it runs from about -120 to +60 deg, with 88% at -120 to -60 deg and 32% at +60 to +90 deg. Particles starting at -120 to -90 deg reach about -54 deg. Not explained. Candidate, not tested: phase motion at low energy in the invented central region. A start at higher energy would separate it.
- Spread is about 15 times T0's. Not explained. Candidate, not tested: a radial oscillation picked up in the first turns, as energy is read at a fixed radius.
- Losses over the 4560 particles: 667 leave the map at a median 0.025 MeV, 631 hit the z aperture at a median 0.047 MeV, 105 reach the turn limit at a median 12.9 MeV and carry most of the lost power outside the septum.
- Checks: closed orbit closes to 1.6 nm; scale by coefficient against a rewritten map, 2e-13 MeV after 20 turns; phase motion over 60 turns +2.14 and -2.13 deg against +2.16 and -2.16 deg from T0's update. Phase at first crossing against the set phase: median 0.43 deg, maximum 5.3 deg, growing with z0, as the off-plane field is large near the centre.
- Cost: 0.97 core-h for 38 tunes, 95 s per pass (median, 59 to 135 s), at 499 to 734 turn limits on 8 cores. Stragglers set the length: survivors take a median 196 turns.

### Verdict without the central region (registered 2026-11-02, before running)

**Question.** Do the phase asymmetry of survival, the energy spread and the lost-power gap of the verdict above come from the invented central region?

**Change.** Only the start. `benchmarks/case_v.py --start-mev 1.0`, closed-orbit energy nearest above 1 MeV: 1.0399 MeV, in `benchmarks/runs/case_v_1mev`. Same 38 tunes, same map, same definitions, one deterministic pass per tune. Crest: median fit at the energies above the start, 91.40 deg against 91.33 deg.

**T0's phases at that energy.** The T0 start phases and heights are used unchanged, each particle on the closed orbit at its phase. This keeps all 120 particles and adds no T0 rule. It ignores T0's phase motion and T0's losses between 0.03 and 1.04 MeV, and any adiabatic damping of z, absent from T0. T0's phase change over that range is computed for each tune and reported with the result.

**Metrics.** Each against its value from 0.091 MeV, R = value at 1.04 MeV / value at 0.091 MeV:
- Asymmetry: survival at start phase -120 to -60 deg minus survival at +60 to +90 deg, over all tunes. 0.562 from 0.091 MeV.
- Spread gap: median over tunes of LIONS minus T0 relative spread. 0.0205 from 0.091 MeV.
- Lost-power gap: median over tunes of LIONS minus T0 lost power, primary verdict. 0.090 kW from 0.091 MeV. With T0's rule it is already -0.014 kW, so it is reported and not attributed.

**Attribution, fixed before the data.** Per effect: R at or below 0.25, it comes from the central region; R at or above 0.75, it does not; between, partly, with the central region one contributor. A start at 1 MeV removes every low-energy effect at once, so a low R does not name which one.

Also reported, with no threshold: in-spec agreement with T0 with and without the rule, survival against start phase, median spread and lost power. Lost particles carry at least their energy at loss, so a higher start can raise lost power by itself: reported with the fates of the lost particles.

## Stage B: reduced model

LIONS on the map tabulates revolution frequency, Qr and Qz against energy. The T0 turn loop reads the tables in place of its closed forms.

**Success.** On tunes not used to build the tables, in-spec disagreement with LIONS is below T0's, and a pass stays below the 0.15 s refresh on the laptop. Isochronism, phase history and tune crossings come from the map. Betatron dynamics and the central region do not.

## Surrogate of the LIONS pass (plan, thresholds fixed before any training)

**Inputs.** F18 has 8 knobs. Three are downstream (`deflector_kv`, `steer_x`, `source_ua`). Hidden errors merge with knobs: the field error scales `b0`, the phase error shifts `phase_deg`, and drift enters through the field error. That leaves 5 inputs: `rf_mhz`, effective `b0`, `dee_kv`, effective `phase_deg`, `buncher`.

**Particles are independent.** LIONS has no space charge, so one particle's fate depends on the machine setting and its own start. Phase and buncher only change the starting phase. The surrogate is therefore per particle: (`rf_mhz`, `b0`, `dee_kv`, start phase, other start coordinates that matter) maps to (survives, turns, final energy, phase against radius). `Beam` is assembled from the 120 macroparticles, as T0 does. One LIONS run gives one sample per particle at a machine setting, leaving 3 machine inputs to cover.

**Features.** From T0's update, with G the peak energy gain per turn (continuous approximation, derived here):

`sin phi(E) = sin phi0 + (2 pi / G) * integral of (f_rf / f_rev(E) - h) dE`

Candidate features: `sin phi0`, the mismatch group Delta = (2 pi / G) * that integral, the field scale, the reduced model's prediction. For a constant relative frequency error d, Delta = 2 pi h d (E - E0) / G.

**Hypothesis to test first.** Fate depends on `rf_mhz`, `b0` and `dee_kv` mainly through Delta. Rival: the absolute field changes the extraction radius enough to need more inputs. A feature group is dropped if removing it keeps validation error below the 0.5% reading noise, in physical units. Groups, not single features: Delta correlates with its parts.

**Outputs.**
- Survival: a continuous margin (peak |phi| against the limit, loss energy), then a threshold.
- Final energy, turns: regression on survivors.
- Phase against radius or energy: curves cut short when a particle dies. POD modes are counted on the stage A runs, with gappy POD or a masked coordinate-based model for the cut curves.
- Any latent space is scored after decoding, in physical units. `Beam` is built from 120 macroparticles by T0's binning, and gates are scored on it.

**Steps.**
1. Dataset: Parquet with one row per particle, machine inputs and outcomes with units; metadata holds map hash, deck, commit, interpolation. Written from stage A on, so no run is wasted.
2. Sensitivity: which start coordinates change fate, on a few settings. Removes inputs before sampling.
3. Sampling: Latin hypercube over the search box of `lions.game.tuners`, in steps of 300, 1000, 3000 settings. Then active learning near the survival edge and along optimiser paths, where a tuner spends its time and where uniform sampling is thin.
4. Candidates, same data:
   - references: T0, reduced model, nearest neighbour;
   - baselines: pass level, 5 inputs to scalars and POD coefficients, with boosted trees and a multi-fidelity Gaussian process (reduced model plus LIONS); interpolated table on the surviving features; polynomial chaos, which also gives Sobol indices for step 2;
   - learned: per-particle coordinate-based MLP on the features with a margin head, predicting LIONS minus the reduced model; grey-box ODE, the phase equation stepped in coarse turn intervals plus a small learned correction; DeepONet, settings in one branch, position along the curve in the other;
   - conditional: autoencoder for curves, only if POD and the MLP miss the noise level. Autoencoder alone, then regressor from inputs to latent, then end to end. Latent size from validation reconstruction error;
   - symbolic regression on the residual, only if few features survive.
5. Verification: the tuners run on the surrogate, and the tunes they find run on LIONS.

**Validation set.** 200 Latin hypercube settings and 200 along optimiser paths, on conditions outside 10000 to 10019. Model and hyperparameter selection use this set only.

**Held-out sets, used once, on the selected model and the references.** 200 Latin hypercube settings in the box, and 200 settings along Nelder-Mead and PPO paths on conditions 10000 to 10019.

**Gates.**
- Accuracy: on both held-out sets, transmission and energy errors (median, 95th percentile) below the 0.5% reading noise of F18, and in-spec agreement with LIONS above that of the reduced model.
- Tuning: tunes found on the surrogate pass on LIONS more often than tunes found on T0 and on the reduced model, same 20 conditions, 5 PPO seeds. An agent exploiting model errors fails here.
- Speed: a pass below 15 ms on the laptop CPU, a tenth of the refresh.
- Gradients: trees and tables give poor gradients or none. If one wins for the game, a differentiable model is still needed for sim2real, and the best differentiable candidate is reported apart.
- Stop: if no candidate beats the reduced model at 3000 settings, stop and report it. The PyTorch tracker becomes the route to speed.

**Cost.** 95 s per pass of 120 particles, median over the stage A verdict on CRIANN, turn limit included. At that rate: 300 settings 8 core-h, 1000 settings 26 core-h, 3000 settings 79 core-h, on CC through laptop-driven `sbatch`. No agent runs on CC.

**Limits.** The surrogate is tied to one map and one machine. A new map means new data: on CIME it has to be retrained, and the steps above are the recipe. LIONS bilinear data is excluded (`docs/training/SPEC.md`). The beam-size check against Hermite in `TODO.md` precedes training.

## To do, not run: map-conditioned operator

**Goal.** A digital twin of an old machine, where drift, deviations and ageing are unknown. The map becomes an input, so a new or perturbed map needs no retraining, and the map itself can be inferred from readings.

**Model.** FNO or DeepONet, with the field map and the settings as inputs and the same outputs as above. Other ageing terms (dee voltage calibration, dee alignment, RF phase) are extra inputs, with priors of their own.

**Training maps.** Two routes, to compare:
- a family of synthetic maps: sectors, flutter, spiral, isochronism error;
- augmentation of a given map with prior knowledge: harmonic bumps, shim and trim errors, field-scale drift, centring errors.

Each perturbed map must stay physical. In the median plane, derivatives and `d2bz` are recomputed from the perturbed `Bz` with the vacuum field equations, never perturbed independently. Magnitudes come from a source or are marked illustrative.

**Questions, in order.**
1. Forward: leave-one-family-out error, against retraining a surrogate on the unseen map.
2. Identifiability: hide a perturbed map, give only control-room readings (radial probe, phase, current) and ask what is recoverable. Different perturbations may give the same readings. Report which ones can be told apart and at what noise. This decides whether a twin can be calibrated from routine logs or needs dedicated measurements.
3. Inverse: amortised or Bayesian inference of the perturbation with uncertainty, checked on hidden maps.
4. Real map: a PSI map as a local test only, never in training, so no trained model carries PSI data. CIME when available.

**Cost.** Unknown: settings per map times number of maps, likely above the 79 core-h of one map. To size after the one-map surrogate shows how many settings a map needs. Start only after the one-map gates, and on a user decision.
