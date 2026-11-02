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

**Test set.** The tunes in spec from `benchmarks/runs/case_t_drift`, F18: 20 conditions for PPO seed 0, and those Nelder-Mead reached. The tunes must be regenerated from the pinned code, since the Parquet files hold outcomes only.

**Metric.** In-spec agreement fraction; transmission, extraction energy and phase against radius for T0 and LIONS. Disagreement is reported as the result. There is no success threshold for this stage.

Every LIONS run writes one row to the dataset below.

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

**Cost.** About 32 s per pass of 120 particles: 221 turns to 18 MeV (`benchmarks/case_p.py`) at 145 ms per turn, an exploratory timing on a loaded machine, to be redone idle. Particles that never reach extraction need a turn cap, which adds to this. At that rate: 300 settings 3 core-h, 1000 settings 9 core-h, 3000 settings 27 core-h, on CC through laptop-driven `sbatch`. No agent runs on CC.

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

**Cost.** Unknown: settings per map times number of maps, likely above the 27 core-h of one map. To size after the one-map surrogate shows how many settings a map needs. Start only after the one-map gates, and on a user decision.
