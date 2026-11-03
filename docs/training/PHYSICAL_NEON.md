# A physical game on a CIME-like machine (registered 2026-11-03, before running)

Follows `docs/training/PHYSICAL.md`. Physical work moves from F18 to the NEON mission, on a synthetic machine built from published CIME values. No CIME field map is used. Results describe the synthetic machine, not CIME.

F18 stays as a labelled game mission. A proton machine with a septum is the wrong class for 18 MeV: that class extracts H- by stripping (IAEA TRS-465).

The NEON mission and the coach in the game do not change on this branch.

## Sources

| Quantity | Value | Source |
|---|---|---|
| K, sectors, spiral | 265 MeV, 4, none | Cyclotron compilation, entry C13, 2005 |
| Average field | 0.75 to 1.56 T | C13; Bourgarel et al., Cyclotrons 1998 |
| Largest hill field | 2.1 T | Bourgarel 1998, table 2 |
| Hill width, gaps | 44 deg, hill 0.12 m, valley 0.30 m | C13 |
| Pole diameter, extraction radius | 3.5 m, 1.5 m | C13 |
| RF | 9.6 to 14.5 MHz, harmonics 2 to 5 | C13 |
| Dees | 2 of 40 deg; 60 deg tips in the Muller central region | C13; Lieuvin et al., Cyclotrons 2001 |
| Dee voltage | 100 kV injection, 95 kV extraction, peak to ground, maximum | C13 |
| Injection | Axial, source up to 34 kV; 34 mm (Muller, h = 2, 3), 45 mm (spiral, h = 4, 5) | C13; Lieuvin 2001 |
| Extraction | 2 electrostatic deflectors 17 deg, 80 kV/cm; 2 magnetostatic channels 16 deg | C13 |
| Extraction efficiency | Typical 65%, best 85% | C13 |
| Tunes at 700 A | Scales Qr 0.96 to 1.04, Qz 0.2 to 0.45 | Bourgarel 1998, figure 8 |
| Isochronism error | About 5 G, scale 10 G | Bourgarel 1998, figure 7 |
| Deck, Mo-100 22+, h = 3 | 14.4 MHz, 51.89 kV, half-angle 20 deg, 1080 steps, deflector 1.513 m, 6.99 MV/m, 250 to 267 deg | `examples/don.120.ref_phi10_mars_2017` |

Tunes and the isochronism error are read by eye from a scanned figure. The figure's scales are used as limits, not its curves.

## Machine

NEON as in `scenarios/missions.yaml`: Ne-20 4+, 200 MeV total, 10 MeV/u, h = 3, 60 kV. It gives 13.87 MHz and 1.522 T, inside the published ranges. Ion mass: T0's 20 u, not the real mass, so T0 and LIONS share a design point.

`benchmarks/case_p.py` with `CASE_P_MISSION=NEON`. Field as on F18: isochronised B = <B>(r)(1 + eps cos 4 theta), F = eps^2 / 2 = 0.1, constant with radius.

Checked: hard-edge flutter from the published values. With hill fraction 44/90, hills at 2.1 T and average 1.56 T, valleys are at 1.04 T, and F = f(1 - f)(Bh - Bv)^2 / B^2 = 0.115. T0 uses 0.1. This assumes 2.1 T occurs at the 1.56 T level.

Labelled, not from a source:
- One harmonic only. A hard edge also has 12 theta and above.
- Constant flutter down to the centre.
- Dee gap half width 1.5 deg, as on F18. CIME's gap is not in these sources.
- No central region, no first harmonic, no trim coils.

## Gates, before any verdict

Same limits as F18, except the start energy, which comes from the dee gate:

1. **Map.** Phase slip below 20 deg up to extraction. Qr and Qz real at all design energies.
2. **Dees.** Gain per turn within 5% of T0's `energy_gain`. Verdicts start at the lowest table energy above every failure.
3. **Off axis.** From the start energy at the crest, z0 = 0, 3 mm and 6 mm: reaches 10 MeV/u with |z| at theta = 0 below 15 mm and slip below 20 deg.
4. **Published ranges, new.** At mean radii 0.2 to 1.4 m: |Qr - 1| at most 0.04, Qz within 0.2 to 0.45. Design point inside the field, rigidity and RF ranges.

Reported, no gate: peak field against 2.1 T. A cosine with F = 0.1 peaks near 2.2 T, 5% above, which is a limit of the field form.

F18 regression: the same code at `CASE_P_MISSION=F18` must reproduce 9.44 deg and the gain ratios of `PHYSICAL.md`.

A failed gate stops the work and is reported. No limit moves after results.

## Gate results, CRIANN jobs 3317824 (NEON) and 3317825 (F18), 53db11c

All four gates pass, with the start energy set by the dee gate as registered. Single deterministic tracks, centred ion unless stated.

| Gate | Limit | Result |
|---|---|---|
| Phase slip to extraction | Below 20 deg | 12.45 deg |
| Stability, 30 design energies | Half traces inside -1 to 1 | All |
| Dee gain against T0, 0.0587 to 8.30 MeV/u | 0.95 to 1.05 | 0.9996 to 1.0090 |
| Dee gain against T0, 0.02 MeV/u | 0.95 to 1.05 | 0.696, fails |
| Off axis, z0 = 0, 3, 6 mm | Reaches 10 MeV/u, z below 15 mm, slip below 20 deg | 10.02 MeV/u; z 0, 3.16, 6.33 mm; slip 6.91, 6.92, 6.75 deg |
| \|Qr - 1\|, 18 orbits, 0.2 to 1.4 m | At most 0.04 | At most 0.0225 |
| Qz, same orbits | 0.2 to 0.45 | 0.296 to 0.324 |
| <B> at 1.5 m | 0.75 to 1.56 T | 1.511 T |
| B rho from <B> | At most 2.344 T m | 2.267 T m |
| f_rf | 9.6 to 14.5 MHz | 13.87 MHz |

Verdicts start at 0.0587 MeV/u, as the registered rule gives. Off-axis runs start there.

LIONS Hermite on the written map agrees with Boris on the analytic field: at 10.3 MeV/u, Qz 0.2888 against 0.2890 and f error 9e-7 for both.

Reported, no gate: peak field 2.187 T against the published 2.1 T, 4.1% above. This is the limit of the cosine form registered above.

F18 regression passed: gain, orbit, slip and off-axis tables and `bmap.dat` are identical to the runs before this change.

Found after the data, not explained: T0's closed form gives 1.522 T where the tracked map has <B> = 1.511 T at 1.5 m, 0.7% apart. T0's knob and the map's field need a stated mapping before a verdict.

What these gates do not test: particles off the design radius or momentum, any real CIME map, the central region, the first harmonic, and extraction.

## Found in T0's NEON, for after the talk

Not changed here.
- The brief says the RF can lock on harmonic 2 or 3. At 10 MeV/u, h = 2 needs 9.25 MHz, below CIME's 9.6 MHz. The RF knob spans 9 to 15 MHz, wider than CIME's range.
- `w_injection` is 20 keV/u. A 34 kV source gives at most 6.8 keV/u for Q/A = 0.2.
- `deflector_kv` spans 3 to 30 kV as a voltage with a fitted scale. CIME's deflectors are specified at 80 kV/cm.

## After the gates

Registered separately once the gates are known: the verdict on NEON tunes (coach 12, MOBO 8, constrained 4, TuRBO 1 of 20 on T0), and the deflector geometry from the published values.

## Verdict on the NEON tunes (registered 2026-11-03, before running)

**Question.** How often do T0 and LIONS agree on a NEON tune T0 reads in spec, on this synthetic machine?

**Test set.** Every run that read in spec on the original box in case T3: coach 12, MOBO 8, constrained BO 4, TuRBO 1. The T3 table holds no field error, so tunes are regenerated with the same seeds on CRIANN. Evaluations, final knobs and the true read must match the table, or the run stops. The 4 tunes T3 found truly out of spec are judged too, and counted apart. A tune on a harmonic other than 3 is not judged, since the map and crest are built for 3; the count is reported.

**As F18, `PHYSICAL.md`.** One deterministic pass per tune, 120 macroparticles, Hermite, 360 steps per turn. Same field scale, phase conversion, T0 start phases and heights, z aperture, map edge, T0 turn limit. `cos(phi) < 0.2` applied only in the second verdict. Septum, slit, deflector, steerers, quadrupole and spec terms from T0's closed forms, so those columns are not physics. Start at 0.0587 MeV/u, from the dee gate. No space charge in LIONS.

**Field mapping.** Field scale = `b0 (1 + field error)` / T0's nominal `b0`. Scale 1 is the map as built, isochronous at T0's frequency. T0's 1.522 T and the map's mean field are not compared directly: the map holds the frequency at each energy, so the design point is the frequency.

**Extraction, two definitions, both run.**
- `mean`, primary, F18's rule: first theta = 0 crossing at or above the closed orbit with mean radius 1.5 m. r_ext = 1.5443 m.
- `energy`: the closed orbit at T0's 10 MeV/u. T0's rho = R is beta c / omega = R, an orbit length over 2 pi. r_ext = 1.5396 m. Derived here from T0's formula.

Measured on the design table: orbit length / (2 pi mean radius) = 1.0031 at all checked energies, scalloping from the flutter. So the mean radius rule gives 10.062 MeV/u at the orbit, 0.62% above target.

**Prediction, before running.** With tolerance 0.3%, the `mean` verdict fails the energy term on most tunes, from the definition alone. F18 did not show this as its tolerance is 3%. If this holds, `mean` does not measure tuning on NEON and `energy` is the informative reading. The decision between them stays with the user.

**Labelled extra, decides nothing.** Share of extracted particles above the entry radius of the Mo deck's deflector, 1.513 m at 250 deg. Another ion on another map, whose azimuth is not tied to these hills. LIONS' deflector has no blade and removes nothing.

**Metric.** In-spec agreement, per method and overall; which term fails; transmission, energy, spread, lost power against T0. No success threshold.

**Doubts stated before running.** Single pass per tune, 25 tunes, 21 truly in spec. The spread term compares a slit closed form with tracked energies. The turn limit stays T0's, not shown converged for this machine.

**Results, 2026-11-03.** CRIANN jobs 3317856 (tunes), 3317857 (`mean`), 3317858 (`energy`), all at 0a85120. `benchmarks/runs/case_v_NEON_{mean,energy}/` stay on CRIANN. All 25 tunes match case T3 in evaluations, final knobs and true read, all on harmonic 3. 12.3 core-h reserved. One deterministic pass per tune.

| Method | Tunes | T0 truly in spec | LIONS `mean` | LIONS `energy` | With T0's cos rule, `mean` / `energy` |
|---|---|---|---|---|---|
| Coach | 12 | 12 | 4 | 0 | 1 / 2 |
| MOBO | 8 | 5 | 0 | 0 | 0 / 0 |
| Constrained BO | 4 | 4 | 0 | 0 | 0 / 0 |
| TuRBO | 1 | 0 | 0 | 0 | 0 / 0 |
| All | 25 | 21 | 4 | 0 | 1 / 2 |

Failing terms over the 25 tunes, more than one per tune possible:

| Term | T0 | `mean` | `energy` |
|---|---|---|---|
| Energy | 4 | 14 | 15 |
| Spread | 0 | 16 | 17 |
| Current too low | 0 | 4 | 5 |
| Current too high | 0 | 2 | 2 |
| Axis | 0 | 4 | 0 |

Median [min, max]:

| Quantity | T0 | `mean` | `energy` |
|---|---|---|---|
| Transmission | 0.355 [0.054, 0.525] | 0.186 [0.033, 0.559] | 0.188 [0.046, 0.553] |
| Energy, MeV | 199.8 [198.0, 200.8] | 200.6 [199.2, 202.0] | 199.4 [198.0, 200.7] |
| Relative spread after the slit | 0.00076 [0.00044, 0.00098] | 0.0027 [0.0006, 0.0167] | 0.0023 [0.0006, 0.0159] |
| Lost power, kW | 0.0062 [0.0036, 0.052] | 0.0075 [0.0045, 0.047] | 0.0073 [0.0045, 0.047] |

**The registered prediction failed.** `mean` fails energy on 14 tunes and `energy` on 15, so moving the radius corrects the median offset and not the count. `energy` is no more informative, and neither definition is shown better.

**The verdict is not robust to the definition.** r_ext moves 4.7 mm and in-spec tunes go from 4 to 0, including all 4 from `mean`.

Read from the tables, not tested:
- Energy errors run from -0.38% to +1.0% against a 0.3% tolerance in `mean`. One crest turn is 0.18 to 0.53% of the energy, median 0.36%, so a turn is about the tolerance.
- Tracked spread on the 25 tunes before the slit, rms, in `mean`: median 1.3%, 0.48 to 3.7%. The spec is 0.1%.
- Coach tunes with the slit near 0.14% keep 0.03 to 0.27 in LIONS against 0.28 to 0.46 in T0, from 0.67 to 0.78 reaching extraction. The slit is T0's closed form on tracked energies.
- Low-energy losses are large. Of 3000 particles, 586 hit the z aperture and 242 leave the map, at median final energies of 0.16 and 0.08 MeV total, against 1.17 MeV at the start. The start energy is within 1.4 crest gains of zero, in the invented central region. 428 are still circulating at the turn limit, median 39 MeV.
- Survival by start phase: 0.93 at -30 to 0 deg, 0.87 at 0 to 30 deg, 0.29 at -120 to -60 deg, 0.16 at 90 to 180 deg.
- T0's cos rule moves little here (4 to 1, 0 to 2). On F18 it moved 14 to 35, so this explanation does not carry over.

Labelled extra: share of extracted particles above 1.513 m at 250 deg, median 0.114 in `mean`, 0.011 in `energy`. It decides nothing.

Checks: closed orbit closes to 0.3 um. Phase motion over 60 turns +3.19 and -3.20 deg against +3.24 and -3.24 deg from T0's update. Scale by coefficient against a rewritten map: 0 MeV after 20 turns. Phase at first crossing against the set phase: at most 0.46 deg. Test suite 200 passed at 926358e, before a one-line change to the deck's species.

Not concluded: single passes, a cosine field, no central region, closed-form slit and septum. The turn limit is not shown converged on this machine. None of this says anything about CIME.

## Diagnostic without the central region (registered 2026-11-03, before running)

**Question.** Do the tracked spread and the low-energy losses come from below 1 MeV/u?

**Change.** Start at 1.196 MeV/u, the table energy nearest above 1 MeV/u (`--start-mev 1.0`), crest from the gain rows at and above it. Same 25 tunes, knobs, field errors, draws, both extraction definitions, code otherwise as the verdict. T0 does not start there, so this is a diagnostic and not a new verdict.

**Primary metric.** Median over the 25 tunes of the rms relative energy before the slit, per definition, against 1.3% from the verdict in `mean`.

**Outcomes, set before running.**
- At most 0.2%, the spec's scale: the central region drives the spread.
- At least 0.7%, half the verdict's value: the spread comes mostly from later in acceleration.
- In between: both contribute, not separated.

Also reported: in-spec counts, failing terms, fates and their median energies, energy error, transmission. No threshold on these. Start phases keep their T0 values, so a particle far from the crest still loses energy, now from 24 MeV total.

**Results.** CRIANN jobs 3317996 (`mean`) and 3317997 (`energy`), at cfec373, in `benchmarks/runs/case_v_NEON_{mean,energy}_1mev/`. 3.4 core-h reserved. Single deterministic passes.

| | Verdict, `mean` | Verdict, `energy` | 1.196 MeV/u, `mean` | 1.196 MeV/u, `energy` |
|---|---|---|---|---|
| Rms before the slit, median | 1.315% | 1.295% | 1.119% | 1.270% |
| LIONS in spec, of 21 | 4 | 0 | 0 | 7 |
| With T0's cos rule, of 21 | 1 | 2 | 0 | 15 |
| Energy fails, of 25 | 14 | 15 | 20 | 5 |
| Spread fails, of 25 | 16 | 17 | 15 | 15 |
| Energy error, median | +0.31% | -0.30% | +0.50% | -0.19% |

**Registered outcome.** 1.12% and 1.27%, both above 0.7%: the spread of the whole bunch comes mostly from above the central region.

Found after the data, not registered:
- Spread depends on start phase. Within 30 deg of the crest, median rms before the slit is 0.33% from 0.0587 MeV/u and 0.10% from 1.196 MeV/u, on the `energy` tables. The whole bunch stays near 1.3% in both. So the central region sets the core's spread, and off-crest particles set the bunch's. 37% of extracted particles are outside 30 deg.
- The low-energy candidate of the verdict is refuted. Lost at the z aperture or map edge: 828 from 0.0587 MeV/u, 786 from 1.196 MeV/u, with median final energies of 0.13 and 0.06 MeV total after starting at 24 MeV. Their first phase is a median 115 deg from the crest. They decelerate to the centre. This follows from T0's start phases, not the start energy.
- The registered prediction on extraction holds once the central region is gone: `mean` fails energy on 20 tunes, `energy` on 5. With the central region both fail near 15.
- Of the 21, `energy` from 1.196 MeV/u passes 6 coach tunes and 1 MOBO tune. By slit: 6 of 11 below 0.5%, 1 of 10 at or above. T0's spread never exceeds 0.098%, so optimisers that left the slit open see no need to close it, and LIONS fails 9 of those 10 on spread.
- Coach tune 1 has 88 particles at extraction and zero current in both 1.196 MeV/u runs. T0's slit keeps particles near the mean of all energies, and off-crest outliers move that mean out of the core. A closed-form artefact, not physics.
- Published values, not used in any gate: CIME's longitudinal emittance is 0.15% x 10 deg rms (entry C13), and a semi-interceptive radial probe removes unwanted phases (Lieuvin 2001). NEON's 0.1% spread spec is tighter than that 0.15%. T0's cos rule stands in for a phase selection it does not model.

## Published phase width (registered 2026-11-03, before running)

**Change of primary definition, after seeing results.** For NEON, `energy` is the primary extraction definition from here on. The earlier tables stay as they are. It is derived from T0's rho = R, and from 1.196 MeV/u it fails energy on 5 tunes against 20. NEON and F18 verdicts then use different rules and are not directly comparable.

**Question.** With a published phase width in place of T0's start phases, does LIONS agree with T0 on the 21 tunes?

**Change.** Start phases Gaussian, 10 deg rms (entry C13: 0.15% x 10 deg rms), centred on each tune's `phase_deg` + phase error, as 120 deterministic quantiles. T0's buncher law is not used. Start at 1.196 MeV/u, `energy`, no cos rule in the primary reading. Same tunes, knobs, field errors and heights. T0's columns do not change.

**Limits stated before running.**
- C13 gives no location or definition beyond "RMS", and 10 deg is taken at the start here.
- T0's law is a full RF cycle at injection with a buncher, consistent with a CW source: Lieuvin 2001 reports a buncher catching about 65% of the CW beam. The narrow width models what survives selection, with no selection hardware and no intensity loss.
- Current is therefore inflated, as transmission counts 120 particles that were not cut. Current terms are not meaningful in this run, and the count without them is the primary one.
- Tunes were found by T0 with its own phases. A good result shows where the mismatch comes from, not that the tunes are right.

**Primary metric.** LIONS in spec without the current terms, of T0's 21.

**Outcomes, set before running.**
- At least 17 of 21: the phase width closes the gap.
- At most 11 of 21: it does not, and another cause is needed.
- Between: partial.

**Also reported.** Rms before the slit, against 0.15% published and 1.27% from T0's phases; each failing term; in spec with current terms; energy error.

**Results.** CRIANN job 3318020 at bc503af, `benchmarks/runs/case_v_NEON_energy_1mev_phase10/`. 1.1 core-h reserved. Single deterministic passes. This run's Parquet files predate the metadata fields added at bd386d2: its settings are in the directory name, this section and the job script `~/lions/neon_phase.sbatch`.

**Registered outcome: partial.** LIONS is in spec without the current terms on 12 of T0's 21, between 11 and 17.

| | T0's phases, 1.196 MeV/u, `energy` | 10 deg rms, same |
|---|---|---|
| In spec without current, of 21 | 9 | 12 |
| Coach, of 12 | 8 | 9 |
| MOBO, of 5 | 1 | 2 |
| Constrained BO, of 4 | 0 | 1 |
| In spec with current, of 21 | 7 | 9 |
| Rms before the slit, median (range) | 1.270% (0.256 to 3.558), 25 tunes | 0.125% (0.058 to 4.088), 22 tunes |
| Spread fails, of 25 | 15 | 10 |
| Energy fails, of 25 | 5 | 8 |
| Axis fails, of 25 | 1 | 6 |
| Tunes with no particle at extraction | 0 | 3 |

The median rms, 0.125%, is below the published 0.15%, and against 1.27% from T0's phases. Current terms stay uninterpreted, as registered: 8 too high and 5 too low.

Found after the data, not registered:
- The 9 of T0's 21 that fail. Five on spread alone: coach tunes 8, 9, 10 at 1.05 to 1.10 of the limit, MOBO tune 20 at 1.16, constrained BO tune 13 at 7.9. Three on energy, spread and axis: constrained BO 14 and 15, MOBO 22. One with no particle at extraction: MOBO 18.
- Counting the four marginal spread failures as passes would give 16, still under 17.
- Six tunes with phases of 40 to 90 deg keep at most 0.21 of the beam to extraction: 14, 15, 18, 22, 23, 24. Phase alone does not decide it: tunes 12, 16, 20, 21, also at 41 to 69 deg, keep 0.88 to 1.0, not explained. T0's full RF cycle gives such a tune particles near the crest anyway, so the phase knob costs little there. A narrow bunch has none. Tune 13 is the exception: phase 12 deg, 0.07 reach extraction, not explained.
- Near the crest the agreement is better: the 12 coach tunes all reach extraction with all 120 particles and rms 0.105 to 0.163%.
- Tunes were found with T0's phases, so they are not tuned for a narrow bunch. A method retuning under the narrow law may differ.

## Reduced model on NEON (registered 2026-11-03, before running)

**Question.** Does the table model of `PHYSICAL.md` stage B describe LIONS on the NEON machine, with the published phase width?

**Model.** `lions.game.reduced`, same rules and same checks as stage B: energy from LIONS crest gains, phase from f_rev(p, s) = v(p) / L(p / s), Betatron action, no phase cut, stall, clamp or septum in the turn loop. Radial motion and the central region are absent, as there.

**NEON settings.** Start at 1.196 MeV/u, orbit-length extraction, Gaussian start phases of 10 deg rms from the same quantiles as LIONS. Table `scenarios/neon_synthetic_orbits.parquet`: LIONS Hermite closed orbits at the 30 design energies, no RF, independent of the tunes. Not connected to the arcade.

**Reference.** The 25 tunes of `case_v_NEON_energy_1mev_phase10`, all held out from the tables. T0 unmodified, with its own phases, is the baseline: it is what the game runs.

**Success.** Median absolute error over the tunes, each line tied to a NEON quantity. Tunes where LIONS extracts nothing have no energy and are left out of the energy line, counted.

| Quantity | Threshold | Tied to |
|---|---|---|
| Transmission before the septum | At most 1/120 | One macroparticle, as F18 |
| Mean extraction energy, relative | Below 0.3% | NEON's energy tolerance |
| Relative energy spread | Below 0.1% | NEON's spread limit |
| Phase against radius, rms over bins | Below 6.91 deg | Largest slip of a centred particle from the off-axis gate, as F18 |
| Pass time | Below 0.15 s | Refresh, as F18 |

The model must pass every line and be closer than T0 on each. Thresholds do not change after the data. Reported besides: rule checks, per-particle fates, survival against start phase, turn limit at 1, 2, 4 times T0's, and the same lines on the coach tunes, which sit near the crest.

Current stays labelled and is not scored: no source gives the intensity cut for a 10 deg bunch here. Passing here would allow a retuning test with this model. It would not make the model true of CIME.

**Result: fails by the registered rule, on one line.** Four lines pass. The energy line meets its threshold but T0 is slightly closer. `benchmarks/case_r.py`, CRIANN jobs 3318071 (tables, checks) at 6a61a0f and 3318080 (comparison) at d6471dd. One deterministic pass per tune, against `case_v_NEON_energy_1mev_phase10`. Synthetic map.

| Quantity, median absolute error [min, max] | Threshold | Reduced | T0 | Outcome |
|---|---|---|---|---|
| Transmission before the septum, 25 tunes | 0.0083 | 0 [0, 0.208] | 0.342 [0.05, 0.70] | Pass |
| Mean extraction energy, relative, 22 tunes | 0.3% | 0.058% [0.001%, no beam] | 0.052% [0.001%, 4.6%] | Below threshold, T0 closer: **fail** |
| Relative energy spread, 22 tunes | 0.1% | 0.016% [0.001%, no beam] | 0.022% [0.001%, 4.0%] | Pass |
| Phase against radius, rms, deg, 22 tunes | 6.91 | 0.20 [0.09, 15.7] | 7.04 [2.3, 67.6] | Pass |
| Pass time | 0.15 s | 27 ms median, 112 ms max, one core, `debug` | 41 ms | Pass |

Coach tunes alone, 12 of 12 with beam, give the same outcomes: transmission 0 against 0.346, energy 0.058% against 0.054%, spread 0.008% against 0.013%, phase 0.15 against 6.56 deg.

Three tunes with no LIONS beam are left out of energy, spread and phase, as registered. Three tunes where LIONS extracts 13, 20, 25 of 120 and the model none count as infinite error.

Corrections made after seeing the first printed tables, each committed. They changed the outcomes: at c1e97a6 transmission and spread failed and energy passed; at d6471dd transmission and spread pass and energy fails. The fixes follow the line names and the thresholds did not move, but they are not blind:
- The first comparison dropped missing values per model and wrote a spread of 0 for a tune with no beam. Replaced by one mask, and no beam in a model counts as failed (c1e97a6).
- `Beam.energies` holds what the slit keeps, and the LIONS columns are before the slit. F18 never closes its slit, so it did not show there. All three models are now read with the slit open, as the line is named (d6471dd).

Rules against their checks:

| Rule | Result |
|---|---|
| f_rev(p, s) = v(p) / L(p / s) | Verified: within 2.0e-7 of Boris closed orbits at s = 0.997 and 1.002, 0.56 and 5.6 MeV/u |
| Gain linear in V | Verified: 40 to 60 kV ratio 1.4999 at 4.35 MeV/u, 1.5008 at 8.30 MeV/u |
| Path-length interpolation | At most 4.7e-6 with every other energy held out |
| Betatron action | 3.126 mm against 3.163 mm tracked at z0 = 3 mm, 6.253 against 6.328 at z0 = 6 mm, top band. Within 1.2% |
| Numerical turn limit | Converged: extracted 2148, 2151, 2151 at 1x, 2x, 4x; no tune changes from 2x to 4x. LIONS' own limit was not studied here |

Per particle, of 3000: extracted in both 2142, LIONS only 65, reduced only 6, neither 787.

| Start phase, deg | -120 to -90 | -90 to -60 | -60 to 0 | 0 to 60 | 60 to 90 | 90 to 120 |
|---|---|---|---|---|---|---|
| LIONS | 0.00 | 0.23 | 0.89 | 0.79 | 0.46 | 0.16 |
| Reduced | 0.00 | 0.22 | 0.89 | 0.77 | 0.42 | 0.03 |

Found after the data:
- The energy miss is 0.006 percentage points on single passes, a fiftieth of the tolerance. The registered rule still fails it. Only a user decision can set that aside.
- The largest errors sit on tunes far from the crest with a small surviving beam: 13, 14, 15, 22. LIONS' spread there is 1.6 to 4.1%, against under 0.1% in the model. No radial motion in the model is a candidate, not tested.
- Lost-particle classes differ: LIONS 420 at the z aperture, 135 at the map edge, 238 at the turn limit. Reduced 0, 836, 16. They agree on which particles are lost far more than on how.
- T0's baseline keeps its own wide phases. Its transmission error measures that difference more than its turn loop.
- Current is still not scored. The comparison reads the slit open, so slit behaviour on tracked energies is not tested here.

Test suite: 200 passed at eac9a62, CRIANN job 3318081.

Nothing here makes the model true of CIME. Not connected to the arcade.
