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
