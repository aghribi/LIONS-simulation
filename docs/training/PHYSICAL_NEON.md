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
