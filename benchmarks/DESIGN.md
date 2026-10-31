# LIONS Phase 5: benchmark design note (stage 1, symplectic tracking)

Status: case A run (section 10). Other cases not run.
Scope of stage 1: single-particle and small-bunch tracking in external fields only.
Out of scope until stage 1 is closed: space charge, 3D map reconstruction, anomalies, calibration on CIME data.

## 1. Questions

| Id | Question | Decides |
|---|---|---|
| Q1 | Do LIONS, OPAL-cycl and RF-Track agree on the same fields, at the same time step? | trust in the LIONS integrator and field model |
| Q2 | Is the LIONS one-turn map symplectic, and how does it compare with OPAL LF2/RK4 and RF-Track? | whether "symplectic first" needs a new integrator |
| Q3 | Can the LIONS core (fields + pusher) be made differentiable at an acceptable memory cost? | whether a PyTorch backend is worth building |
| Q4 | Do the codes agree on acceleration through dee gaps (energy gain, phase slip)? | gap model parity before CIME maps |
| Q5 | Does the beam handed over at extraction survive the conversion to Xsuite and Cheetah? | the integration route (stage 1 end) |

## 2. Codes

| Code | Version | Source | Integrators used |
|---|---|---|---|
| LIONS | main at the branch point | this repo | Boris leap-frog (kick, then drift) |
| OPAL-cycl | 2024.1 (tag `OPAL-2024.1.0`) | github.com/OPALX-project/OPAL, container `lions-bench/opal:2024.1.0` | `LF2` (drift half step, Boris kick, drift half step), `RK4` |
| RF-Track | 2.7.4 | PyPI `RF_Track` | `leapfrog`, `rk2` (default), `rkf45` |

Read from the installed sources, not from memory:
- OPAL `Cyclotron::interpolate` and `Cyclotron::apply`: bilinear interpolation in (r, theta);
  Br = z dBz/dr, Btheta = (z/r) dBz/dtheta, Bz = Bz(r, theta). No z^2 term.
  dBz/dr and dBz/dtheta are computed by OPAL from the Bz map (5-point differences, `getdiffs`).
  The map start angle must be 0 ("this does not work if the start point of field map does not equal zero").
- OPAL `CARBONCYCL` map: header rmin [mm], dr [mm], theta_min [deg], dtheta [deg], ntet, nr, then Bz[nrad][ntet].
  Units kG (the FFA and CYCIAE readers convert T to kG; Carbon reads values as given).
- OPAL RF (`TYPE=BANDRF`): list of 3D Cartesian maps (`RFMAPFN`), each scaled by `ESCALE * cos(2 pi RFFREQ t + RFPHI)`, superposed when `SUPERPOSE` is true.
- RF-Track: `Volume` + `Bunch6dT` for time tracking (mm, MeV/c, time in mm/c).
  Fields through `UserField.get_field(x, y, z, t)` (mm, mm/c -> V/m, T), or `Static_Magnetic_FieldMap` and `RF_FieldMap` on 3D meshes.

## 3. Field parity

Each code must see the same field. Rule: every converted field is sampled from LIONS's own field evaluation, so the conversion cannot add a second field model.

| Item | LIONS | OPAL | RF-Track | Parity action |
|---|---|---|---|---|
| Bz map | (r, theta) grid, 0..360 deg incl. endpoint, T | CARBONCYCL, kG, no endpoint, start 0 | 3D mesh or callback | drop the 360 deg column, T -> kG, check the sign convention on case A |
| Off-plane B | z dBz/dx, z dBz/dy from the file, - z^2/2 d2Bz/dr2 | first order, derivatives computed from Bz | as given | parity cases use d2Bz = 0; the z^2 term is measured separately in LIONS alone |
| dBz/dr | read from file | 5-point difference of Bz | as given | analytic maps, so the OPAL difference error is O(h^4); stated in the report |
| Gap E | potential maps, `V cos(w t + phi)` per dee | 3D E map per dee, `ESCALE cos(...)` | `UserField` or `RF_FieldMap` | sample LIONS E on a Cartesian grid, one map per dee |
| RF time convention | cos(w t + phase) | cos(2 pi f t + RFPHI) | `RF_FieldMap` convention unverified | check on a one-gap test before case C |

## 4. Conversion functions

One module, `benchmarks/convert.py`, one function per handoff, each validated on case A before use:
- `bmap_to_opal_carbon(lions_map) -> file`
- `efield_to_opal_map(lions_state, grid) -> file` (format chosen after a format check: `FM3DDynamic` text or `FM3DH5Block`)
- `bunch_to_opal(bunch)`, `bunch_from_opal(...)`
- `bunch_to_rftrack(bunch) -> Bunch6dT`, `bunch_from_rftrack(...)`
- stage 1 end: `bunch_to_xsuite(bunch, s0)`, `bunch_to_cheetah(bunch, s0)`

## 5. Cases

All runs record: code, version, case, dt (as steps per turn), integrator, metric, value, unit.
Single particles are deterministic. Bunch cases use 5 fixed seeds, spread reported.
Tolerances are fixed here, before any run. A change after results is a new analysis, documented as such.
No quantitative CIME requirement was found (repo, iRIS project folder; the ag-rag index returned errors), so the tolerances below are numerical targets. They are replaced by CIME requirements when these exist.

### A. Uniform Bz, no E (sanity, conversions)
- Field: Bz = 1 T, proton, 1 MeV (as `tests/test_run.py`).
- Metrics: radius of the circle fitted to one turn of positions vs p/(qB); relative change of |p|; phase per step.
- Sweep: 90, 180, 360, 720, 1440 steps per turn (fixed-step integrators); tolerance 1e-6, 1e-9, 1e-12 (adaptive ones).
- Expected, derived here: Boris positions lie on a regular polygon of circumradius rho sqrt(1 + (w dt / 2)^2), so the radius error is (w dt)^2 / 8 to leading order (2.4e-6 at 1440 steps).
- Success:
  - fixed-step codes: observed order 2 (LIONS Boris, OPAL LF2) or 4 (OPAL RK4), within 0.2;
  - LIONS fitted radius equal to rho sqrt(1 + (w dt / 2)^2) within 1e-10 relative;
  - |p| constant within 1e-12 relative for the Boris-type integrators;
  - adaptive codes (RF-Track `rkf45`, `rk2`): error at or below the requested tolerance.
- Correction made before any run: the first draft asked for 1e-8 at 1440 steps, which a second-order integrator cannot reach (2.4e-6 expected).
### B. Synthetic 4-sector static field (equilibrium orbits, tunes)
- Field: Bz = B0 (1 + eps cos 4 theta), B0 = -1 T (LIONS sign for protons), eps = 0.3, d2Bz = 0;
  analytic dBz/dx, dBz/dy for LIONS. Map: r = 0.05 to 0.6 m, dr = 5 mm, 1 deg.
  eps = 0.3 gives a flutter of eps^2 / 2 = 0.045, so Qz of order 0.2 with no radial gradient.
- Energies: 1, 2, 3, 4, 5 MeV protons (rho from 0.14 to 0.32 m at 1 T).
- Section: theta = 0 (y = 0, x > 0). State at the section from 4-point Lagrange interpolation in time of
  full-step positions and momenta (LIONS momenta are at half steps: the full-step value is the mean of
  the two neighbours).
- Equilibrium orbit: Newton on (r, pr) at the section, one-turn map, finite-difference Jacobian.
- Tunes: particles offset by 1e-4 m in r or in z from the equilibrium orbit, 256 turns; the frequency of
  the Hann-windowed turn-by-turn signal, refined by maximising the windowed DFT (golden section).
  OPAL tune mode as a second instrument if its output can be read.
- Steps per turn: 360, 720, 1440; the comparison uses 1440 and the convergence across the three.
- Success: |Delta Q| < 1e-4 between codes at 1440 steps; equilibrium radius within 1e-6 relative.
### C. Acceleration through two dees (gap model)
- Field: case A or B plus two tanh gaps (left and right dee), V_dee set per dee.
- Metrics: energy gain per turn, RF phase per turn (phase slip), orbit centre drift, over 50 turns.
- Success: energy gain per turn within 1e-4 relative between codes; phase per turn within 0.01 deg.
- Prerequisite: RF time-convention check (section 3).

### D. Symplecticity (Q2)
- Map: Poincare map at theta = 0 in case B, 4D (r, pr, z, pz) at fixed energy, theta as the independent variable.
- Canonical coordinates: with the gauge A_r = (z^2 / 2r) dBz/dtheta, A_z = 0, the kinetic pr and pz equal the canonical momenta up to O(z^2). Derived in this session, first order in z; small amplitudes only.
- Jacobian J by central differences (steps 1e-6 m and 1e-6 in p), metrics: eps_S = max|J^T S J - S| and |det J - 1|.
- Crossing of theta = 0 interpolated with cubic Hermite between steps; a dt sweep separates the integrator error (scales with dt) from the interpolation floor.
- Long-term: energy drift and phase-space area over 1e4 turns in cases A and B.
- Expectation, to be tested, not assumed: Boris and LF2 keep det J = 1 to round-off; RK4 drifts in energy. Whether Boris is symplectic here is the question.

### E. Differentiability spike (Q3), LIONS only
- Throwaway port of the field interpolation and the Boris step to PyTorch, float64, cases A and B.
- Gradients of final (r, pr, phase) with respect to (B scale, V_dee, phi0, initial r, pr), compared with central finite differences.
- Memory and time against 1e3, 1e4, 1e5 steps, with and without `torch.utils.checkpoint`.
- Success: gradients within 1e-6 relative of finite differences; memory at 1e5 steps reported with and without checkpointing.
- Decision gate: the numbers go to the user before any backend is built.

### F. Extraction handoff (Q5), deferred
- Deferred: the CIME deck holds one electrostatic deflector (r = 1.513 m, 250 to 267 deg) and no transfer-line element, so there is nothing after extraction to compare yet. Kept for the step that adds the extraction line.
- LIONS bunch at the deflector exit, converted to Xsuite and Cheetah at a plane s0, tracked through a drift and a quadrupole in each code.
- Success: centroid and rms sizes agree within 1e-6 relative between Xsuite and Cheetah, and with LIONS tracked through the same drift.

## 6. Outputs
- One Parquet table per case: columns above, units in the schema metadata, plus code versions, git commit, machine.
- A report page generated from the tables (Phase 6 docs).
- Phase 3 byte-level guard (`phase3_golden.py`) folded into the performance part.

## 7. Findings from the install checks (single instrument, not yet reproduced elsewhere)
- RF-Track `leapfrog` in a `Volume`, pure B field (ConstantField and UserField give the same numbers): |p| grows by sqrt(1 + (w dt)^2) per step, observed dr/r = +5.636e-2 per turn at 360 steps and +1.380e-2 at 1440, matching (1 + (2 pi / n)^2)^(n/2) - 1 = 5.635e-2 and 1.380e-2. The magnetic kick is applied explicitly, without the Boris rotation. Excluded for cyclotron tracking; `analytic` (exact per step in a constant field) and `rkf45` with tight tolerances are used instead.
- RF-Track `rk2` and `rkf45` are adaptive: the error follows `odeint_epsabs`/`odeint_epsrel` (defaults 1e-3 and 0), not `dt_mm`. Default tolerances give |p| errors of 1e-5 per turn.
- OPAL `CARBONCYCL` values are read with `fscanf("%16lE")`: 16 characters per value at most, so 10 significant digits. A longer number is split silently (first test gave alternating 1 and 0 kG).
- OPAL sets dt = 1 / (RFFREQ x STEPSPERTURN); RFFREQ must be the RF frequency even with no RF map.
- OPAL `-trackOrbit.dat` prints 9 significant digits: metrics below 1e-8 need the H5 phase-space dumps.

## 8. Open checks before running
1. `RF_FieldMap` time dependence and `direction` argument, on a one-gap test.
2. OPAL RF map format (`FM3DDynamic` text vs `FM3DH5Block`) and its axis order.
3. OPAL sign convention: a +10 kG map gives Bz = -1 T (`bz = -bzint`); a proton started at PHIINIT = 0 turns counter-clockwise (checked on case A).
4. OPAL tune-mode output format.
5. RF-Track `UserField` speed: a Python callback per particle per step limits it to single-particle cases; a 3D mesh is needed for bunches.

## 9. Order
A, then D and E in parallel (they need only the case A and B fields), then B, C. F when the extraction line is modelled.

## 10. Case A results and amendments

Run: `uv run --group bench python benchmarks/case_a.py`, output `benchmarks/runs/case_a/case_a.parquet`.
Deterministic single particle, one instrument per code. Observed orders from successive halvings of dt.

| Code, integrator | radius error vs rho | abs(p) change | phase error |
|---|---|---|---|
| LIONS Boris | order 2.00 (6.1e-4 to 2.4e-6) | 0 | order 2.00 (-4.1e-4 to -1.6e-6) |
| OPAL LF2 | round-off (1e-15) | round-off | order 2.00, same values as LIONS |
| OPAL RK4 | order 5.0 (3.6e-8 to 3.2e-14) | order 5.0 | order 4 down to a floor near 5e-11 |
| RF-Track analytic | round-off | round-off | round-off |
| RF-Track rkf45, rk2 | set by the tolerance, not dt | set by the tolerance | set by the tolerance |

Against the criteria of section 5:
- LIONS: order 2, fitted radius equal to rho sqrt(1 + (w dt / 2)^2) within 1e-15, abs(p) exact. Pass.
- OPAL LF2: abs(p) constant within 4e-16. Pass. The radius criterion does not discriminate: with drift-kick-drift, full-step positions are the midpoints of the Boris polygon edges, on the inscribed circle of radius rho exactly (derived after the run).
- OPAL RK4: the radius criterion (order 4) was wrong. RK4 on circular motion scales the amplitude by 1 - (w dt)^6 / 144 per step, so radius and abs(p) errors are order 5; (2 pi)^6 / (144 n^5) = 7.06e-11 at n = 360, observed 7.05e-11 (derived after the run). The phase is order 4 as expected. The floor near 5e-11 at 1440 steps is not explained.
- RF-Track adaptive integrators: errors at or below the requested tolerance. Pass.

Amendments, made after seeing the first results (a new analysis, not a pass/fail change):
- Added the phase metric: the radius metric is blind for LF2. The phase lag of LIONS Boris and OPAL LF2 matches (n - 1)(w dt - 2 atan(w dt / 2)) / (2 pi), 2.54e-5 at 360 steps.
- Harness bug found and fixed before tabulating: `on_step` also runs at step 0, so the first LIONS run recorded n + 1 points; it showed as a first-order phase error.

## 11. Case B, exploratory runs (not the registered analysis)

1 MeV proton, 4-sector field on the 1 deg x 5 mm map, closed orbit at the theta = 0 section.

| Code | 360 steps | 720 | 1440 |
|---|---|---|---|
| LIONS Boris, r_eo [m] | 0.1467722813 | 0.1467734041 | 0.1467638793 |
| OPAL RK4, r_eo [m] | 0.1467728119 | 0.1467604137 | 0.1467624669 |
| LIONS, det(M_rr) - 1 | -6.6e-5 | +2.8e-4 | +2.6e-5 |
| LIONS, det(M_zz) - 1 | +5.4e-9 | +6.9e-10 | -1.8e-10 |

- Neither code converges smoothly in dt; RK4 does not show order 4. The radial one-turn map has det errors of 1e-5 to 3e-4 that do not shrink with dt; the vertical block is clean.
- Competing explanations: (H1) the bilinear interpolation in theta makes Bz continuous with kinks every map column, and the step grid aliases with the map grid; (H2) an error in the benchmark harness; (H3) an error in the codes.
- Discriminating run: Bz = B0 (1 - 0.5 (r - 0.15)), linear in r and constant in theta, so bilinear interpolation is exact and the field is smooth. Closed orbit known in closed form (p / q = |B0| r (1 - 0.5 (r - 0.15))).
  LIONS: r_eo relative error +4.13e-5, +1.03e-5, +2.58e-6 (order 2.00); det(M_rr) - 1 = -1.0e-7, -1.1e-8, -1.2e-9; Qr = 0.96341, Qz = 0.26803, both matching sqrt(1 - n) and sqrt(n) for n = 0.0718.
  This rules out H2 for LIONS and supports H1. H3 is not tested for OPAL on the smooth field yet.
- Consequence: a 1e-6 orbit tolerance is out of reach on a bilinear 1 deg map for any integrator, and one-turn Jacobians straddle kinks. The same holds for gradients in case E.
- Harness bugs found on the way and fixed: the radius array used by `compute_fields` was stale in the harness (it is refreshed by `check_box`); full-step momenta are now a half Boris kick of p(k - 1/2), the exact inverse of the start, so the recorded map is conjugate to the Boris map.
- LIONS bug found on the way, fixed in its own commit: theta = 2 pi exactly for y < 0 with x / r rounding to 1, one column past the map (IndexError on 721-column maps).

## 12. Case B amendment (registered after section 11, before the B0 runs of OPAL and RF-Track)

Section 11 showed that the case B criteria cannot be met on a bilinear map by any integrator. Case B is split in two.
Already seen before this registration: the LIONS B0 values at 1 MeV in section 11. Nothing else.

### B0. Smooth radial gradient (integrator parity, exact references)
- Field: Bz = B0 (1 - k (r - r0)), B0 = -1 T, k = 0.5 /m, r0 = 0.15 m, d2Bz = 0,
  dBz/dx = -k B0 x / r, dBz/dy = -k B0 y / r. Same map grid as case B (5 mm, 1 deg).
  Bz: bilinear interpolation is exact in LIONS and OPAL. Off the plane, OPAL derives Br = z dBz/dr (5-point
  difference, exact here) and Btheta = 0. LIONS interpolates dBz/dx and dBz/dy bilinearly in theta: the
  horizontal field follows the chord of (cos theta, sin theta), off by at most k |B0| |z| dtheta^2 / 8
  (9.5e-8 T at z = 5 mm), which shifts Qz by about 1e-5 relative (derived here).
  RF-Track: `UserField` with the closed form, checked against LIONS `compute_bfield` at random points
  (r in the map, z up to 5 mm) before any run: Bz within 1e-13 T and Bx, By within the chord bound are a
  precondition.
- Energies: 1, 2, 3, 4, 5 MeV protons. Section theta = 0, as case B.
- Exact references, derived here: closed orbit at pr = 0 with p / q = |B0| r (1 - k (r - r0));
  field index n = k r / (1 - k (r - r0)); linear tunes Qr = sqrt(1 - n), Qz = sqrt(n)
  (azimuthal symmetry, field linear in z off the plane, static B so gamma is constant).
  Tunes are compared as fractional parts folded into [0, 0.5], the range of the trace and of the FFT.
- Codes: LIONS Boris, OPAL LF2 and RK4 at 360, 720, 1440 steps per turn;
  RF-Track `rkf45` at tolerances 1e-9 and 1e-12 (absolute and relative), sampled at 1440 points per turn.
- Tunes: trace of the 2x2 blocks of the one-turn Jacobian (central differences, 1e-6 m and 1e-6 p / rho).
  Second instrument for LIONS at 1440 steps: FFT tune of 256 turns at 1e-4 m offsets (section 5).
- Success:
  - fixed-step codes: the closed-orbit error converges at order 2 (Boris, LF2) or at least 3.8 (RK4)
    between 720 and 1440 steps, measured where the error exceeds 1e-10 relative;
  - Richardson extrapolation of r_eo from 720 and 1440 steps with the nominal order: within 1e-7 relative
    of the exact radius (field, conversion and section interpolation, independent of dt);
  - tunes within 1e-4 of the exact values at 1440 steps, for every fixed-step code;
  - RF-Track at 1e-12: r_eo within 1e-8 relative, tunes within 1e-5 (finite differences of an adaptive
    integrator carry noise of order tolerance / step);
  - LIONS FFT tunes within 1e-4 of the Jacobian tunes.
- Dropped for fixed-step codes: the 1e-6 relative orbit criterion of section 5. A second-order integrator
  gives 2.6e-6 relative at 1440 steps here (LIONS in section 11), as in case A.

### B. 4-sector field on the bilinear map (report only, no pass/fail)
- Field, energies, codes and steps as in section 5; tunes from the one-turn Jacobian.
- Reported: r_eo, Qr, Qz, det(M) - 1 per code and step count; spread across codes at fixed steps and
  across steps for each code.
- Map spacing study, LIONS only: 1 deg (361 columns) and 0.5 deg (721 columns, the LIONS maximum), 5 mm radial.
  Prediction of H1: the spread across steps falls with the azimuthal spacing, by about 4 from 1 deg to 0.5 deg
  if it follows the O(dtheta^2) interpolation error. A spread that does not fall weakens H1.

## 13. Case B0 and B results

Run: `uv run --group bench python benchmarks/case_b.py`, output `benchmarks/runs/case_b/case_b.parquet`.
Deterministic single particles, one instrument per code (two for LIONS tunes). Errors relative to the exact values.

### B0
Field parity: Bz within 3.3e-16 T, Bx and By within 8.1e-8 T (bound 9.5e-8 T).

| Code | r_eo error 360 / 720 / 1440 (1 to 5 MeV) | order | Richardson | Qr, Qz error at 1440 |
|---|---|---|---|---|
| LIONS Boris | 3.8e-5 to 4.1e-5 / 9.5e-6 to 1.0e-5 / 2.4e-6 to 2.6e-6 | 2.000 | 1.5e-11 to 2.9e-11 | 2e-8, -3.2e-6 to -5.3e-6 |
| OPAL LF2 | 3.0e-6 to 7.7e-6 / 7.6e-7 to 1.9e-6 / 2.1e-7 to 5.0e-7 | 1.88 to 2.12 | -3.0e-8 to +2.8e-8 | 2e-8, 1.8e-7 to 4.2e-7 |
| OPAL RK4 | constant in steps: +2.2e-8, -2.9e-8, -3.0e-8, +2.8e-8, +1.6e-8 | n/a | same as LF2 | 3e-9, 7e-9 |
| RF-Track rkf45 | -1.8e-12 to -3.3e-12 | | | 6e-10, 6e-13 |

Against the criteria of section 12:
- Order: LIONS and OPAL LF2 pass. OPAL RK4 fails by the letter: its error does not change with the step
  (changes of 3e-10 between 360 and 1440), so no order can be measured above the 1e-10 threshold.
- Richardson within 1e-7: pass for all fixed-step codes.
- Tunes within 1e-4 at 1440 steps: pass for all. LIONS Qz is low by 3e-6 to 5e-6, consistent with the
  chord estimate of section 12.
- RF-Track at 1e-12: r_eo within 3.3e-12, tunes within 6e-10. Pass. Tolerances 1e-9 and 1e-12 give identical
  numbers: the integrator step is capped by the sampling interval (period / 1440), where the local error is
  already below 1e-12.
- LIONS FFT tunes against Jacobian tunes: within 3e-7. Pass.

Unexplained:
- OPAL carries an energy-dependent offset of 1.6e-8 to 3.0e-8 in r_eo, the same in LF2 and RK4, with changing sign.
  It is the RK4 error and the LF2 Richardson residual. Single-precision field storage is ruled out (`BfieldData`
  holds doubles, OPAL-2024.1.0 source). Not tested: input parsing of RINIT/PRINIT, the internal
  dBz/dr evaluation, the momentum conversion.
- OPAL LF2 errors are 5 to 14 times smaller than LIONS Boris and grow with energy, where LIONS decreases.
  Both are second order. Kick-drift-kick and drift-kick-drift differ at O(dt^2); not analysed further.

### B (report only)
- The field does not depend on r, so orbits at different energies are scaled copies: Qr, Qz and det(M) - 1
  are the same at all energies for a given code and step count (observed). The five energies add no information.
- Fractional Qr is 0.005 to 0.007 (Qr close to 1). Spread across steps per code: r_eo 1e-5 to 5e-5 m, Qr 1.1e-3
  to 1.7e-3, Qz 5e-5 (Boris, LF2) and 1e-6 (RK4). Spread across codes at 1440 steps: r_eo 2e-6 to 4e-6 m,
  Qr 5e-4 to 9e-4. det(M) - 1 from 3e-5 to 6e-4, no trend with steps.
- Map spacing (LIONS): the spread across steps at 0.5 deg over the spread at 1 deg is 1.5 for r_eo and 0.93 to 1.03
  for Qr and Qz. H1 predicted about 0.25. The prediction failed: the step-to-step scatter is not set by the
  azimuthal map spacing, so H1 as stated is weakened.
- Rival explanations still open:
  - (H4) conditioning: with Qr close to 1, the closed-orbit Newton matrix M - I has det of about
    (2 pi 0.006)^2 = 1.4e-3, so map errors of 1e-8 become orbit errors of 1e-5, and the trace gives
    dQ = d(trace) / (4 pi sin 2 pi Q), an amplification of about 2 at Q = 0.006.
  - (H1') kinks in theta at any spacing: a second-order integrator crossing derivative jumps; finer maps give more,
    smaller kinks, and the total need not fall.
  - Neither explains the det(M) - 1 of 1e-4 alone; B0 gives 1e-9 with the same harness.
- Discriminating runs, not registered yet: RF-Track on the closed-form 4-sector field (no map, no kinks); the same
  4-sector field with a radial gradient that moves Qr away from 1.

## 14. Case B discriminating runs (registered before running)

Field family: Bz = B0 g(r) f(theta), g = 1 - k (r - r0), f = 1 + eps cos 4 theta, B0 = -1 T, r0 = 0.15 m, eps = 0.3.
Off the plane, first order in z: dBz/dx = B0 (-k f cos theta + 4 eps g sin 4 theta sin theta / r),
dBz/dy = B0 (-k f sin theta - 4 eps g sin 4 theta cos theta / r) (k = 0 gives the case B expressions).
RF-Track uses the closed form; parity with LIONS is checked at the map nodes in theta, random r, z up to 5 mm:
Bz within 1e-13 T (exact there); Bx, By within the linear-in-r interpolation error of g(r) / r,
|z| dr^2 eps |B0| (1 + k r0) / r_i^3 with r_i the inner node of the cell (derived here; found when the first
parity check, written for 1e-13 T on all components, failed at 4.9e-5 T).

### B-smooth: no map
- RF-Track `rkf45`, tolerance 1e-12, closed-form field with k = 0, 1 MeV only (the energies of case B are
  scaled copies), sampling at 360, 720, 1440 points per turn.
- If the map kinks are the source of the scatter (H1'), RF-Track shows none: r_eo spread across sampling
  below 1e-9 m, Qr spread below 1e-6, det(M) - 1 below 1e-8. LIONS-like scatter (above 1e-6 m) rejects H1'
  and points at the harness or at the conditioning alone.

### B-gradient: Qr away from 1
- k = 0.5 /m, same map grid as case B (bilinear in r exact, kinks in theta kept). Energies 1 to 5 MeV.
  LIONS Boris, OPAL LF2 and RK4 at 360, 720, 1440 steps; RF-Track `rkf45` at 1e-12 with the same sampling.
- If the near-integer tune is the amplifier (H4), the r_eo spread across steps times det(M_rr - I) is of the
  same size as in case B at the same energy, so the spread itself falls by the ratio of det(M_rr - I).
  A spread that does not fall rejects H4 as the amplifier. det(M) - 1 is reported; neither hypothesis
  predicts a change.

## 15. Case B discriminating runs, results

Run: `uv run --group bench python benchmarks/case_b2.py`, output `benchmarks/runs/case_b2/case_b2.parquet`.
Field parity at the map nodes: Bz within 1.1e-15 T, Bx and By at 0.92 of the bound.

### B-smooth (RF-Track, closed form, 1 MeV, Qr fractional 0.0060 as in case B)
- Spread across sampling (360, 720, 1440): r_eo 1.9e-9 m, Qr 2.0e-7, Qz 2.8e-10. det(M) - 1: 1.0e-7, 1.2e-8, 9.0e-10.
- Against the registered thresholds: Qr passes; r_eo is 1.9e-9 m against 1e-9 m; det(M) - 1 passes at 1440
  only and falls with the sampling (the section interpolation). The scatter is 4 orders of magnitude below
  LIONS and OPAL on the map (1e-5 m, 1e-3 in Qr) with the same near-integer tune.

### B-gradient (Qr fractional 0.03 to 0.10, det(M_rr - I) 0.037 to 0.40, 26 to 280 times case B)
- RF-Track: r_eo spread 4e-11 to 1.7e-9 m, Qr spread below 1.8e-7, det(M) - 1 at 1440 steps 1e-9 to 1e-8.
- LIONS: r_eo spread 1.4e-6 to 1.7e-5 m (case B: 9.5e-6 to 2.1e-5 m), Qr spread 2e-4 to 2.5e-3.
  OPAL LF2: 7.7e-7 to 1.1e-5 m; RK4: 3.0e-7 to 2.1e-6 m. det(M) - 1 from 8e-6 to 1.4e-2 for the map codes.
- Spread times det(M_rr - I): case B 1.3e-8 to 6.2e-8 m; B-gradient 8e-8 to 6.6e-6 m. Not the same size.
  The spread falls by 1.3 (LIONS) to 55 (OPAL RK4), short of the conditioning ratio.

### Reading
- H4 (near-integer tune as the amplifier) is rejected by its own prediction: the spread does not follow the
  conditioning, and RF-Track is clean at the same near-integer tune.
- The map is the source: same field, same harness, same tunes, clean without the map and scattered with it.
  H1' (kinks at any spacing) is consistent with all observations so far; it is not shown directly.
  A candidate mechanism, not tested: with kinks, the discrete one-turn map is only piecewise C1 in the initial
  conditions (its derivative jumps each time a step point crosses a kink line), so finite-difference Jacobians
  and Newton iterations see noise at every step size and h. This is the same reason AD gradients on the
  bilinear map are expected to disagree with finite differences in case E.
- What would change this: a smooth (bicubic or spline) map interpolation in LIONS giving a clean B-gradient
  run. That changes LIONS results; it stays a decision for the user and O. Kamalou.
- Map bias against the closed form (1 deg, 5 mm), B-gradient: Qz low by 4e-5 to 6e-5 in all three map codes; Qr differs
  by up to 1e-3 at 1440 steps.

## 16. Case D amendment (registered before running)

Changes to section 5 D, from sections 13 to 15:
- Case B is dropped from D: on the bilinear map the discrete one-turn map scatters at 1e-4 in det(M) - 1 for
  reasons unrelated to symplecticity. D uses fields without kinks in the plane: B0.
- The one-turn Jacobian does not discriminate. A symplectic integrator gives an exactly symplectic step map,
  but the section crossing (4-point Lagrange in time) is not symplectic and its error is O(dt^4), the order of
  the RK4 truncation error. At 1440 steps both are below the finite-difference floor (det(M) - 1 of 1e-9 in
  section 13). The B0 Jacobian metrics of section 13 (symplectic_error, det_error) are reported as they are;
  no new short-term runs.
- The discriminating test is long-term: interpolation errors do not feed back into the tracking, integration
  errors that are not symplectic accumulate.

### D2. Long-term invariant in B0
- B0 is azimuthally symmetric in the plane (bilinear map exact in Bz, z = 0), so the canonical angular momentum
  is an exact invariant of the continuous motion. In LIONS units (p = beta gamma), derived here:
  L = x py - y px + (q / (m c)) B0 ((1 + k r0) r^2 / 2 - k r^3 / 3).
- One particle, 1 MeV, z = pz = 0, r = r_eo + 1e-3 m, pr = 0 at the theta = 0 section; 1e4 turns at 360 steps
  per turn (largest step, strongest effect).
- Codes: LIONS Boris, OPAL LF2 and RK4, RF-Track `rkf45` at 1e-12 sampled 360 times per turn.
- Metrics, from L and abs(p) at every sample: secular relative change of L over the 1e4 turns (least-squares
  slope times the duration, over abs(L0)); max abs(L - L0) / abs(L0); relative change of abs(p).
- Classification rule: secular change below 1e-9 is "no drift detected" (round-off bound for 3.6e6 steps).
- Predictions: Boris and LF2 conserve abs(p) exactly in a pure B field and show no drift in L, with a bounded
  oscillation of O((w dt)^2). RK4 shrinks abs(p) by (2 pi)^6 / (144 n^5) = 7.1e-11 per turn (section 10),
  7e-7 over 1e4 turns, and L drifts at the same order. RF-Track drifts at the level of its tolerance.
- Off-plane motion is excluded: the first-order off-plane field of LIONS and OPAL is not divergence free at
  O(z^2), so no exact invariant exists for the model at finite z.

## 17. Case E amendment (registered before running)

- Port: LIONS `compute_bfield` (bilinear map, z terms) and `boris_push` with E = 0, to PyTorch, float64, CPU
  (MPS has no float64). Throwaway, in `benchmarks/case_e.py`, run with `uv run --group bench --with torch`.
- Validation first: the port reproduces LIONS positions and momenta at every step over one turn (360 steps,
  1 MeV, cases A, B0 and B) within 1e-13 relative. Nothing else runs until this holds.
- Outputs: r, pr and the unwrapped phase atan2(y, x) after N steps.
- Parameters: map scale s (multiplies every map array), initial r0, initial pr0.
  V_dee and phi0 are deferred: the spike has no electric field, and the gap model is benchmarked in case C first.
- Gradients: reverse-mode autograd against central finite differences, fixed steps 1e-6 (relative) for s,
  1e-7 m for r0, 1e-7 for pr0; a sweep of 1e-4 to 1e-9 is reported alongside, not used to choose.
- Success: agreement within 1e-6 relative on A and B0, for N = 360 and 3600 steps.
  B (bilinear 4-sector map): report only. Prediction from section 15: autograd and finite differences disagree
  at every step size, at a level comparable to the det(M) - 1 scatter (1e-4).
- Cost: N = 1e3, 1e4, 1e5 steps, 1 and 1000 particles, with and without `torch.utils.checkpoint` over segments of
  sqrt(N) steps. Wall time of forward plus backward and peak RSS, one subprocess per configuration.
- Decision gate: the numbers go to the user before any differentiable backend is built.

## 18. Case D2 results

Run: `uv run --group bench python benchmarks/case_d.py`, output `benchmarks/runs/case_d/case_d.parquet`.
One particle, one run per code, 1e4 turns at 360 steps (3.6e6 steps).

| Code, integrator | L secular change | max abs(L - L0) / abs(L0) | abs(p) change |
|---|---|---|---|
| LIONS Boris | -4.6e-11 | 1.1e-6 | +4.9e-15 |
| OPAL LF2 | -3.6e-12 | 7.8e-8 | -1.0e-13 |
| OPAL RK4 | -1.47e-6 | 1.47e-6 | -7.18e-7 |
| RF-Track rkf45, 1e-12 | +2.0e-8 | 2.0e-8 | +9.7e-9 |

- Boris and LF2: no drift detected (below 1e-9), abs(p) conserved to round-off; L oscillates with a bounded
  amplitude (1.1e-6 for Boris, 7.8e-8 for LF2, the same as over the first 20 turns of the smoke test).
- RK4: abs(p) falls by 7.18e-7, predicted 7.1e-7; L drifts monotonically by 1.5e-6.
- RF-Track: monotonic drift of 2e-8 in L and 1e-8 in abs(p), at the level set by the tolerance accumulated
  over 3.6e6 steps.
- Answer to Q2 on this test: the LIONS Boris step behaves as a symplectic integrator (no secular drift of an
  exact invariant over 1e4 turns, bounded oscillation), like OPAL LF2; RK4 does not. One test, planar, one field,
  one amplitude: consistent with symplectic behaviour, not a proof. The LIONS oscillation of L is 14 times
  the LF2 one; not analysed.

## 19. Case E results (decision gate)

Run: `uv run --group bench --with torch python benchmarks/case_e.py`, outputs `case_e.parquet` and
`case_e_cost.parquet` in `benchmarks/runs/case_e/`. torch 2.15.0, CPU, float64. Run twice: gradients identical,
cost within 10% except the first 1e3-step checkpointed point (6.1 s, then 2.0 s).

### Port
Positions and momenta at every step of one turn equal LIONS within 1.5e-16 (A), 3.0e-16 (B0), 1.0e-15 (B)
relative. Pass.

### Gradients (autograd against central differences, registered steps)

| Case | max relative difference, 360 steps | 3600 steps |
|---|---|---|
| A, uniform | 7.5e-5 | 1.5e-3 |
| B0, smooth gradient | 2.8e-8 | 3.7e-8 |
| B, 4-sector bilinear map | 2.3e-8 | 7.3e-3 |

- B0: pass.
- A: fails the 1e-6 criterion on d(phase)/d(pr0) and d(phase)/d(r0) only. These gradients are small
  (3.0e-5 and 1.1e-3 at 360 steps: the uniform field is isochronous), and the round-off of a central difference
  on a phase of 6.3 rad with h = 1e-7 is of order 3e-8 absolute. The sweep agrees: the difference falls as h grows
  (1.7e-6 at h = 1e-4). Every other entry is within 1e-6. Reading, made after the run: the reference is noisy,
  not the autograd value. A criterion mixing relative and absolute error would pass; it was not registered.
- B: the prediction of section 17 failed in part. At 360 steps autograd and finite differences agree (2.3e-8).
  At 3600 steps they differ by 1e-4 to 7e-3 for h from 1e-7 to 1e-4, and agree within 2e-8 for h of 1e-8 and 1e-9.
  This fits the piecewise-C1 picture of section 15: autograd returns the local derivative of the discrete map;
  finite differences with h above the spacing of the derivative jumps average over them. Example at 3600 steps:
  dr/ds = -0.9824 (autograd) and -0.9897 (h = 1e-6). On a bilinear map, a gradient-based optimiser over
  parameter changes above about 1e-7 sees a derivative that differs from the local one by up to 0.7% after 10 turns.

### Cost (forward plus backward, peak RSS including about 0.4 GB of imports)

| Steps | Particles | no checkpoint | checkpoint every sqrt(N) steps |
|---|---|---|---|
| 1e3 | 1 | 0.44 s, 0.42 GB | 2.0 s, 0.52 GB |
| 1e3 | 1000 | 0.64 s, 0.94 GB | 2.2 s, 0.55 GB |
| 1e4 | 1 | 4.9 s, 1.8 GB | 16 s, 2.2 GB |
| 1e4 | 1000 | 8.1 s, 5.9 GB | 19 s, 2.3 GB |
| 1e5 | 1 | skipped, 18 GB extrapolated | 184 s, 5.6 GB |
| 1e5 | 1000 | skipped, 59 GB extrapolated | 212 s, 6.0 GB |

- The memory guard (skip a configuration extrapolated past half the RAM) was added after registration.
- Without checkpointing, memory is about 140 kB per step for one particle: the autograd graph of about 60 small
  operations per step dominates, not the tensor data. Checkpointing helps with 1000 particles (5.9 to 2.3 GB at
  1e4 steps) and not with one.
- Time per step (forward, recompute, backward) is about 1.8 ms at 1e5 steps, set by Python dispatch of small
  operations; 1000 particles cost little more than one. A compiled step (`torch.compile`, or JAX `lax.scan`)
  is the obvious lever; not measured.
- CIME scale for reference [UNVERIFIED for a CIME deck]: the step counts per turn and turn counts of a CIME run
  set N; at 360 steps per turn, 1e5 steps is about 280 turns.
