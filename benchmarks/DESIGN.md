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
- Scale: 1e5 steps at 360 steps per turn is about 280 turns. The steps per turn and turn count of the CIME deck
  are not checked here.

## 20. Case E2: speed spike (registered before running)

Question: can the PyTorch port of case E beat plain LIONS throughput by 10 times, and is float32 usable?
- Variants of the case E port: eager (case E), `torch.compile` of the step (fields and Boris push), each in float64
  and float32. CPU on this machine. GPU (IN2P3) once access is set up: same script, reported separately.
- Field: case B map (1 deg), 1 MeV protons, 3600 steps.
- Throughput: particle-steps per second, forward only, at 1, 1e3 and 1e5 particles; forward plus backward
  (gradient of the summed phase with respect to the map scale) at 1e3 particles. Compile time reported apart.
  Median of 3 repetitions after one warm-up.
- Baseline: plain LIONS (`lions.run`, output files included) at the same particle counts, same machine.
- Go criterion for building the backend (section 19 recommendation): forward throughput of the best float64
  variant at least 10 times plain LIONS at 1e3 or 1e5 particles, on the hardware where it would run.
- float32: maximum relative difference of r and absolute difference of phase against float64 after 3600 steps,
  cases B0 and B. Report only: its eligibility is judged later against the calibration needs.

Amendments to section 20, before the reruns:
- First CPU run (2026-10-31): the float32 compiled variant ran at eager speed with no compile time, and the
  compiled-against-eager check gave exactly 0: torch fell back to eager after its recompile limit. Harness bug,
  fixed by resetting the compiler cache per variant. Those numbers are discarded; the float64 rows are kept
  for reference only and the CPU run is repeated.
- Devices: `--device cpu|mps|cuda`. MPS runs the float32 variants only (no float64 on MPS), so it cannot meet the
  go criterion; it is reported as a float32 throughput. The plain LIONS baseline runs on the CPU of the same node.
- Sites: CC-IN2P3 (H100, V100) and Jean-Zay (A100 or H100), from the site clones in the container of
  `environment/`; CRIANN when reachable.
- Accuracy reference, amended before the site and MPS runs: eager float64 on the CPU for every device; each device
  variant (float64 eager and compiled, float32 eager; float32 only on MPS) is compared with it. Reason: a smoke test
  gave MPS float32 against CPU float64 a difference of 5.4e-6 in r after 360 steps, against 1.3e-7 after 3600 steps
  for CPU float32.
- V100 excluded: torch 2.15.0+cu132 has no kernels for sm_70 ("no kernel image is available", CC gpu_v100, job
  62922802). GPU runs: CC H100 and Jean-Zay A100.

## 21. Case G: benchmarks on external machine maps (approved 2026-10-31)

Questions: (G1) how large is the map-interpolation error of LIONS against an exact answer, and does a smooth
interpolation remove the noise floor of section 15? (G2) do LIONS closed orbits and tunes on a real machine agree
with an independent code (FIXPO) as well as OPAL does?

Sources (fetched into the git-ignored `benchmarks/external/`, pinned, sha256 recorded in each run):
- FROM_ORBIT, TRIUMF, GPL-3.0, gitlab.triumf.ca/tplanche/from-orbit at d338aa4. Generates isochronous median-plane
  maps from prescribed closed-orbit shapes and computes the tunes from the orbit geometry. Map generation uses its
  core module only (numpy, scipy, ndsplines); CYCLOPS, closed source, is not needed. Isolated environment (it imports
  `numpy.trapz`, absent from recent numpy).
- PSI Ring 590 MeV, OPALX-project.github.io at 7ed8895, `Examples/Cyclotron/`: map S03AV (TYPE RING: per radius,
  Bz and three azimuthal derivatives; r 1900 to 4700 mm by 20 mm, 45 deg sector by 1/3 deg, 8-fold), FIXPO static
  equilibrium orbits `ic.dat` (E, r, pr; 72 to 530 MeV by 2 MeV) and tunes `refsol.dat`. No licence: internal use
  only, never committed or redistributed without PSI's permission.

### G1. Manufactured solutions (FROM_ORBIT)
- Maps: the published 3-sector (ideal) and 10-sector examples, whose tunes are published against CYCLOPS
  tracking (Cyclotrons 2019 FRB01 and later). A 4-sector, CIME-like variant only afterwards, labelled as our own
  construction, if the published cases behave.
- Truth: the prescribed closed-orbit shapes r(theta) and constant revolution frequency (isochronism) at each
  orbit; FROM_ORBIT geometric tunes. The fidelity of the geometric tunes is checked first with RF-Track on a
  smooth (spline) interpolation of the generated map; that check sets the usable precision of the truth.
- Variables: map spacing in theta (LIONS limit 721 columns), interpolation (LIONS bilinear as it is; bicubic in
  the benchmark harness only, LIONS unchanged), steps per turn.
- Metrics: closed-orbit radius error against r(theta) at theta = 0, revolution-frequency deviation, Qr and Qz
  against the truth, det(M) - 1, spread across steps.
- Codes: LIONS, OPAL (CARBONCYCL), RF-Track (spline interpolation, reference).

### G2. PSI Ring against FIXPO
- Conversion: RING format to a full 360 deg LIONS map (Bz, d2Bz/dr2, dBz/dx, dBz/dy) by 8-fold replication at
  the native 1/3 deg (1081 columns: the LIONS map limit is raised, no resampling), and to OPAL (native RING reader).
  Radial derivatives from the map by a stated scheme; azimuthal derivatives from the file.
- Energies: the 230 FIXPO energies for LIONS; every 20 MeV for OPAL and RF-Track.
- Metrics: equilibrium-orbit r and pr at the FIXPO reference azimuth, Qr and Qz, against FIXPO.
- Criterion: LIONS minus FIXPO no larger than OPAL minus FIXPO plus the interpolation floor measured in G1.
- Unknowns to settle from the files before running: the azimuth of the FIXPO orbits, map units (kG implied by
  the OPAL reader), whether trim coil TC1 is in the FIXPO data.

- PSI files: internal benchmarks only; A. Ghribi asks PSI colleagues for permission to redistribute.

### Out of scope here
Acceleration in the PSI Ring (RF cavities): it needs the gap model of case C. Later.

Case G1 implementation notes (before any tracking):
- The map is the exact field of the orbit family (FROM_ORBIT `def_B`: root finding for a(r, theta), then
  B = brho(a) curv(a, theta) / r_inf), with dBz/dr and dBz/dtheta by central differences (1e-6 m, 1e-6 rad).
  FROM_ORBIT's own `write_field_map` fits a(r, theta) with a smoothing spline: on the 3-sector example at 1 deg it
  departs from the exact field by up to 9.7e-5 (median 5.7e-6, 118 nodes), too large for a reference. Observed once.
- Orbits: indices 1 to 3 of the example grid (3-sector: a = 0.14, 0.24, 0.35; 19, 58, 130 MeV H2+). The centre
  (a = 0.004) is degenerate and the outermost orbit (a = 0.425) bounds the region where the field is defined, so
  no map can surround it.
- Map: r 0.40 to 2.00 m by 1 cm, 360 points per 120 deg sector (1/3 deg); coarser maps are subsets of its nodes.
  LIONS takes -Bz (counterclockwise travel); d2Bz = 0, as in OPAL (no effect on linear tunes).
- Environment: `benchmarks/g1_maps.py` is a uv script with pinned numpy 1.26.4 and scipy 1.13.1 (FROM_ORBIT uses
  `numpy.trapz`).

## 22. Case E2 results (CPU, Apple GPU, A100; H100 pending)

Runs: M3 CPU and MPS locally (`benchmarks/runs/case_e2_cpu`, `case_e2_mps`), A100 on Jean-Zay (job 1791641, image
lions-sim-58392b29f72a, torch 2.15.0+cu132, code 501c88f). Median of 3 after a warm-up; plain LIONS timed on the CPU
of the same node. Throughput ratio to plain LIONS on the same node:

| Hardware | best float64, 1e3 particles | best float64, 1e5 particles | float32 compiled, 1e5 | compiled forward+backward, 1e3 |
|---|---|---|---|---|
| M3 CPU | 6.7 | 15.7 | 15.1 | 0.71 s |
| Apple GPU (MPS) | no float64 | no float64 | 150 | fails (InductorError) |
| A100 (Jean-Zay) | 5.0 | 267 | 267 | 2.49 s |

- Go criterion (section 20): met. Compiled float64 reaches 267 times plain LIONS at 1e5 particles on the A100
  (6.6e8 particle-steps/s), and 15.7 times on the M3 CPU. At 1e3 particles no device reaches 10.
- A100 wall time is 0.52 to 0.55 s for 3600 steps whatever the particle count (1 to 1e5): launch-bound, so the
  throughput keeps rising with the batch; 1e5 is not the ceiling.
- float32 gives no speed gain on the A100 (same wall time) and loses accuracy: r 1.2e-5, phase 1.8e-4 rad after
  10 turns against CPU float64, as on MPS; CPU float32 gives 1.3e-7. GPU float64 against CPU float64: 2e-14 in r,
  1.6e-13 rad in phase, compiled or not.
- Gradients at 1e3 particles are faster on the CPU (0.71 s) than on the A100 (2.49 s): small batches are latency-bound.
- The plain LIONS baseline is 2.3 to 2.6 times slower on the Jean-Zay node than on the M3.

### G1 step 1: precision of the FROM_ORBIT references (registered before running)
- Instrument: RF-Track `rkf45`, tolerance 1e-12, on a quintic spline (r, theta; periodic in theta) of the exact
  master map, with Br = z dBz/dr and Btheta = (z / r) dBz/dtheta from splines of the exact derivative grids.
  Sampled 1440 and 2880 times per turn.
- Interpolation error of the spline, estimated before tracking: the spline built on every second node, evaluated
  at the skipped nodes against the master; the master's own spline error is below that estimate.
- Metrics per orbit (3-sector, a = 0.14, 0.24, 0.35): closed-orbit r and pr at theta = 0 against the prescribed
  orbit; one-turn time against 1 / f_rev (isochronism); Qr, Qz (one-turn Jacobian, folded) against FROM_ORBIT.
- Reading: no pass/fail. The RF-Track to FROM_ORBIT differences, if larger than the interpolation and tolerance
  effects, are the precision of the references and become the floor for every code in G1.

G1 step 1 results (3-sector, `benchmarks/runs/case_g1/step1/case_g1_step1.parquet`). Spline on every second node:
8.7e-7 at the skipped nodes (full-grid spline about 64 times smaller for a quintic). Identical at 1440 and 2880
samples per turn:

| a (Ek, H2+) | r_eo error | pr_eo error / p | turn time error | Qr RF-Track - FROM_ORBIT | Qz RF-Track - FROM_ORBIT |
|---|---|---|---|---|---|
| 0.14 (19 MeV) | -6e-14 | -4e-13 | -7e-14 | -1.30e-5 | +1.36e-6 |
| 0.24 (58 MeV) | -2.0e-9 | +1.9e-9 | -2.7e-10 | -2.05e-5 | +1.69e-6 |
| 0.35 (130 MeV) | +2e-13 | -1e-13 | +5e-14 | +6.3e-7 | +7.5e-7 |

- The prescribed closed orbits and isochronism are reproduced to 1e-13 (2e-9 at a = 0.24): the map, the LIONS-frame
  conversion (sign, direction, theta origin) and the references agree.
- FROM_ORBIT tunes are good to 2e-5 (Qr) and 2e-6 (Qz); the difference does not move with sampling, so it belongs
  to its tune formula. G1 references: the prescribed orbits and isochronism, and the RF-Track tunes on the spline.

### G1 step 2: LIONS on the FROM_ORBIT 3-sector map (registered before running)
- Map: the master nodes at 1/3, 1 and 4 deg (every 1, 3, 12 columns of the master), r 0.40 to 2.00 m by 1 cm,
  360 deg with the endpoint; Bz = -B(FROM_ORBIT), dBz/dx and dBz/dy from the exact derivative grids, d2Bz = 0.
- Particle: the FROM_ORBIT mass and charge (H2+, 2 x 938.272 MeV, charge 1), given to LIONS as an atomic mass
  such that AMU (A - m_e q) c^2 equals it; RF frequency = f_rev of the references, harmonic 1.
- Steps per turn: 360, 720, 1440. Orbits a = 0.14, 0.24, 0.35.
- Metrics: closed-orbit r and pr at theta = 0 against the prescribed orbit; one-turn time against 1 / f_rev;
  Qr, Qz against the RF-Track tunes of step 1 (2880 samples); det(M) - 1; spread across steps.
- Prediction (from sections 13 to 15, stated before running): with the bilinear map, a step-dependent scatter;
  at 1/3 deg, closed-orbit errors of 1e-6 or more and tune errors of 1e-4 to 1e-3. Much better than that means the
  case B mechanism needs a second look.
- Added before the full step 2 run, after the 4 deg map made Newton diverge on the first orbit (r error -10%,
  residual 2e-3 m, then a particle left the map): a configuration whose closed-orbit residual exceeds 1e-9 m, or
  that raises, is recorded as failed with its error, and the run continues.

G1 step 2 results (`benchmarks/runs/case_g1/step2/case_g1_step2.parquet`; the 1/3 and 1 deg rows reproduce the
first, interrupted run exactly). Errors against the references, 360 / 720 / 1440 steps per turn, over the three orbits:

| Map | r_eo | pr_eo / p | turn time | Qr | Qz | det(M) - 1 |
|---|---|---|---|---|---|---|
| 1/3 deg | 4.8e-5 to 5.8e-5 / 1.5e-5 to 1.7e-5 / 4.1e-6 to 6.9e-6 | -5.5e-5 to -6.7e-5 / -2.0e-5 to -2.4e-5 / -1.0e-5 to -1.3e-5 | 2.5e-6 to 3.7e-6 at 1440 | 3.5e-5 to 5.4e-4, no trend in steps | -2.0e-5 to -7.5e-5 | 9e-7 to 3.7e-3 |
| 1 deg | 5.9e-5 to 6.9e-5 / 1.5e-5 to 3.5e-5 / 4.8e-6 to 2.6e-5 | -9.3e-5 to -1.1e-4 / about -6e-5 / about -5e-5 | 2.2e-6 to 1.5e-5 at 1440 | 2.2e-4 to 1.6e-3 | -1.1e-4 to -4.2e-4 | 1.5e-4 to 1.6e-2 |
| 4 deg | closed orbit not found (9 of 9: residual 7e-4 to 9e-3 m, or particle lost) | | | | | |

- Prediction of step 2 confirmed: closed-orbit errors of 4e-6 to 7e-5, tune errors up to 1.6e-3, det(M) - 1 erratic.
- The closed orbit converges more slowly than second order at 1/3 deg (ratios 3.3 and 2.4 per halving) and stalls at
  1 deg (pr error near -5e-5 at all steps): a map floor that grows with the azimuthal spacing.
- Isochronism: the turn time is off by 2e-6 to 1.5e-5 at 1440 steps on an exactly isochronous field. Over 500 turns
  at harmonic h this is an RF phase slip of order 2 pi h x 500 x 1e-5 rad (about 1.8 deg per unit h): relevant for
  phase studies on CIME maps [estimate, not simulated].
- Not yet separated: integrator error and interpolation error. Step 3 does that with the same Boris step on a
  smooth interpolation.

### G1 step 3: integrator against interpolation (registered before running)
- Instrument: a numpy Boris tracker with the LIONS step (same rotation formulas, positions in m, momenta beta gamma,
  momentum at half steps, start by a backward half kick, full-step momenta by a forward half kick), with a
  pluggable field.
- 3a, validation: with the LIONS bilinear field (same nodes, dBz/dx and dBz/dy maps, same weights) it must
  reproduce LionsTracker on the 1 deg map: one turn from the reference orbit at 360 steps within 1e-12 relative
  (theta from atan2 here, arccos in LIONS). Nothing else runs until this holds.
- 3b: the same tracker and node sets (1/3 and 1 deg) with bicubic and quintic splines (periodic in theta) of
  Bz, dBz/dr, dBz/dtheta; steps 360, 720, 1440; metrics and references as step 2.
- Reading: if the smooth fields converge at second order to the references (as B0 did in section 13), the step 2
  floor and scatter belong to the bilinear interpolation, and the bicubic numbers measure what a smooth LIONS
  interpolation would give. If the scatter stays, the cause is elsewhere.

G1 step 3 results (`benchmarks/runs/case_g1/step3/case_g1_step3.parquet`, code 197648c).
- 3a passed: the numpy Boris tracker on the LIONS bilinear field reproduces LionsTracker over one turn to 4.7e-15.
- 3b: absolute errors against the references, range over the three orbits and both splines (bicubic, quintic),
  360 / 720 / 1440 steps per turn. The 1/3 and 1 deg maps give the same numbers except det(M) - 1:

| Metric | 360 | 720 | 1440 | ratio per halving |
|---|---|---|---|---|
| r_eo | 4.5e-5 to 5.7e-5 | 1.1e-5 to 1.4e-5 | 2.8e-6 to 3.6e-6 | 4.00 |
| pr_eo / p | 4.7e-5 to 5.9e-5 | 1.2e-5 to 1.5e-5 | 2.9e-6 to 3.7e-6 | 4.00 |
| turn time | 3.5e-5 to 3.9e-5 | 8.7e-6 to 9.8e-6 | 2.2e-6 to 2.4e-6 | 4.00 |
| Qr | 4.5e-6 to 1.4e-5 | 1.2e-6 to 3.5e-6 | 2.9e-7 to 8.7e-7 | 3.7 to 4.0 |
| Qz | 2.2e-6 to 5.5e-5 | 5.3e-7 to 1.4e-5 | 1.2e-7 to 3.6e-6 | 3.9 to 4.5 |
| det(M) - 1, 1/3 deg | 4.6e-7 to 7.6e-7 | 7.5e-8 to 1.0e-7 | 1.1e-8 to 2.0e-8 | 3.2 to 370 (both maps) |
| det(M) - 1, 1 deg | 5.8e-7 to 2.3e-6 | 7.1e-8 to 2.9e-7 | 8.0e-10 to 1.5e-8 | |

- Bicubic and quintic agree to 4e-4 relative in orbit and turn time, and to 8 to 14% in the tunes (differences of
  order 1e-7). The remaining error is that of the Boris step: the spline interpolation error is below it.
- Against step 2: the same step on the same nodes, with the bilinear field replaced by a spline, removes the floor
  and the scatter. Tune errors fall from up to 1.6e-3 to at most 5.5e-5 (360 steps), det(M) - 1 from up to 1.6e-2
  to at most 2.3e-6, and the 1 deg map behaves as the 1/3 deg one. The interpolation is the only change, so on this
  map the step 2 floor and scatter are attributed to the bilinear interpolation. At 360 steps and 1/3 deg the
  closed-orbit errors of step 2 and step 3 are close (4.8e-5 to 5.8e-5 against 4.5e-5 to 5.7e-5): there the
  integrator dominates.
- det(M) - 1 does not converge at a clean order and differs by up to a factor 3 between the two splines. Not
  separated: the finite-difference Jacobian (1e-6 m), spline derivative noise, and the Boris step being exactly
  volume-preserving only in a uniform field. It stays at or below 2.3e-6.
- Isochronism with a smooth field is limited by the step: 2.2e-6 at 1440 steps. The phase-slip estimate of step 2
  therefore holds in part for any field interpolation at that step count [estimate, not simulated].
- Scope: one map (3-sector, ideal), three orbits, one turn. Changing the LIONS interpolation changes results and
  breaks Fortran bit-identity: decision for A. Ghribi and O. Kamalou.

### G1 10-sector: steps 1 to 3 on the second published example (registered before running)
- Example: `example/10-sector/10sector-2h.py` (N = 10, orbit shape with a second harmonic, protons, r_inf = 28 m),
  run as is up to its field-map section. Orbits: the interior indices 1 to 4 (a = 0.8401 to 0.9094, 932 to
  1630 MeV, r 22.8 to 26.3 m); indices 0 and 5 bound the family. FROM_ORBIT tunes: Qr 2.78, Qz 2.72 to 2.73
  (fractional 0.22 and 0.27 after folding, away from 0 and 0.5).
- Master: r 22.40 to 26.60 m by 2 cm (211 nodes, LIONS limit 250), 108 points per 36 deg sector (1/3 deg; 1081
  columns on the full turn, the LIONS limit). Exploratory pre-check, before registering: root finding converges
  without warning from 22.0 to 27.3 m. No rectangular map lies inside the orbit family (largest r of orbit 0:
  23.69 m, smallest r of orbit 5: 25.40 m), so the map corners use FROM_ORBIT's spline extrapolation of the orbit
  parameters beyond orbits 0 and 5. The tracked orbits lie inside the family, so the references hold on them. The
  map edges are at least 0.35 m from any tracked orbit; the extrapolated field degrades beyond 26.8 m (peak 1.12 T
  at 26.5 m, 0.93 T at 27.0 m).
- Steps 1, 2 and 3 as registered for the 3-sector map, unchanged: instruments, metrics, references, 360 / 720 /
  1440 steps per turn, 1/3, 1 and 4 deg maps, failure rule. Outputs in `benchmarks/runs/case_g1/10sector/`.
- Differences that change the expected numbers: the field period is 36 deg with an 18 deg harmonic, so 360 steps
  per turn give 18 steps per harmonic period (120 per period on the 3-sector map), and the 4 deg map has 4.5 nodes
  per harmonic period. The radial node step is 8e-4 of the radius (1e-2 on the 3-sector map).
- Prediction: step 1 reproduces the prescribed orbits to 1e-9 or better; step 2 errors larger than on the 3-sector
  map at equal steps per turn, and the 4 deg map fails; step 3 converges at second order with both splines, with
  larger coefficients than on the 3-sector map. If step 3 does not converge at second order, the 3-sector reading
  (the floor belongs to the bilinear interpolation) does not carry over to this map.

G1 10-sector results (`benchmarks/runs/case_g1/10sector/step{1,2,3}`, code cd316ad; master 211 x 108, no root
warnings). Absolute errors against the references, range over the four orbits.

Step 1 (RF-Track on the quintic spline):
- Closed orbit and isochronism: 5e-12 to 6e-11 on the outer orbits (a = 0.8401, 0.9094), 5.6e-8 to 2.3e-7 on the
  middle ones (a = 0.8632, 0.8863), identical at 1440 and 2880 samples: not integration. Prediction (1e-9 or better)
  missed on two orbits of four.
- Spline error, every-second-node estimate: 9.1e-4 relative to the local Bz, which is ill-posed here (Bz crosses
  zero in the valleys, -0.48 to +1.12 T); 1.3e-5 T absolute, 1.0e-5 of max |Bz| (3-sector: 1.7e-7). Consistent with
  the middle-orbit errors [exploratory check, after the run; the reason for the difference between orbits is not
  established].
- Tunes RF-Track minus FROM_ORBIT: 8e-6 to 4e-5 (3-sector: up to 2e-5). det(M) - 1: 5e-8 to 1.3e-6 at 2880 samples.

Step 2 (LIONS, bilinear) and step 3 (numpy Boris, bicubic and quintic), 360 / 720 / 1440 steps per turn:

| Metric | step 2, 1/3 deg | step 3, 1/3 deg | step 3 ratio per halving | step 2, 1 deg | step 3, 1 deg |
|---|---|---|---|---|---|
| r_eo | 8.7e-5 to 1.3e-4 / 1.9e-5 to 4.8e-5 / 3.7e-6 to 2.8e-5 | 8.4e-5 to 1.6e-4 / 2.1e-5 to 3.9e-5 / 5.2e-6 to 9.7e-6 | 4.0 | 7.1e-5 to 4.8e-4 / 2.1e-6 to 7.7e-5 / 1.6e-5 to 1.7e-4 | 4.5e-6 to 1.1e-5 at 1440 |
| turn time | 7.2e-5 to 9.2e-5 / 1.2e-5 to 3.5e-5 / 1.9e-7 to 2.4e-5 | 7.0e-5 to 1.1e-4 / 1.8e-5 to 2.6e-5 / 4.4e-6 to 6.7e-6 | 4.0 to 4.1 | 4.5e-5 to 1.8e-4 at 1440 | 2.4e-6 to 9.4e-6 at 1440 |
| Qr | 2.9e-3 to 1.4e-2 / 5.6e-4 to 6.3e-3 / 7.5e-4 to 3.2e-3 | 4.4e-4 to 3.6e-3 / 1.1e-4 to 9.1e-4 / 2.8e-5 to 2.3e-4 | 3.7 to 4.2 | up to 1.1e-1 / 6.1e-3 / 9.5e-3 | 2.2e-5 to 4.1e-4 at 1440 |
| Qz | 9.7e-4 to 4.2e-3 / 4.0e-4 to 1.4e-3 / 3.9e-5 to 1.8e-3 | 6.8e-5 to 1.5e-3 / 1.5e-5 to 3.6e-4 / 9.1e-7 to 9.1e-5 | 4 on three orbits; sign change near 5e-6 on a = 0.8632 | up to 1.7e-2 at all steps | 2.0e-5 to 4.4e-4 at 1440 |
| det(M) - 1 | 5.2e-2 to 2.4e-1 / 7.3e-2 to 1.1e-1 / 7.0e-3 to 5.3e-2 | 1.4e-4 to 5.3e-4 / 4.0e-6 to 9.8e-5 / 1.3e-6 to 2.1e-5 | 3.1 to 36 | up to 5.0 | up to 5.6e-5 at 1440 |

- 3a passed: 2.4e-14. 4 deg map (step 2): 11 of 12 configurations failed (closed orbit not converged, or the
  particle left the map); the one that converged has Qr and Qz errors of 0.22 and 0.27.
- 1/3 deg: the splines converge at second order in orbit, turn time and Qr, and agree with each other to 4e-3 in
  orbit and 11% in Qr. Qz on a = 0.8632 changes sign between 720 and 1440 steps at the 1e-6 to 1e-5 level, near the
  precision of the step 1 references for the middle orbits. Prediction held.
- Same step, same nodes, bilinear against spline at 1440 steps: Qr 7.5e-4 to 3.2e-3 against 2.8e-5 to 2.3e-4,
  det(M) - 1 7e-3 to 5.3e-2 against 1.3e-6 to 2.1e-5, and step 2 does not converge in the tunes. The bilinear floor
  of the 3-sector map is reproduced on this map, larger.
- 1 deg (18 nodes per 18 deg harmonic period): the spline results no longer converge cleanly (ratios 3.1 to 8.7 in
  turn time, irregular in the tunes) and bicubic and quintic differ by up to a factor 5 in the tunes; Qz stays near
  4e-4 at 1440 steps. The interpolation error is no longer negligible. Prediction (second order with both splines)
  failed at 1 deg. The 3-sector reading "bicubic at 1 deg removes the floor" does not carry over: what matters is
  the node count per period of the field harmonics, and 1/3 deg (54 nodes per period) is needed here.
- The Boris step error is larger than on the 3-sector map at equal steps per turn (r_eo 3 times, Qr up to 250
  times at 360 steps): at 1440 steps the integrator alone leaves Qr errors up to 2.3e-4 and turn-time errors of
  4e-6 to 7e-6.
- Scope: one map, four orbits, one turn each; no acceleration.

### G2 step 1: conventions of the PSI Ring files (registered before running)
Settled from the files and the OPAL 2024.1.0 source (`Cyclotron::getFieldFromFile_Ring`, `getdiffs`, `TrimCoil`,
`TrimCoilMirrored`), before any tracking:
- `bfield.dat` (S03AV): header 1900 mm, 20 mm, 0 deg, -3 (1/3 deg); 141 radii, 135 azimuths per 45 deg sector,
  8-fold. A first record of 1089 values is skipped by OPAL ("useless"); then one block per radius: Bz, dBz/dtheta,
  d2Bz/dtheta2, d3Bz/dtheta3 (135 values each). Values about 6 kG at 2.1 m: kG, consistent with the average field
  Brho / r = 0.58 T of a 72 MeV proton there (order-of-magnitude check).
- OPAL does not use the file derivatives: `getdiffs` recomputes all theta derivatives from Bz by 5-point
  differences, with stencils clamped inside the 45 deg sector (one-sided near its edges, no periodic wrap), and
  dBz/dr by 5-point differences in r.
- `ic.dat`: E [MeV], r [mm], pr; the OPAL example (`cyclotron1.in`) starts each orbit at PHIINIT = 0 with
  PRINIT = pr / 1000 in units of m0 c, while the OPAL manual names the column P_r [mrad]. The two readings differ
  by beta gamma (20% at 530 MeV, where pr = 212). FIXPO orbits end at r = 4300 mm.
- Trim coil TC1 (`cyclotron1.in`, TYPE PSI-BFIELD-MIRRORED, RMIN 4350, RMAX 4470 mm, SLPTC 0.6): OPAL 2024.1
  converts RMIN and RMAX to m and passes r in m, so SLPTC acts per metre and the profile is nearly flat: from the
  formula, about +5.9e-3 kG at every radius from 2 to 4.5 m (1e-3 of the field). With r in mm the field would stay
  within 4350 to 4470 mm. Possible units problem in OPAL [from the source, not run].

Runs:
- 1a: file dBz/dtheta against a periodic 5-point difference of Bz (per deg and per rad): max and median absolute
  difference, to fix the units of the file derivatives.
- 1b: closed orbits at theta = 0 with the numpy Boris tracker (validated in G1 step 3a) on a quintic spline of Bz
  (kG to T), periodic in theta, 1440 steps per turn, at 72, 150, 300, 450 and 530 MeV (rows of `ic.dat`), both
  directions of travel. Compare r and pr with `ic.dat` under both pr readings.
- 1c: OPAL (container 2024.1.0) closed-orbit finder at 72 and 530 MeV, with and without TC1: shift of the
  closed-orbit radius and tunes.
- Success: one direction and one pr reading match `ic.dat` at the five energies within 1 mm in r and 1% in pr,
  and the other combinations miss by at least ten times that. Otherwise stop and report. 1c reads as: a shift at
  72 MeV (far from 4350 mm) confirms the flat TC1 profile; none refutes it. G2 tracking (LIONS, OPAL, RF-Track
  against FIXPO) is registered after step 1, without TC1 if 1c confirms the flat profile.

G2 step 1 results (`benchmarks/runs/case_g2/run2/step1/case_g2_step1.parquet`, code 7c57e70; a first run stopped at
1c because the OPAL input lacked `OPTION, VERSION`, its 1a and 1b output was identical).
- 1a: the file dBz/dtheta is per radian: median difference to a periodic 5-point difference of Bz 7.8e-3 kG/rad
  (6.0 if read per degree), max 2.4 kG/rad at the sector edges against a peak of 568 kG/rad.
- 1b (numpy Boris, quintic spline, 1440 steps): closed orbits at the `ic.dat` radii to 0.004 to 0.065 mm at the five
  energies. pr read as 1e-3 m0c: within 4e-7 to 2.5e-4 relative; read as mrad: off by 8.3% to 150%. The map azimuth
  increases along the motion: the reversed direction gives the same r and the opposite pr (mirror). Settled:
  `ic.dat` pr is in 1e-3 m0c at theta = 0, travel with increasing map azimuth. The registered criterion asked the
  other readings to miss by ten times the tolerance (10%); the mrad reading misses by 8.3% at 450 MeV, so the
  criterion is not met as written, while the match is 300 times closer than the alternative at every energy.
- 1c (OPAL closed-orbit finder): without TC1, r = 2.1315 m and pr = -2.392e-4 at 72 MeV, 4.30024 m and 0.212152 at
  530 MeV (`tunes.dat` gives r to 0.01 mm at most). TC1 changes nothing at 72 MeV and moves r by +0.66 mm, Qr by
  +3e-5 and Qz by +8e-4 at 530 MeV (r = 4.30 m, 5 cm inside RMIN). The flat TC1 profile read from the source is
  refuted for the closed-orbit-finder path: the field stays near the coil there. Withdrawn. Not checked: the
  tracking path, and the size of the 530 MeV shift.
- OPAL minus FIXPO (`refsol.dat`), no TC1: Qr +5.3e-3 and Qz -6.7e-3 at 72 MeV, +4.8e-3 and -7.6e-3 at 530 MeV.
- OPAL 2024.1 needs two changes to the published example input: FMLOWE, FMHIGHE in GeV (the example gives MeV; all
  energies then fail), and an element name other than `Ring` (protected keyword).

### G2 step 2: LIONS and the spline tracker on the PSI Ring map (registered before running, approved 2026-11-01)
Amendments to section 21, decided before any G2 tracking:
- TC1 off in every code. FIXPO's trim-coil content is unknown, the FIXPO orbits end 5 cm inside TC1, and its only
  measured effect (step 1c) is the sensitivity quoted with the results: +0.66 mm in r and +8e-4 in Qz at 530 MeV.
- Criterion: |LIONS - FIXPO| <= |OPAL - FIXPO| + |LIONS - spline| per metric and energy, where the last term is
  the interpolation floor measured on this map in this step (same nodes, same Boris step), not the G1 floors (which
  differ by a factor 6 between the two G1 maps). Evaluated in step 3, once OPAL runs.

Step 2:
- LIONS map: the 45 deg sector replicated 8 times at the native 1/3 deg (1081 columns, endpoint included), r 1.90 to
  4.70 m by 2 cm (141). Bz in T with the LIONS sign (-0.1 x file, travel with increasing map azimuth, step 1b);
  dBz/dtheta from the file (per rad, step 1a); dBz/dr by 5-point differences in r (OPAL's scheme, stencils clamped
  at the map edges); dBz/dx, dBz/dy from these; d2Bz = 0.
- Spline tracker: the numpy Boris step (G1 step 3) with quintic splines, periodic per sector, of the same Bz,
  dBz/dr and dBz/dtheta node values.
- Energies: the 230 `ic.dat` rows (72 to 530 MeV). 1440 steps per turn; RF 50.65 MHz, harmonic 6; protons.
- Metrics: closed orbit r and pr at theta = 0 against `ic.dat` (r in mm, pr in 1e-3 m0c); Qr, Qz (one-turn
  Jacobian, folded into [0, 0.5]) against `refsol.dat` folded the same way; det(M) - 1; LIONS minus spline for each.
  Tunes between 386 and 430 MeV (Qr crosses 1.5, folded Qr near 0.5) are reported but flagged: the trace formula is
  ill-conditioned there.
- Failure rule as G1 step 2 (residual above 1e-9 m, or an exception, recorded as failed).
- Prediction: spline tracker within 0.1 mm of `ic.dat` (step 1b gave 0.065 mm at 1440 steps without vertical field);
  LIONS minus spline of order 1e-5 to 1e-4 relative in r and 1e-4 to 1e-3 in the tunes (bilinear at 1/3 deg on a
  field with 8-fold spiral sectors, between the two G1 maps); OPAL minus FIXPO tunes about 5e-3 (step 1c), so the
  step 3 criterion is expected to hold unless LIONS has a defect beyond interpolation.

### G2 step 3: OPAL and RF-Track, and the criterion (registered before running)
- Energies: every tenth `ic.dat` row (72 to 512 MeV by 20 MeV) and 530 MeV: 24 energies.
- OPAL 2024.1 closed-orbit finder, RING reader, no trim coil, 1440 steps (as step 1c), one run per energy at the
  exact `ic.dat` energy: r, pr at theta = 0 (`tunes.dat`, r to 0.01 mm at most), Qr, Qz.
- RF-Track `rkf45`, tolerance 1e-12, 2880 samples per turn (G1 step 1) on the step 2 quintic spline: closed orbit,
  Qr, Qz, det(M) - 1.
- Criterion (section 21 as amended in step 2): at each of the 24 energies and for r, pr, Qr, Qz,
  |LIONS - FIXPO| <= |OPAL - FIXPO| + |LIONS - spline|. Reported per metric as the number of energies where it holds,
  the half-integer band (386 to 430 MeV) counted separately.
- Reading of RF-Track: a second smooth-field instrument; RF-Track minus spline tracker measures the Boris step error
  at 1440 steps on this map.

### G2 step 2b: LIONS one-turn determinant, discriminating run (registered before running)
Observation from step 2, before its write-up: LIONS det(M) - 1 = 0.19 (median over energies, up to 0.42) and Qr
error 1.2e-2 (median), against 2.6e-5 and 9e-5 for the spline tracker on the same node values.
Hypotheses: H1, the bilinear field model on the steep spiral edges; H2, a LIONS-specific defect (map conversion or
LIONS itself); H3, an artefact of the finite-difference Jacobian on a piecewise-bilinear field.
- Runs at 152, 292 and 452 MeV (`ic.dat` rows 40, 110, 190): LIONS, the numpy Boris on the LIONS bilinear field
  (G1 step 3a: reproduces LIONS to 1e-14), the spline tracker. Metrics: det of the r block, of the z block, of M;
  Qr, Qz. LIONS also with Jacobian steps 1e-7, 1e-6, 1e-5 m.
- Reading: numpy bilinear equal to LIONS supports H1 and refutes H2; a determinant that moves with the Jacobian
  step by more than 10% supports H3.

G2 step 2b results (`benchmarks/runs/case_g2/step2b/case_g2_step2b.parquet`): numpy Boris on the bilinear field
equals LIONS to all printed digits (H2 refuted); the determinant moves by at most 2% between Jacobian steps 1e-7 and
1e-5 m (H3 refuted). The effect is in the radial block only: det(M_r) - 1 = +0.25, +0.41, -0.17 at 152, 292,
452 MeV; z block 3e-10 to 7e-7; spline 2e-6 to 3e-4. H1 holds and splits: H1a, the bilinear field itself; H1b, the
Boris step (volume-preserving, not symplectic) on a field whose gradient jumps at every cell edge. The exact flow in
any continuous median-plane Bz preserves area, so under H1b an accurate integrator gives det near 1.

### G2 step 2c: integrator against field model on the bilinear field (registered before running)
- RF-Track `rkf45`, tolerance 1e-12, 2880 samples, on the LIONS bilinear field (numpy, G1 step 3a), at 152, 292,
  452 MeV; metrics as step 2b.
- Reading: det(M_r) near 1 (below 1e-3) with RF-Track supports H1b (the non-symplectic Boris step on a kinked
  field); det as large as LIONS supports a field-model effect beyond the exact-flow argument.
