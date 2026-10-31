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
- Field: Bz = B0 (1 + eps cos 4 theta), analytic derivatives, d2Bz = 0.
- Metrics: equilibrium orbit radius at theta = 0, Qr and Qz, at 5 energies.
- Method: two-particle method (reference on the equilibrium orbit, one offset by 1e-4 m in r or z), FFT of the offset over 200 turns; OPAL tune mode as a second instrument.
- Success: |Delta Q| < 1e-4 between codes at the converged step; equilibrium radius within 1e-6 relative.

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
