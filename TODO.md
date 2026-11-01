# To do

## Benchmarks
- [ ] Contact the OPALX authors about the radial modulus sensitivity (manual, `physics/linear-transfer-maps.qmd`, open there).
  Case H and G3 bear on it. Later: goes with the PSI permission request, since the results use their PSI map.
- [ ] Report the RF-Track leapfrog bug.
- [ ] CIME checks: needs `cime700rec.1804.pbob1.lio.dat`.
- [ ] G3c, space charge at CIME currents: needs the CIME map and a CIME current and bunch. Not registered.
- [ ] Case C: check the RF time convention of each code on a one-gap test.

## Untested hypotheses
- [ ] Autograd Jacobian of the PyTorch bilinear port against the finite-difference one (DESIGN.md section 24, erratum 4).
- [ ] Convergence order of det(M) with the step count, beyond 2880 steps per turn.
- [ ] Fixed point and det(M_r) against the launch step phase, on one orbit (DESIGN.md section 25).
- [ ] Why LIONS det(M_r) is near the inverse of the OPALX one at equal field (r = +0.94, DESIGN.md section 25, step 2d).
- [ ] Why a closed-orbit launch reaches 3 to 4 mm at 252, 312 and 332 MeV and stays below 0.2 mm elsewhere (step 2d).
- [ ] Beam size and emittance on the bilinear field against Hermite: over 1000 turns or more, on a bunch, with acceleration. Amplitudes at
  0.1 to 1 mm change by 0.24 to 10 times in 100 turns (G3 step 2). Before any surrogate is trained on LIONS bilinear data.

## Differentiable backend
- [ ] Sim-to-real: needs CIME parameters and probe data from A. Ghribi.

## Training game (`docs/training/SPEC.md`)
- [ ] Operator tuning procedure: CIME operations team, not written yet (D4).
- [ ] Sim-to-real magnitudes from CIME data (D3).

## Decisions for A. Ghribi and O. Kamalou
- [ ] Python adds a vertical spread for `espace='6d'`; the Fortran only for `'3d'` (`lions_2.0_ep.f90` line 498).
- [ ] `mass_c` is indexed by family number in both codes.
- [ ] Second-harmonic dBz/dy sign (README, known issue 10).
- [ ] Push `lecture-notes` and `training-spec` to the fork? Their text discusses PSI-derived results.

## Housekeeping
- [ ] Delete `benchmarks/runs/case_g2/{probe_opal,probe_opal_gev,try_opal}` and the older `step1` (left for A. Ghribi).
