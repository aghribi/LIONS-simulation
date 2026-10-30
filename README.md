# LIONS

Ion tracking in separated-sector cyclotrons, from injection in the central region to extraction.
Python port of the Fortran 90 code by P. Bertrand (GANIL, 1998-2004, v2 March 2017), with a PyQt5 GUI.

<!-- docs-start -->

## Physics

- Relativistic Lorentz equation, Boris leap-frog pusher, fixed time step `dt = harm / (fhf * npaspt)`.
- Magnetic field: median-plane map in (r, theta): Bz, d2Bz/dr2, dBz/dx, dBz/dy. Bilinear interpolation, first-order expansion off the median plane (Bx = z dBz/dx, By = z dBz/dy, Bz - z²/2 d2Bz/dr2). Full 360° map or 1/8 map extended by symmetry. Optional first and second harmonic added to the map.
- Electric field: 3D potential map in the central region (Cartesian, binary Fortran file), normalised single-gap potential in the outer region (text file, frame rotated by the dee half-angle). Both scaled by `V_dee cos(omega_HF t + dphi)`, left and right dees set separately.
- Extraction: electrostatic deflectors, magnetic dipoles, extraction grids.
- Initial bunch: 6D sigma matrix (r, pr, z, pz, phi, dp/p), sampled through Cholesky.
- Vertical sorter (oscillating Ez), ion traps (pillars, lines) read from a marker file.
- Not modelled: space charge, full 3D magnetic map, stripping and charge exchange, more than 2 dees.

## Install

```bash
uv sync                 # numpy, PyQt5, matplotlib, pyqtgraph
uv sync --extra gl      # + PyOpenGL for the 3D trajectory view
```

## Run

GUI:

```bash
uv run python lions_gui.py
```

Command line. The engine reads its answers from stdin, in this order: byte swap for the central potential file (`y`/`n`), input deck, vertical sorter (two lines), first bidim number (`0` for none).

```bash
printf 'n\ndon.120.ref_phi10_mars_2017\n1000 2000 0 0\n0 0 0 0 0 0\n0\n' | uv run python lions.py
```

The input deck format is described line by line in section 4 of `DOCUMENTATION_LIONS.md`.

## Files

| File | Role |
|---|---|
| `lions.py` | Engine: deck and map readers, field interpolation (`calcbe`), Boris pusher (`pousse`), losses (`piege`), outputs (`info`), bidims, main loop (`lions`) |
| `lion_m.py` | `LionData`: all state (scalars, per-particle arrays, field maps), port of the Fortran module |
| `lions_worker.py` | QThread that runs the engine for the GUI, feeds `input()` from the form, captures `print()` |
| `lions_gui.py`, `lions_input_form.py` | Main window and parameter form (reads and writes decks) |
| `lions_plots.py`, `lions_field_viewer.py`, `lions_bidim_viewer.py` | 2D/3D plots, field map viewer, bidim viewer |
| `lions_2.0_ep.f90`, `lion_m.f90` | Fortran reference (`gfortran -o lions_2.0_ep lion_m.f90 lions_2.0_ep.f90`) |
| `lions_2.0_ep_windows.f90` | Same code, reads the central potential and writes saved bunches as big-endian |
| `don.120.ref_phi10_mars_2017` | Example deck (CIME). Its maps and sigma file are not in the repository |
| `generate_pdf.py` | PDF of `DOCUMENTATION_LIONS.md` (needs `fpdf2`) |

## Units

Internal state is SI: positions in m, time in s, B in T, E in V/m. Momenta are normalised, `p = beta gamma`.
Deck inputs mix units: r0 and z0 in m, theta0 and phases in deg, pr0 and pz0 in mrad, energy in MeV/u, extraction grids in mm.
Text outputs are in mm, mrad, kV/cm and MeV/u.

## Known issues

From a full read of `lions.py` against the Fortran source (2026-10-30). None fixed yet.

1. **B-map scaling differs from the Fortran.** The Python engine divides the map by max|Bz| and multiplies by `b_factor` (T), read from an optional line after the three map files. If the line is missing, `b_factor = 1.0` and the field peaks at 1 T. The Fortran uses the map as read. The GUI hides this: it auto-detects Bmax when the deck has no `b_factor`. The example deck has no such line, so a command-line run of it changes the physics.
2. **Random bunch differs from the Fortran.** Python seeds the legacy global RNG with 42 inside `remp`. The Fortran seeds the gfortran generator (its seed loop sets all 12 entries to 12). Bunches cannot match particle by particle. Python vs Fortran checks have to be statistical, or start from a single particle or a loaded bunch.
3. **Saved bunch files are not Fortran-compatible.** The Fortran writes `ficsto` with record markers and reads `ficloa` with them, honouring the byte swap. Python writes and reads raw native doubles and ignores the swap.
4. **Losses are counted only when traps are defined.** As in the Fortran, `piege` returns at once without a marker file, so particles beyond zmin/zmax are never flagged lost. The marker file is read only when bidims are on.
5. **Interactive I/O in the engine.** `init`, `initgra` and `lions` call `input()`. `lions_worker.run` repeats the main loop of `lions`. Output files sit in a module-level dict (`files_dict`). There is no function that takes parameters and returns results.
6. **Speed.** Per-particle Python loops in `calcbe`, `calcb` and `bid`.
7. **Duplicated readers.** `lions_field_viewer.py` has its own field-map readers, separate from `init`.
8. **No tests**, and no Python vs Fortran comparison on a reference case.

<!-- docs-end -->

## Documentation

Full technical documentation (French): `DOCUMENTATION_LIONS.md`. Sphinx site with the same content and the API:

```bash
uv sync --group docs
uv run sphinx-build -b html docs docs/_build/html
```

## Credits

- Fortran code: P. Bertrand (GANIL).
- Python port and GUI: O. Kamalou (GANIL). Original repository: https://github.com/omarino1979/LIONS-simulation.

## License

MIT, see `LICENSE`.
