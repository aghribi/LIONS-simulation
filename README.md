# LIONS

Ion tracking in separated-sector cyclotrons, from injection in the central region to extraction.
Python port of the Fortran 90 code by P. Bertrand (GANIL, 1998-2004, v2 March 2017), with a PyQt5 GUI.

<!-- docs-start -->

## Physics

- Relativistic Lorentz equation, Boris leap-frog pusher, fixed time step `dt = harmonic / (rf_frequency * steps_per_turn)`.
- Magnetic field: median-plane map in (r, theta): Bz, d2Bz/dr2, dBz/dx, dBz/dy. Bilinear interpolation, first-order expansion off the median plane (Bx = z dBz/dx, By = z dBz/dy, Bz - z²/2 d2Bz/dr2). Full 360° map or 1/8 map extended by symmetry. Optional first and second harmonic added to the map.
- Electric field: 3D potential map in the central region (Cartesian, binary Fortran file), normalised single-gap potential in the outer region (text file, frame rotated by the dee half-angle). Both scaled by `V_dee cos(omega_rf t + phase)`, left and right dees set separately.
- Extraction: electrostatic deflectors, magnetic dipoles, extraction grids.
- Initial bunch: 6D sigma matrix (r, pr, z, pz, phi, dp/p), sampled through Cholesky.
- Vertical sorter (oscillating Ez), ion traps (pillars, lines) read from a marker file.
- Not modelled: space charge, full 3D magnetic map, stripping and charge exchange, more than 2 dees.

## Install

```bash
uv sync                 # numpy, PyQt5, matplotlib, pyqtgraph
uv sync --extra gl      # + PyOpenGL, PyQtWebEngine, plotly for the 3D trajectory view
```

## Run

GUI:

```bash
uv run lions-gui
```

Command line. The program reads its answers from stdin, in this order: byte swap for the central potential file (`y`/`n`), input deck, vertical sorter (two lines), first bidim number (`0` for none).

```bash
cd examples
printf 'n\ndon.120.ref_phi10_mars_2017\n1000 2000 0 0\n0 0 0 0 0 0\n0\n' | uv run lions
```

Python, no prompt, no bidims. `run` returns the final `State`; output files named in the deck are written as usual.
The state is grouped in components (`run`, `rf`, `field`, `bunch`, `extraction`, ...), see `docs/glossary.md`.

```python
import lions

d = lions.run('don.120.ref_phi10_mars_2017', swap=False)
b = d.bunch
b.x, b.y, b.z, b.px, b.py, b.pz, b.alive   # per particle, SI, p = beta gamma
d.rf.v_right, d.field.bz, d.extraction.deflectors   # other components

# on_step is called after every step; returning False stops the run
d = lions.run('don.120.ref_phi10_mars_2017', on_step=lambda d: d.run.step < 1000)
```

Decks as data: `read_deck` returns a `Deck` dataclass (file units, English names), `write_deck` writes one that the Fortran also reads. `run` and `load` take a `Deck` or a path.

```python
deck = lions.read_deck('don.120.ref_phi10_mars_2017')
deck.v_dee_right = deck.v_dee_left = 52e3   # V
deck.families[0].phi0_deg = -245.0
d = lions.run(deck)
lions.write_deck(deck, 'don.modified')
```

`vertical` takes the keyword arguments of `lions.set_vertical_sorter` (lengths in m, angles in deg); the default is no vertical sorter.
State names and their Fortran counterparts: `docs/glossary.md`.

The input deck format is described line by line in `docs/guide/02_inputs.md`.
The deck ends with two optional entries, `[b_factor] [ficloa ficsto]`.
Without `b_factor`, the magnetic map is used as read, as in the Fortran.
With it, the map is scaled so that max|Bz| = `b_factor` (T).
The GUI always writes `b_factor`.

## Layout

| Path | Content |
|---|---|
| `src/lions/state.py` | `State` and its components (`Run`, `RF`, `FieldMap`, `Region`, `Bunch`, element lists, ...) |
| `src/lions/deck.py` | `Deck` dataclass, `read_deck`, `write_deck` (engine and GUI) |
| `src/lions/io.py` | `load`: State from a deck; field map and potential readers |
| `src/lions/fortran_text.py` | Fortran list-directed text reading |
| `src/lions/fields.py` | B and E at the particles, field harmonics |
| `src/lions/pusher.py` | Boris pusher |
| `src/lions/apertures.py` | Box, central and outer regions, deflectors, dipoles, grids, traps |
| `src/lions/beam.py` | Initial 6D distribution, phase-space coordinates |
| `src/lions/diagnostics.py` | z, xy and correlation outputs |
| `src/lions/bidim.py` | 2D histograms (bidims) |
| `src/lions/cli.py` | Tracking loop (`track`), scripted run (`run`), interactive program (`main`) |
| `src/lions/gui/` | Main window, parameter form, worker thread, plots, field and bidim viewers |
| `fortran/` | Fortran reference. `lions_2.0_ep_windows.f90` reads the central potential and writes saved bunches as big-endian |
| `examples/` | CIME decks. Their maps and sigma file are not in the repository |
| `docs/` | Sphinx sources, guide pages, `generate_pdf.py` (PDF of the guide) |
| `scripts/lions-gui.bat` | Windows launcher for `uv run lions-gui` |

## Units

Internal state is SI: positions in m, time in s, B in T, E in V/m. Momenta are normalised, `p = beta gamma`.
Deck inputs mix units: r0 and z0 in m, theta0 and phases in deg, pr0 and pz0 in mrad, energy in MeV/u, extraction grids in mm.
Text outputs are in mm, mrad, kV/cm and MeV/u.

## Known issues

From a read of the Python port against the Fortran source (2026-10-30).

1. **Random bunch differs from the Fortran.** Python seeds the legacy global RNG with 42 inside `sample_unit_ball`. The Fortran seeds the gfortran generator (its seed loop sets all 12 entries to 12). Bunches cannot match particle by particle. Python vs Fortran checks have to be statistical, or start from a single particle or a loaded bunch.
2. **Saved bunch files are not Fortran-compatible.** The Fortran writes `ficsto` with record markers and reads `ficloa` with them, honouring the byte swap. Python writes and reads raw native doubles and ignores the swap.
3. **Losses are counted only when traps are defined.** As in the Fortran, `apply_traps` returns at once without a marker file, so particles beyond zmin/zmax are never flagged lost. The marker file is read only when bidims are on.
4. **`ficv1 = 'no'` still prompts** for the central grid (`xmin1 xmax1 ymin1 ymax1 nx1 ny1`), so `run` waits on stdin in that case.
5. **Speed.** Per-particle Python loops in `compute_fields`, `compute_bfield` and `fill_bidims`.
6. **No bounds check on the central grid.** A particle with r < rmax1 outside the pot1 grid gets a negative index: NumPy wraps it to the far end of the array, the Fortran reads outside it. Both give meaningless fields. The pot1 grid must cover r < rmax1.
7. **Correlation outputs only partly ported.** With `iflcor != 0`, the Fortran also writes the sigma and grid matrices, the window-filtered matrices (two extra prompts at the end) and the inflector output to `ficmat`, and the dc matrices to `ficdc`. The Python opens these files and leaves them empty.
8. **`6d` bunch differs from the Fortran.** The Fortran adds the vertical offsets (dz, dpz) of the sampled bunch only when `espace = '3d'` (`lions_2.0_ep.f90` line 498). The Python adds them for `3d` and `6d`. The CIME example deck uses `6d`: the Fortran keeps all particles at z0, the Python spreads them vertically. Which behaviour is intended is open.
9. **Correlations wrong for families after the first, in both codes.** The orbit radius uses `mass_c[k]` (Fortran `am0c(k)`), the particle array indexed by the family number instead of the family's first particle.
10. **Untested paths.** Multiple families, the vertical sorter, harmonics, `ficmat`/`ficdc`. No comparison on CIME maps yet.

## Tests

```bash
uv run pytest
```

Inputs are synthetic, written per test into a temporary directory (`tests/synthetic.py`).

- `tests/test_run.py`: 1 MeV proton in a uniform Bz, zero E. Orbit radius p/(qB) within 1e-4, energy constant within 1e-12, Boris phase lag 2 atan(w dt/2) per step within 1 %, `b_factor` scaling, deck tail parsing, `on_step` stop.
- `tests/test_fortran.py`: the same inputs through `fortran/lions_2.0_ep.f90` (gfortran, `-ffp-contract=off`) and `python -m lions`, both run as programs with the same stdin. Six cases: outer-region dees with a 4-sector field and off-plane terms; the same field as a 1/8 map; central-region potential with Ez; deflector, magnetic dipole and two extraction grids; 20-particle loaded bunch with correlations; the same bunch with bidims and a trap. Final bunches agree within 1e-12 (bit-identical on macOS arm64); text outputs (z, xy, correlations, grid crossings, trapped list) agree within their printed precision. Skipped without gfortran.
- `tests/test_gui.py`: headless main window, deck round trip through the input form, worker run, field viewer on both map formats.

<!-- docs-end -->

## Documentation

Technical guide: `docs/guide/` (inputs, outputs, coordinates, algorithms, bunch generation, GUI, worked example).
Sphinx site with the guide, the glossary and the API: https://aghribi.github.io/LIONS-simulation/ (built from `main`). Local build:

```bash
uv sync
uv run sphinx-build -b html docs docs/_build/html
```

PDF of the guide, written to `docs/_build/pdf/LIONS.pdf`:

```bash
uv run --with fpdf2 python docs/generate_pdf.py
```

## Credits

- Fortran code: P. Bertrand (GANIL).
- Python port and GUI: O. Kamalou (GANIL). Original repository: https://github.com/omarino1979/LIONS-simulation.

## License

MIT, see `LICENSE`.
