# Overview

LIONS tracks ions in a separated-sector cyclotron, from injection in the central region to extraction.
It integrates the relativistic Lorentz equation with a Boris (leap-frog) pusher in static magnetic maps and RF electric potential maps.
The Python package `lions` is a port of the Fortran 90 code by P. Bertrand (GANIL).

## History

- Fortran 90 version: P. Bertrand, GANIL, December 1998 to December 2004.
- Version 2: March 2017 (random seed, byte swap of binary files).
- Python port and PyQt5 GUI: O. Kamalou, GANIL, 2026.

## Physics content

- Magnetic field: median-plane map in (r, theta) with Bz, d2Bz/dr2, dBz/dx, dBz/dy. Bilinear interpolation, first-order expansion off the median plane. Optional first and second field harmonics added to the map.
- Electric field: two potential maps.
  - Central region (r < `central.rmax`): potential on a Cartesian (x, y) grid, binary Fortran file, with d2V/dz2 for Ez.
  - Outer region: normalised single-gap potential in a frame rotated by the dee half-angle, text file.
  - Both are scaled by `V_dee cos(omega_rf t + phase)`, with separate voltage and phase for the right (x >= 0) and left (x < 0) dees.
- Extraction: electrostatic deflectors, magnetic gradient dipoles, extraction grids that record crossing particles.
- Initial bunch: particles distributed around a reference particle from a 6x6 sigma matrix (Cholesky factorisation).
- Vertical sorter: an oscillating Ez in an (r, theta) sector.
- Ion traps: circles read from a marker file; particles inside are lost.
- Not modelled: space charge, full 3D magnetic map, stripping and charge exchange, more than two dees.

## Use cases

1. Check that a set of machine parameters accelerates and extracts a given ion.
2. Tune RF frequency, dee voltages and phases, field level (`b_factor`, dB/B) for transmission and beam quality.
3. Size electrostatic deflectors, magnetic dipoles and extraction grids.
4. Follow phase-space correlations (r-pr, z-pz, phi-dW/W) during acceleration.
5. Study the effect of an oscillating vertical field (vertical sorter).

## Files needed for a run

| File | Deck field | Content |
|---|---|---|
| Input deck | (passed to `load`) | Physical and numerical parameters, see [Inputs](02_inputs.md) |
| Magnetic field map | `ficb` | Bz, d2Bz/dr2, dBz/dx, dBz/dy on (r, theta), text |
| Central potential | `ficv1` | Binary Fortran file, or `no` |
| Outer potential | `ficv2` | Text file, normalised gap potential |
| Sigma file | `ficsig` | 6x6 sigma matrices, only with `iflcor != 0` |
| Marker file | (prompt) | Lines and circles for the bidims, traps. Read only when bidims are on |
| Bunch file | `ficloa` | Saved bunch to start from, optional |

The CIME decks in `examples/` name maps and a sigma file that are not in the repository.

## Package layout

| Module | Content | Main names |
|---|---|---|
| `lions.state` | All data of a run | `State` and its component classes, array size constants |
| `lions.deck` | Input deck as a dataclass, reader and writer | `Deck`, `Family`, `Deflector`, `Dipole`, `Grid`, `FieldHarmonic`, `read_deck`, `write_deck` |
| `lions.io` | State from a deck, field map and potential readers | `load`, `read_bfield_map`, `read_central_potential`, `read_outer_potential`, `prompt` |
| `lions.fortran_text` | Fortran list-directed text reading | `FortranFile`, `fortran_float`, `parse_floats`, `parse_ints` |
| `lions.fields` | B and E at the particles, field harmonics | `compute_fields`, `compute_bfield`, `add_harmonics`, `add_vertical_sorter_field`, `add_extraction_fields` |
| `lions.pusher` | Boris pusher | `boris_push` |
| `lions.apertures` | Region and element masks, traps | `check_box`, `check_central`, `check_outer`, `check_deflector`, `check_dipole`, `check_grid`, `apply_traps` |
| `lions.beam` | Initial 6D distribution, phase-space coordinates | `distribute_bunch`, `cholesky_upper`, `sample_unit_ball`, `phase_space_coords`, `project_on_radial_line`, `projected_phase_space` |
| `lions.diagnostics` | z, xy and correlation outputs, end-of-run block, bunch save | `write_outputs` |
| `lions.bidim` | 2D histograms in the GANIL binary format, marker file, traps setup | `BidimWriter`, `setup_bidims`, `fill_bidims` |
| `lions.cli` | Tracking loop, scripted run, interactive program | `track`, `run`, `main`, `set_vertical_sorter`, `read_vertical_sorter`, `write_trapped` |
| `lions.constants` | Physical constants | `C_LIGHT`, `E_CHARGE`, `MEV`, `AMU`, `ELECTRON_MASS_U`, `SQRT8` |
| `lions.gui` | PyQt5 GUI | `app.LionsMainWindow`, `input_form.InputFormPanel`, `worker.SimulationWorker`, `plots`, `field_viewer`, `bidim_viewer` |

`lions` exports `State`, `load`, `run`, `track` and `set_vertical_sorter`.
State attribute names and their Fortran counterparts are listed in the [glossary](../glossary.md).
Signatures are in the [API reference](../api.md).

### Array sizes

Fixed sizes inherited from the Fortran module, defined in `lions.state`:

| Constant | Value | Limit on |
|---|---|---|
| `CENTRAL_NX`, `CENTRAL_NY` | 500, 300 | Central potential grid |
| `OUTER_NX`, `OUTER_NY` | 50, 1000 | Outer potential grid |
| `MAP_NR`, `MAP_NTHETA` | 250, 721 | Magnetic map grid (r, theta) |
| `MAX_BIDIMS` | 10 | Bidim histogram scalings |
| `MAX_BARY` | 10800 | Orbit-centre barycentre buffer |

Each field harmonic takes at most 10 coefficients.
Families, extraction elements and trap circles are lists, with no fixed size.
`load` creates the per-particle arrays as `Bunch(n_particles)`.

## Call sequence

`cli.main` (command `lions`) and the GUI worker run the same sequence.
`cli.run` skips the prompts, the vertical sorter and the bidims.

```
main()
  banner()                          prints the banner and the byte-swap prompt
  State()
  run.swap_bytes  <- stdin          y / n
  load(d, deck)                     deck, field map, harmonics, b_factor, dB/B,
                                    potentials, bunch distribution, output files
  read_vertical_sorter(d)           two lines from stdin
  setup_bidims(d, BidimWriter())    first bidim number, bidim prompts, marker file, traps
  track(d, bidim)
  write_trapped(d.bunch, d.traps)   optional list of trapped particles

track(d, bidim, on_step)
  step 0:   check_box -> compute_fields -> write_outputs -> on_step
  for step = 1 .. n_steps:
            apply_traps             losses (only when traps are defined)
            check_box               r, in_box
            compute_fields          B, E central, E outer, vertical sorter,
                                    deflectors, dipoles, grids
            boris_push              momenta, gamma, positions
            t = step * dt + t0
            write_outputs           turn count, z / xy / correlation files
            fill_bidims             if bidims on and (step % bidims.every == 0 or step == 1)
            on_step(d)              stops the loop if it returns False
  finally:  close bidims and output files
```

The fields used in step `n` are evaluated at `t = (n - 1) dt + t0`; `t` is advanced after the push.
