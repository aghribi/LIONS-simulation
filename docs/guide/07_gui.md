# GUI

```bash
uv run lions-gui
```

On Windows, `scripts/lions-gui.bat` runs the same command from the repository folder.

## Modules

| Module | Content |
|---|---|
| `lions.gui.app` | `LionsMainWindow`: toolbar, parameter dock, plot area, console dock, status bar; `main` |
| `lions.gui.input_form` | `InputFormPanel`: parameter form in 16 collapsible sections (`CollapsibleSection`); deck read and write |
| `lions.gui.worker` | `SimulationWorker` (`QThread`): runs the engine with preloaded stdin answers |
| `lions.gui.plots` | `PlotWidget`, `PlotCanvas`: six live plots (pyqtgraph); `Plot3DWindow`: 3D view |
| `lions.gui.field_viewer` | `FieldViewerDialog`: maps of Bz, d2Bz, dBz/dx, dBz/dy and of the two potentials |
| `lions.gui.bidim_viewer` | `BidimViewerDialog`, `BidimViewerWidget`, `read_bidim`: `.bid` viewer (matplotlib) |

## Toolbar

| Action | Effect |
|---|---|
| Run | Writes the form to a deck and starts the worker |
| Stop | Sets `stop_flag`; the run stops after the current step |
| Reset | Clears plots, 3D view and console |
| Unzoom | Resets the zoom of all plots |
| Export | Saves the plots as PNG, PDF or SVG |
| Bidims | Opens the bidim viewer on the `*.bid` files next to the deck |
| Fields | Opens the field viewer; files can be chosen in the dialog |
| 3D | Opens the 3D trajectory view |

The 3D action is shown only when `pyqtgraph.opengl` imports (`uv sync --extra gl` installs PyOpenGL).
The 3D window itself draws with Plotly in a Qt WebEngine view and needs the `plotly` and `PyQtWebEngine` packages, which the project does not declare; without them it shows the install hint.

## Parameter form

Sections, in order:

| # | Section | Content (deck fields) |
|---|---|---|
| 1 | General | byte swap, deck path, Load data file button |
| 2 | Description | `corps` |
| 3 | RF / Cyclotron | `npaspt`, `fhf`, `pothf`, `teta1`, `harm`, `dbsb` |
| 4 | Radius / Turns | `r1th`, `ntours` |
| 5 | Dees | `dphid`, `dphig`, `pothfd`, `pothfg` |
| 6 | Harmonics | `ncoh`, coefficients and phases |
| 7 | Z geometry | `zmin`, `zmax`, `rmax1` |
| 8 | Extraction | `rex`, `tet1ex`, `tet2ex`, `eex`, deflectors, dipoles, grids |
| 9 | Ion(s) | families: `nbrion`, `amat`, `a`, `qn`, initial conditions; total `n` |
| 10 | Correlations | `iflcor`, `ficsig`, `espace`, diagonal sigma, `ficcor ficmat ficdc` |
| 11 | Simulation | `nbdt`, `ifreq` |
| 12 | Output files | `ficz`, `ficxy` |
| 13 | Field files | `ficb`, `ficv1`, `ficv2`, `b_factor` |
| 14 | Vertical sorter | stdin answers 4 and 5 |
| 15 | Bidim / Diagnostics | stdin answers 6 to 15 |
| 16 | Post-Simulation | stdin answers 18 and 19 |

Answers 4 to 19 refer to the stdin table in [Inputs](02_inputs.md#standard-input).
Choosing a magnetic map sets `b_factor` to the max|Bz| of the file.
The diagonal sigma fields are described in [Bunch generation](06_bunch.md#gui).

Buttons under the form:

- Save params, Load params: all form values as JSON.
- Save data file: writes the form as a deck (`InputFormPanel.write_data_file`).

`InputFormPanel.load_data_file` reads a deck with the same line structure as `lions.io.load` and emits `file_loaded`.
The form has no fields for `yomega` and `ficloa ficsto`: it keeps the values of the deck last loaded and writes them back.

## Run sequence

1. `on_run` checks that the deck, `ficb` and `ficv2` are set, otherwise it shows a warning and stops.
2. `get_input_sequence` writes the sigma file (when `iflcor != 0`) and the form to `_lions_gui_temp.dat` next to the deck, and builds the stdin answers in the order of `lions.cli.main`. `get_post_answers` adds the trapped-list answers.
3. `SimulationWorker` changes to the deck folder (paths in the deck are relative to it), replaces `builtins.input` with the answer queue (`InputQueue`) and `sys.stdout` with `OutputInterceptor`, then runs `banner`, `load`, `read_vertical_sorter`, `setup_bidims`, `track` and `write_trapped`. Input, stdout and working folder are restored at the end.
4. The `on_step` hook of `track` emits `progress` about 200 times per run, emits plot data at step 0, at step 1, every `ifreq` steps and at the last step, and returns `False` when Stop was pressed.
5. A 100 ms timer flushes the buffered plot data and console text.

Worker signals:

| Signal | Arguments | Meaning |
|---|---|---|
| `progress` | int, int | current step, `nbdt` |
| `iteration_data` | dict | plot data, one dict per family |
| `console_output` | str | printed text |
| `simulation_finished` | | normal end |
| `simulation_error` | str | exception message and traceback |
| `input_needed` | str | the engine asked more answers than queued |

## Plots

Six plots, three rows of two, one curve per family plus a scatter of the live particles when a family has more than one:

| Plot | x | y |
|---|---|---|
| XY trajectory | X (mm) | Y (mm), deflectors drawn |
| Energy vs turn | turn | W (MeV/u) |
| R - theta | theta (deg) | R (mm) |
| Z - Pz phase space | Z (mm) | Pz (mrad) |
| R - dphi | R (mm) | dphi (deg) |
| dphi vs turn | turn | dphi (deg) |

Curves follow the first particle of each family.
`dphi = omega_rf t - harm * theta_total - 90` (deg) and Pz is `asin(pz / p) * 1000` (mrad).
Controls under the plots limit the radius range and filter dphi.

## Viewers

**Field viewer.** `LionsMainWindow.on_view_fields` pre-fills the map names by skipping a fixed number of deck lines; the count does not follow the deck format, so the names are wrong for most decks (the CIME example included) and the files have to be chosen by hand. The dialog plots Bz, d2Bz, dBz/dx, dBz/dy and the two potentials and their d2pot, as 2D maps or surfaces, with a choice of colour map. It reads the files with `lions.io.read_bfield_map`, `read_central_potential` (falls back to the other byte order) and `read_outer_potential`.

**Bidim viewer.** Reads `.bid` files with `read_bidim` and shows them as 2D histograms on a log scale. Colour maps: hot, viridis, plasma, inferno, jet, gray. Bins: 128, 256, 512, 1024.

## Errors in the GUI

| Situation | Behaviour |
|---|---|
| No deck, or no `ficb` / `ficv2` | Warning dialog, run not started |
| Run already active | Warning dialog |
| Exception or `sys.exit` in the engine | `simulation_error`, error dialog |
| Engine asks an answer not in the queue | `input_needed`, input dialog; Cancel sends `0`; no answer within 300 s raises an error |
| Deck load fails | Warning dialog with the message |
| Plot data error | Ignored, the run continues |
| Window closed during a run | Stop requested, waits up to 5 s |
