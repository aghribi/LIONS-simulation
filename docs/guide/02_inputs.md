# Inputs

## Input deck

The input deck is a text file. `lions.deck.read_deck` reads it into a `Deck` dataclass (file units, English names), `lions.deck.write_deck` writes one, and `lions.io.load` sets up the `State` from a `Deck` or a path.
The GUI form reads and writes decks through the same two functions.
`write_deck` quotes file names (Fortran list-directed input ends a value at an unquoted `/`) and writes numbers with all their digits, so a read, write, read cycle is exact.
Field names are those of the Fortran code: they define the file format.
Each field fills the `State` attribute given in the tables below.

### Reading rules

- Blank lines are skipped.
- Most numeric lines are read as one line split on blanks. Extra values on a line are ignored.
- File names on lines 16 to 22 and the three correlation outputs are read as Fortran list-directed tokens: they can share a line or span several lines, and may be quoted with `'` or `"`.
- Numbers accept the Fortran exponent: `1.50d0`, `3.14D-4` (`lions.io.fortran_float`).
- Angles in the deck are in degrees. `load` converts them to rad, except where the table says otherwise.

### Line by line

| Line | Fields | Units | State attribute | Condition | Example |
|---|---|---|---|---|---|
| 1 | `corps` | text | `run.title` | | `"don.forw.H3 i700 posbmax "` |
| 2 | `npaspt` | - | `run.steps_per_turn` | | `1080` |
| 3 | `r1th ntours` | m, - | `run.design_r_final`, `run.design_turns` | | `1.50d0 250` |
| 4 | `fhf pothf teta1 harm dbsb` | Hz, V, deg, -, - | `rf.frequency`, (unused), `rf.dee_half_angle`, `rf.harmonic`, `field.db_over_b` | | `14400000.d0 -51890.d0 20.d0 3.d0 +1.0d-4` |
| 5 | `dphid dphig` | deg | `rf.phase_right`, `rf.phase_left` | | `0 0` |
| 6 | `pothfd pothfg` | V | `rf.v_right`, `rf.v_left` (sign reversed) | | `51890.d0 51890.d0` |
| 7 | `ncoh(1)` | - | `field.harmonic_n_coef[0]` | | `0` |
| 7b | `coh(1..ncoh(1), 1)` | see below | `field.harmonic_coef[:, 0]` | `ncoh(1) != 0` | |
| 7c | `phihar(1)` | deg | `field.harmonic_phase[0]` (kept in deg) | `ncoh(1) != 0` | |
| 8 | `ncoh(2)` | - | `field.harmonic_n_coef[1]` | | `0` |
| 8b | `coh(1..ncoh(2), 2)` | see below | `field.harmonic_coef[:, 1]` | `ncoh(2) != 0` | |
| 8c | `phihar(2)` | deg | `field.harmonic_phase[1]` (kept in deg) | `ncoh(2) != 0` | |
| 8d | `yomega` | m | `rf.gap_offset` | `teta1 < 0` | `0.01d0` |
| 9 | `zmin zmax rmax1` | m | `box.zmin`, `box.zmax`, `central.rmax` | | `-.5d0 .5d0 .164d0` |
| 10 | `rex tet1ex tet2ex eex` | -, -, -, - | `extraction.mode`, `extraction.cut_bunch`, `extraction.grids_on`, `extraction.grid_start` | | `-1. -1. -1. -1.` |
| 10b | `ndf` | - | `len(extraction.deflectors)` | `rex < 0` | `1` |
| 10c | `rdf edf xdfi ydfi xdfo ydfo tetdfi tetdfo` | m, V/m, m, m, m, m, deg, deg | `deflector_*` | `ndf` lines | `1.513 6990000 0 0 0 0 250 267` |
| 10d | `ndm` | - | `len(extraction.dipoles)` | `rex < 0` | `0` |
| 10e | `rdm bdm xdm ydm tetdmi tetdmo` | m, T/m, m, m, deg, deg | `dipole_*` | `ndm` lines | |
| 10f | `ngex` | - | `len(extraction.grids)` | `rex < 0` and `tet2ex > 0` | `2` |
| 10g | `tetagex xgex ygex xgexc ygexc` | deg, mm, mm, mm, mm | `extraction.grids[k].angle_deg`, `extraction.grids[k].x`, `extraction.grids[k].y`, `extraction.grids[k].xc`, `extraction.grids[k].yc` (m) | `ngex` lines | |
| 10h | `ficdebgex` | text | `extraction.grid_file_prefix` | `rex < 0` and `tet2ex > 0` | `grid_ext` |
| 11 | `nbrfam` | - | `len(bunch.families)` | | `1` |
| 11b | `n` | - | `bunch.n_particles` | `nbrfam > 1` | |
| 12 | `nbrion amat a qn` | -, u, -, e | `bunch.families[k].size`, `bunch.families[k].atomic_mass`, `bunch.families[k].mass_number`, `bunch.families[k].charge_state` | once per family, before line 13 | `100 100.36496d0 100.36496d0 22.d0` |
| 13 | `r0 teta0 z0 pr0 pz0 w0mev phi0` | m, deg, m, mrad, mrad, MeV/u, deg | `bunch.families[k].r0`, `bunch.families[k].theta0`, `bunch.families[k].z0`, `bunch.families[k].pr0`, `bunch.families[k].pz0`, `bunch.families[k].w0`, `bunch.families[k].phi0_deg` | once per family | `0.0331618 266.8973517 0.d0 426.2811142 0.d0 5.6650d-3 -250.00d0` |
| 14 | `iflcor` | - | `bunch.correlations` | | `1` |
| 14b | `ficsig` | text | `files.sigma` | `iflcor != 0` | `sigma6d_no_corr.dat` |
| 14c | `espace` | text | `bunch.phase_space` | `iflcor != 0` | `6d` |
| 14d | `ficcor ficmat ficdc` | text | `files.correlation`, `files.matrix`, `files.dc_matrix` | `iflcor != 0` | `no no no` |
| 15 | `nbdt ifreq` | - | `run.n_steps`, `run.output_every` | | `350000 30` |
| 16 | `ficz` | text | `files.z_output` | | `exp.for.z.ref_phi10_mars_2017` |
| 17 | `ficxy` | text | `files.xy_output` | | `exp.for.ref_phi10_mars_2017` |
| 18 | `ficb` | text | `files.bfield` | | `'cime700rec.1804.pbob1.lio.dat'` |
| 19 | `ficv1` | text | `files.central_potential` | | `'pot3d.impair.30pilbleu30.620.dat'` |
| 20 | `ficv2` | text | `files.outer_potential` | | `'newgap612.dat'` |
| 21 | `b_factor` | T | `field.b_factor` | optional | `1.56d0` |
| 22 | `ficloa ficsto` | text | `files.bunch_load`, `files.bunch_save` | optional, default `no no` | `no bunch.bin` |

Lines 12 and 13 repeat for each family, in the order 12, 13, 12, 13, ...
Lines 14b and 14c are read as whole lines: quotes are kept as part of the name.
Line 10h is read the same way.

### Run parameters

- `npaspt`: time steps per RF period times `harm`, that is per revolution of an isochronous orbit. The time step is `dt = harm / (fhf * npaspt)` (s). A negative `npaspt` gives a negative `dt` and tracks backwards; the turn counter then starts at `-nbdt / npaspt` and counts down.
- `r1th`, `ntours`: design final radius (m) and number of turns. Read and stored, not used by the tracking.
- `nbdt`: number of time steps. Tracked time `nbdt * dt` (s).
- `ifreq`: output period in steps. Outputs are written at steps that are multiples of `ifreq`, at step 1 and at the last step (see [Outputs](03_outputs.md)).

Example: `fhf` = 14.4 MHz, `harm` = 3, `npaspt` = 1080 give `dt` = 0.1929 ns and a revolution period of 208.3 ns. `nbdt` = 350000 covers 67.5 us, 324 turns. `ifreq` = 30 writes 36 lines per turn.

### RF and dees

- `fhf`: RF frequency (Hz). `omega_rf = 2 pi fhf` (rad/s), `omega_rev = omega_rf / harm`.
- `pothf`: read and not used. The dee voltages come from line 6.
- `teta1`: dee half-angle (deg). The outer-region potential frame is rotated by this angle. A negative value marks a skewed gap: `load` makes it positive and reads `yomega` (line 8d).
- `harm`: harmonic number. `harmonic_sign` = -1 for an even harmonic, +1 for an odd one. It sets the field symmetry between quadrants in the outer region.
- `dbsb`: relative correction dB/B, applied to Bz and d2Bz/dr2 after the map is read.
- `dphid`, `dphig`: phase offsets of the right (x >= 0) and left (x < 0) dees (deg, stored in rad).
- `pothfd`, `pothfg`: right and left dee voltages (V). Stored with the sign reversed (`v_dee_right = -pothfd`).

### Field harmonics (lines 7 and 8)

`ncoh(1)` and `ncoh(2)` give the number of coefficients (at most 10) of the first and second harmonic added to the magnetic map.
Zero means none, and lines 7b, 7c (or 8b, 8c) are absent.
The coefficients are on one line; the phase (deg) is on the next.
`lions.fields.add_harmonics` adds them to Bz, d2Bz, dBz/dx and dBz/dy on every map node (see [Algorithms](05_algorithms.md)).
With x, y in m and Bz in T, `coh(k)` of the first harmonic is in T/m^(2k+1) and of the second in T/m^(2k+2), k counted from 0.

### Skewed gap (line 8d)

Present only when `teta1 < 0`.
`yomega` (m) shifts the gap: the outer-region frame uses `|y| - yomega` (see [Coordinates](04_coordinates.md)).

### Vertical box and regions (line 9)

- `zmin`, `zmax` (m): vertical limits of the tracking box. With traps defined, a particle outside [zmin, zmax] is lost.
- `rmax1` (m): central / outer region boundary. For r < `rmax1` the central potential is used, beyond it the outer potential.

### Extraction (line 10)

- `rex < 0` enables extraction elements; lines 10b to 10e follow. `rex >= 0`: no deflector, dipole or grid.
- `tet1ex <= 0`: bunch cut. A particle enters the first deflector when r > `rdf(1)` (`extraction_start_step = 0`).
  `tet1ex > 0`: whole bunch. All particles may enter the first deflector from step `extraction_start_step = int(tet1ex + 0.5)`.
- `tet2ex > 0`: extraction grids are defined, lines 10f to 10h follow.
- `eex`: with grids, `grid_start_step = int(eex + 0.5)` is the first step at which grids are tested.

`load` also stores `tet1ex` and `tet2ex` multiplied by pi/180 in `cut_bunch` and `grids_on`; only their sign is used afterwards.

**Deflectors (line 10c).**
`rdf` (m) is the entry radius of the first deflector; it is not used for the others.
`edf` (V/m) is the field magnitude inside the deflector.
The entrance line goes through (`xdfi`, `ydfi`) with direction angle `tetdfi`; the exit line through (`xdfo`, `ydfo`) with angle `tetdfo`.
`load` stores the angles in rad and the unit vectors (`deflector_in_ux`, `deflector_in_uy`) = (cos, sin) of the angle, same for the exit.

**Dipoles (line 10e).**
`bdm` (T/m) is a field gradient: inside the dipole, Bz increases by `(d - rdm) * bdm`, d being the distance (m) to (`xdm`, `ydm`).
The entrance and exit lines both go through (`xdm`, `ydm`), with angles `tetdmi` and `tetdmo`.

**Grids (lines 10g, 10h).**
Values are in mm, converted to m on reading.
(`xgex`, `ygex`) is a point of the grid line, (`xgexc`, `ygexc`) a second point.
`load` computes `grid_length` = distance between them and the unit vector (`grid_ux`, `grid_uy`) from (`xgexc`, `ygexc`) to (`xgex`, `ygex`).
`tetagex` (deg) is only copied to the grid file header.
Each grid writes the file `ficdebgex.1`, `ficdebgex.2`, ... (`grid_files`).

### Ions (lines 11 to 13)

- `nbrfam`: number of families (at most 100). With more than one family, line 11b gives the total particle count `n`, and the sum of `nbrion` must equal it. With one family, `n = nbrion`.
- `amat`: atomic mass (u). The rest mass is `m0 = AMU * (amat - 0.0005486 * qn)` (kg), the mass of the missing electrons removed.
- `a`: mass number. Energies per nucleon use it: `w0 = w0mev * a` (MeV, stored in J).
- `qn`: charge state (e).
- `r0`, `teta0`, `z0`: initial position of the family (m, deg, m). All particles of the family start there before the bunch distribution.
- `pr0`, `pz0`: initial radial and vertical angles (mrad, stored in rad). The momentum direction in the median plane is `alpha = pi/2 + teta0 - pr0`; `px = p cos(alpha) cos(pz0)`, `py = p sin(alpha) cos(pz0)`, `pz = p sin(pz0)`, with `p = beta gamma`.
- `w0mev`: kinetic energy per nucleon (MeV/u). Family 0 also sets `w0_mev_per_u`.
- `phi0`: initial RF phase (deg). The start time is `t0 = phi0 (rad) / omega_rf`, from family 0.

### Correlations (line 14)

`iflcor != 0` distributes the particles of each family around its first particle with the sigma matrices of `ficsig`.
`espace` is `6d` or `3d`.
`ficcor`, `ficmat`, `ficdc` name the correlation outputs; `no` turns one off.
`iflcor = 0`: lines 14b to 14d are absent and all particles of a family start on the same point.
See [Bunch generation](06_bunch.md).

### Output files (lines 16, 17)

`ficz` and `ficxy` name the z and xy outputs of particle 0. `no` turns one off.

### Field files (lines 18 to 20)

`ficb`: magnetic map. `ficv1`: central potential, or `no` (the program then prompts for a central grid and uses a zero potential). `ficv2`: outer potential.

### Optional tail (lines 21, 22)

The deck may end with `[b_factor] [ficloa ficsto]`.
The first token is taken as `b_factor` if it parses as a number; the next two tokens are `ficloa ficsto`.

- No `b_factor`: the magnetic map is used as read, as in the Fortran.
- `b_factor` given: the map (Bz, d2Bz, dBz/dx, dBz/dy) is scaled so that max|Bz| = `b_factor` (T).
- `ficloa`: bunch file to load after the distribution, `no` for none. Read only when `iflcor != 0`.
- `ficsto`: bunch file to write at the end of the run, `no` for none.

The GUI always writes `b_factor`.

## Magnetic field map (`ficb`)

Text file, list-directed:

```
line 1   title
line 2   rmin rmax nr ntheta                  (m, m, -, -)
block 1  Bz(r, theta)                         nr rows of ntheta values (T)
block 2  d2Bz/dr2(r, theta)                   nr rows of ntheta values (T/m^2)
block 3  dBz/dx(r, theta)                     nr rows of ntheta values (T/m)
block 4  dBz/dy(r, theta)                     nr rows of ntheta values (T/m)
```

Row i is the radius `rmin + i * dr`, `dr = (rmax - rmin) / (nr - 1)`.
`lions.io.read_bfield_map` handles two layouts:

- **Full map** (`ntheta != 91`): columns cover 0 to 360 deg, `dtheta = 360 / (ntheta - 1)` deg.
- **1/8 map** (`ntheta = 91`): columns cover 45 deg down to 0 deg in 0.5 deg steps (the column order is reversed on reading). The map is extended to 721 columns over 360 deg by the 8-fold symmetry. Blocks 3 and 4, read as (a, b), become dBz/dx = (a + b) / sqrt(2) and dBz/dy = (a - b) / sqrt(2) on 0 to 45 deg before the extension.

Limits: `1 < nr <= 250`, `1 < ntheta <= 721`.
A 1/8 map without blocks 3 and 4 is accepted with dBz/dx = dBz/dy = 0: Bx = By = 0 and the motion stays in the median plane.
A full map that ends early raises `ValueError`.

Processing order in `load`: read, add harmonics, scale to `b_factor` (if given), multiply Bz and d2Bz by `(1 + dbsb)`.
LIONS uses Bz < 0 for ions turning counter-clockwise (theta increasing).

## Central potential (`ficv1`)

Fortran unformatted sequential file: each record is framed by 4-byte length markers.

```
record 1            xmin xmax ymin ymax (float64, mm)   nx ny (int32)
records 2..nx+1     pot(i, 0..ny-1)      (float64)
records nx+2..2nx+1 d2pot(i, 0..ny-1)    (float64, optional)
```

`pot` is the potential for 1 V on the dee (dimensionless). `d2pot` is d2pot/dz2 (1/m^2).
The grid header is in mm and is stored in m.
`swap_bytes = 'y'` reads the file as big-endian, otherwise little-endian.
Limits: `1 < nx <= 500`, `1 < ny <= 300`; a header outside them raises `ValueError`, often a sign of the wrong byte order.
Without the `d2pot` records, d2pot = 0 and Ez = 0 in the central region.
The grid must cover the whole disc r < `rmax1`: there is no bounds check.

With `ficv1 = no`, `load` prompts for `xmin1 xmax1 ymin1 ymax1 nx1 ny1` (m, m, m, m, -, -) and uses a zero potential.

## Outer potential (`ficv2`)

Text file, list-directed:

```
xmin xmax ymin ymax nx ny          (m, m, m, m, -, -)
pot(0, 0..ny-1)                    nx rows
...
d2pot(0, 0..ny-1)                  nx rows
...
```

`pot` is the single-gap potential normalised to 1 V (0 to 1). `d2pot` is d2pot/dz2 (1/m^2).
The grid is in the gap frame (x', y') defined in [Coordinates](04_coordinates.md).
Limits: `1 < nx <= 50`, `1 < ny <= 1000`.

## Sigma file (`ficsig`)

One 6x6 matrix per family, six lines of six values each, one matrix row per line, Fortran exponents allowed.
The matrix is always 6x6, for `espace = 6d` and `3d`.
Format, units and examples are in [Bunch generation](06_bunch.md).

## Bunch file (`ficloa`, `ficsto`)

Raw native float64, no record markers:

```
x y z      for each particle (m)
px py pz   for each particle (beta gamma)
```

`ficsto` is written at the end of the run: positions at step `nbdt - 1`, momenta at step `nbdt`.
`ficloa` is read after the sigma distribution, only when `iflcor != 0`.
These files are not compatible with the Fortran ones, which carry record markers.

## Marker file

Read by `lions.bidim.setup_bidims` when bidims are on; its name is asked on stdin (`no` for none).
It draws lines and circles in the trajectory bidim and defines ion traps.

```
line 1-3     ignored
droite                         straight lines
<count>
x1 y1 x2 y2                    count lines (mm)
cercle                         circles and arcs
<count>
xc yc r theta1 theta2 v flag   count lines (mm, mm, mm, deg, deg, unused, -)
fintyp                         end of file
```

`droite` and `cercle` blocks can appear in any order and number.
`droite` lines take plain numbers only; `cercle` lines accept Fortran exponents.
Each line or arc is drawn with 360 points.
Circle `flag`: below -10, the circle is also a trap (`trap_x`, `trap_y`, `trap_r`, in m); above -9000, it is drawn.
Any other keyword stops the program.

## Standard input

`lions` (`lions.cli.main`) reads its answers from stdin, one per line, in this order:

| # | Answer | Units | State attribute | Condition |
|---|---|---|---|---|
| 1 | byte swap of the central potential, `y` or `n` | | `run.swap_bytes` | |
| 2 | input deck path | | `run.deck_file` | |
| 3 | `xmin1 xmax1 ymin1 ymax1 nx1 ny1` | m, -, - | `central_*` | `ficv1 = no` |
| 4 | `rmin_ver rmax_ver teta1_ver teta2_ver` | m, m, deg, deg | `sorter.rmin`, `sorter.rmax`, `sorter.theta1`, `sorter.theta2` | |
| 5 | `pot_ver h_ver phi_ver anuz coupnuz dpotspot` | V, m, deg, -, m, - | `sorter.voltage`, `sorter.gap`, `sorter.phase`, `sorter.nuz`, `sorter.nuz_z_cut`, `sorter.dv` | |
| 6 | first bidim number, `0` for none | | `bidims.first` | |
| 7 | bidim period (steps) | | `bidims.every` | bidims on |
| 8 | `c` (orbit centres) or `b` (barycentres) | | `bidims.centre_mode` | bidims on |
| 9 | bidim box `xmin xmax ymin ymax` | m | `bidim_ax*`, `bidim_ay*` | bidims on |
| 10 | `dphi_max dw/w_max` | deg, - | | bidims on, n > 1 |
| 11 | `dr_max dpr/p0_max` | mm, mrad | | bidims on, n > 1 |
| 12 | `dz_max dpz/p0_max` | mm, mrad | | bidims on, n > 1 |
| 13 | first step stored in the bidims | | `bidims.start_step` | bidims on, n > 1 |
| 14 | barycentre window `xmin xmax ymin ymax` | m | | bidims on |
| 15 | marker file, `no` for none | | | bidims on |
| 16 | last step with traps active | | `traps.last_step` | traps defined |
| 17 | `0`: trapped particles stay; `1`: they follow particle 0 | | `traps.follow` | traps defined |
| 18 | write the list of trapped particles, `o` or `n` | | | after the run, traps defined |
| 19 | file name for the list | | | answer 18 is `o` |

`1000 2000 0 0` and `0 0 0 0 0 0` disable the vertical sorter: it is off whenever `h_ver = 0`.

Example, no bidim:

```bash
cd examples
printf 'n\ndon.120.ref_phi10_mars_2017\n1000 2000 0 0\n0 0 0 0 0 0\n0\n' | uv run lions
```

`lions.run(deck_file, swap=False, vertical=None, on_step=None)` takes the same inputs as arguments and runs without bidims.
`vertical` is a dict of `set_vertical_sorter` keywords, in the order of answers 4 and 5: `rmin`, `rmax` (m), `theta1`, `theta2` (deg), `voltage` (V), `gap` (m), `phase` (deg), `nuz`, `nuz_z_cut` (m), `dv`.
It asks nothing, except answer 3 when `ficv1 = no`.

## Input errors

| Situation | Behaviour |
|---|---|
| Missing deck, map, potential or sigma file | `FileNotFoundError` |
| Sum of `nbrion` differs from `n` | `ValueError` |
| Map grid outside 250 x 721 | `ValueError` |
| Full map ends before block 4 | `ValueError` |
| 1/8 map without blocks 3 and 4 | Message, dBz/dx = dBz/dy = 0 |
| Central potential header outside 500 x 300 | `ValueError` (check the byte swap) |
| Central potential with fewer than nx rows | `ValueError` |
| Central potential without d2pot | Message, d2pot = 0 |
| Outer potential grid outside 50 x 1000 or file ends early | `ValueError` |
| `b_factor` on a zero map | `ValueError` |
| Wrong bidim count or bit number in `BidimWriter.open` | Message, `sys.exit(1)` |
| Unknown keyword in the marker file | Message, `sys.exit(1)` |
