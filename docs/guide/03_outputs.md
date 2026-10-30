# Outputs

All text outputs follow particle 0, the reference particle of family 0, unless stated.
They are written by `lions.diagnostics.write_outputs`, called at step 0 and after every step.
Files named in the deck are opened in write mode and overwritten.

## Output steps

- The turn counter `turn` is updated at every step: it increases by one when the azimuth of particle 0 decreases (forward tracking), and decreases by one when it increases (backward tracking, `npaspt < 0`).
- `theta_total_deg = theta_deg + 360 (turn - 1)` (deg).
- A line is written to the z, xy and correlation files when `step % ifreq == 0` or `step == 1`. Step 0 writes one too.
- Steps `nbdt - 1` and `nbdt` add the end-of-run values (see below).

## z file (`ficz`)

Header: two text lines (program name, then a label for the vertical motion).
One line per output step:

| Column | Quantity | Units |
|---|---|---|
| step | Time step | - |
| turn | Turn counter | - |
| theta | Azimuth of particle 0 | deg |
| x, y, z | Position | mm |
| r | Radius | mm |
| Bx, By, Bz | Magnetic field at the particle | T |
| Ez | Vertical electric field | kV/cm |
| tpz | Vertical angle `asin(pz / p)` | mrad |
| nuz | Local vertical tune estimate | - |
| Vver | Vertical sorter voltage seen, `-Ez * h_ver` | V |

`nuz` is computed only when `|z| > coupnuz` (m), otherwise 0:

```
nuz^2 = | -1/(m0 gamma) (harm/omega_rf)^2 (q / z) (Ez + c/gamma (px By - py Bx)) |
```

## xy file (`ficxy`)

Header: program name, deck title, then one value per line: atomic mass (u), A, q (e), r0 (mm), theta0 (deg), z0 (mm), pr0 (mrad), pz0 (mrad), w0 (MeV/u), phi0 (deg), fhf (Hz), harm, right and left dee voltages as in the deck (V), dB/B, `b_factor` (T, only when given), particle count, step count, steps per turn, dt (ns).
Then a column header line and one line per output step:

| Column | Quantity | Units |
|---|---|---|
| step | Time step | - |
| turn | Turn counter | - |
| theta | Azimuth | deg |
| x, y | Position | mm |
| xc, yc | Centre of curvature | mm |
| bz | Bz | T |
| es | E along the transverse momentum (tangential) | kV/cm |
| er | E across the transverse momentum (radial) | kV/cm |
| r | Radius | mm |
| tpr | Radial angle `asin(pr / p)`, `pr = (x px + y py) / r` | mrad |
| w | Kinetic energy per nucleon | MeV/u |
| dphi | RF phase lag | deg |

Definitions, with `p` in beta gamma and `p_xy = sqrt(px^2 + py^2)`:

```
xc   = x + py m0 c / (q Bz)
yc   = y - px m0 c / (q Bz)
es   = (Ex px + Ey py) / p_xy / 1e5
er   = (Ex py - Ey px) / p_xy / 1e5
w    = m0 c^2 (gamma - 1) / (MEV * A)
dphi = omega_rf t (deg) - harm * theta_total_deg - 90
```

### End-of-run block

At the last step the xy file receives:

```
---------------------------------
start point for reverse tracking:
r0 theta0 z0 pr0 pz0 w0mev phi0 :
<r> <theta> <z> <pr>              values saved at step nbdt - 1
<pz> <w> <phi>
start point to continue tracking:
<r> <theta> <z> <pr>              values at step nbdt
<pz> <w> <phi>
input for DEDUCE_SIGMA:
amat a qn w0mev :
<amat> <a> <qn> <w>
magnetic field (absolute value):
<-Bz>
fhf harmonic:
<fhf> <harm>
```

Units are those of deck line 13: r and z in m, theta and phi in deg, pr and pz in mrad, w in MeV/u, so the two start points can be pasted into a deck.
phi is `omega_rf t` in deg, not reduced modulo 360. -Bz is in T.

## Correlation files (`ficcor`, `ficmat`, `ficdc`)

Written when `iflcor != 0` and at least one of the three names is not `no`, at the same steps as the z and xy lines, one line per family.

`ficcor` columns (header line written by `load`):

| Column | Quantity | Units |
|---|---|---|
| tr | Turn counter | - |
| r | Radius of the family reference particle | mm |
| teta | Its azimuth | deg |
| dr | Full radial width, `sqrt(8 sigma_11)` | mm |
| ds | Full longitudinal width, `sqrt(8 sigma_55)` | mm |
| dphi | Full phase width | deg |
| epsr | Radial emittance, `8 sqrt(det) / p` | mm mrad |
| dwsw | Full relative energy spread dW/W | - |
| wmpa | Energy per nucleon | MeV/u |
| dw | Full energy spread | MeV/u |
| dwdphi | Longitudinal emittance | keV/u deg |
| dz | Full vertical width, `sqrt(8 sigma_33)` | mm |
| epsz | Vertical emittance | mm mrad |

The factor 8 converts rms values into the full extent of a uniform 6D ball (`lions.constants.SQRT8`).
Formulas are in [Algorithms](05_algorithms.md#correlations).
The sigma matrix used here is the second moment over the live particles of the family, divided by `n_alive`.

`ficmat` and `ficdc` are opened (`ficdc` with a header line) and nothing else is written to them: the matrix outputs of the Fortran are not ported.

## Extraction grid files (`ficdebgex.N`)

One file per grid, opened by `load`.
Header: a label line, then `tetagex xgex ygex xgexc ygexc` (deg, m, m, m, m); a particle line `n A amat q`; an RF line `npaspt fhf harm`.
Then one line per crossing:

| Column | Quantity | Units |
|---|---|---|
| i | Particle number, from 1 | - |
| step | Time step | - |
| x, y, z | Position | m |
| px, py, pz | Momentum | beta gamma |
| Bz | Magnetic field | T |

## Bidim files (`N.bid`)

2D histograms in the GANIL binary format, written by `lions.bidim.BidimWriter`.
File names are the bidim numbers: `numbid.bid`, `numbid+1.bid`, ...
Each axis has 13 bits (channels 0 to 8191).

Header, 351 bytes, native byte order:

| Bytes | Content |
|---|---|
| 4 | bidim number (int32) |
| 10 | type `bidim` |
| 12 + 12 | date and time (3 int32 each, zero) |
| 10 + 10 | place and author |
| 25 | title |
| 25 + 25 | x and y labels |
| 8 | x and y channel maxima (2 int32) |
| 10 | format `Binaire` |
| 4 x 20 | xmin, xmax, ymin, ymax as text |
| 120 | blank |

Data: blocks of 4096 int32 events, each `ix * 65536 + iy`.
A partly filled last block ends with the marker `9999 * 65536 + 9999`.

Histograms, numbered from the first bidim number `nb`:

| Number | Content | Axes |
|---|---|---|
| nb | Trajectory of particle 0, plus marker-file lines and circles | x, y in the bidim box (m) |
| nb + 1 | Orbit centres (`c`), or their barycentre over one turn (`b`) | x, y in the box (`c`) or barycentre window (`b`) |
| nb + 2 | Trajectories of all particles | x, y in the bidim box |
| nb + 3 | dphi vs dW/W | +-dphi_max (deg), +-dw/w_max |
| nb + 4 | dr vs dpr/p0 | +-dr_max (mm), +-dpr/p0_max (mrad) |
| nb + 5 | dr vs dphi | +-dr_max (mm), +-dphi_max (deg) |
| nb + 6 | step vs dr | 0 to nbdt, +-dr_max (mm) |
| nb + 7 | dz vs dpz/p0 | +-dz_max (mm), +-dpz/p0_max (mrad) |
| nb + 8 | step vs dz | 0 to nbdt, +-dz_max (mm) |

With one particle, 5 files are opened and only nb and nb + 1 are filled.
With more than one particle, 17 files are opened; nb + 3 to nb + 8 are filled only with a single family.
Events are stored from step `bidim_start_step`, at step 1 and every `bidim_every` steps.
Out-of-range positions are dropped for nb to nb + 2 and clipped to the edge for the others.

## Bunch file (`ficsto`)

See [Inputs](02_inputs.md#bunch-file-ficloa-ficsto).

## Trapped-particle list

With traps defined, `lions.cli.write_trapped` offers to write the trapped particles: one line per particle, `count index` (index from 1).

## Console

The program prints the banner, the prompts, the reading steps (map type, harmonics, `b_factor` scaling, extraction mode, bunch distribution) and error messages.
The GUI shows the same text in its console.
