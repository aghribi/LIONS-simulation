# Worked example

`examples/don.120.ref_phi10_mars_2017`, a 100Mo22+ bunch in CIME (GANIL).
The field maps and the sigma file it names are not in the repository.

## Deck

```
 "don.forw.H3 i700 posbmax "
1080
1.50d0 250
14400000.d0 -51890.d0 20.d0 3.d0 +1.0d-4
0 0
51890.d0 51890.d0
0
0
-.5d0 .5d0 .164d0
-1.  -1. -1.  -1.
1
1.513  6990000 0 0 0 0 250 267
0

1
100 100.36496d0 100.36496d0 22.d0
0.03316181934576   266.8973516728 0.d0   426.2811141832   0.d0 5.6650d-3   -250.00d0
1
sigma6d_no_corr.dat
6d
no no no
350000 30
exp.for.z.ref_phi10_mars_2017
exp.for.ref_phi10_mars_2017

'cime700rec.1804.pbob1.lio.dat'
'pot3d.impair.30pilbleu30.620.dat'
'newgap612.dat'
```

## Reading

| Line | Value | Meaning |
|---|---|---|
| 2 | 1080 | steps per turn |
| 3 | 1.50 m, 250 | design final radius and turns (not used) |
| 4 | 14.4 MHz, harmonic 3, dee half-angle 20 deg, dB/B = +1e-4 | `pothf` = -51890 V is ignored |
| 5 | 0, 0 | dee phase offsets (deg) |
| 6 | 51890 V, 51890 V | right and left dee voltages |
| 7, 8 | 0, 0 | no field harmonics |
| 9 | -0.5 m, 0.5 m, 0.164 m | vertical box, central region radius |
| 10 | -1 -1 -1 -1 | extraction on (`rex < 0`), bunch cut (`tet1ex <= 0`), no grid (`tet2ex <= 0`) |
| 10b, 10c | 1 deflector | entry radius 1.513 m, 6.99 MV/m, both lines through the origin, at 250 deg and 267 deg |
| 10d | 0 | no dipole |
| 11 | 1 | one family |
| 12 | 100 ions, 100.36496 u, A = 100.36496, q = 22 | |
| 13 | r0 = 33.16 mm, theta0 = 266.90 deg, z0 = 0, pr0 = 426.28 mrad, pz0 = 0, 5.665 keV/u, phi0 = -250 deg | |
| 14 | `iflcor = 1`, `sigma6d_no_corr.dat`, `6d`, no correlation output | bunch of 100 around particle 0 |
| 15 | 350000 steps, output every 30 | |
| 16, 17 | z and xy outputs | |
| 18 to 20 | field map, central and outer potentials | no `b_factor`: the map is used as read |

Derived values:

| Quantity | Value |
|---|---|
| dt = harm / (fhf npaspt) | 0.1929 ns |
| Revolution period harm / fhf | 208.3 ns |
| Tracked time 350000 dt | 67.5 us |
| Turns 350000 / 1080 | 324 |
| Output lines per turn 1080 / 30 | 36 |
| Rest mass `AMU (100.36496 - 0.0005486 * 22)` | 1.66640e-25 kg |
| Initial beta | 3.49e-3 |

With an isochronous field, `beta = omega_rev r / c` gives 10.8 MeV/u at r = 1.50 m.
This is an estimate from the deck, not a tracking result.

## Run

Command line, no bidim, no vertical sorter (the map files must be present in `examples/`):

```bash
cd examples
printf 'n\ndon.120.ref_phi10_mars_2017\n1000 2000 0 0\n0 0 0 0 0 0\n0\n' | uv run lions
```

Python:

```python
import lions

d = lions.run('don.120.ref_phi10_mars_2017', swap=False)
b = d.bunch
b.x, b.y, b.z, b.px, b.py, b.pz, b.alive   # per particle, SI, p = beta gamma
```

GUI:

1. `uv run lions-gui`.
2. In General, choose the deck and press Load data file: the form fills.
3. Check the parameters, press Run.

## Expected behaviour

- XY trajectory: a spiral from the injection radius (33 mm) outwards.
- R - theta: the radius grows in steps, one per gap crossing.
- Energy vs turn: the energy per nucleon grows from 5.665 keV/u.
- The z and xy files get one line every 30 steps, and the end-of-run block at step 350000.
