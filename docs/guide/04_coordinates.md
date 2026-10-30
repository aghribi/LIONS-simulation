# Coordinates, units and constants

## Laboratory frame

- (x, y): median plane of the cyclotron, origin on the machine axis (m).
- z: vertical, normal to the median plane (m).
- (r, theta): polar coordinates in the median plane. theta is in [0, 2 pi), measured from x towards y:

```
theta = acos(x / r)            if y >= 0
theta = 2 pi - acos(x / r)     if y < 0
```

Particles turn counter-clockwise (theta increasing) for Bz < 0.
The right dee covers x >= 0, the left dee x < 0.

## Gap frame (outer region)

The outer potential map is stored in a frame (x', y') rotated by the dee half-angle `teta1` and shifted by the skewed-gap offset `yomega`:

```
x' =  |x| cos(teta1) + (|y| - yomega) sin(teta1)
y' = -|x| sin(teta1) + (|y| - yomega) cos(teta1)
```

The map covers one quadrant. The field computed in (x', y') is rotated back by `teta1` and the other quadrants follow from the symmetries given in [Algorithms](05_algorithms.md#outer-region-electric-field).

## Units in memory

The `State` holds SI units, angles in rad, momenta normalised, unless the attribute name says otherwise (`run.theta_deg`, `bunch.families[k].phi0_deg`, `bunch.w0_mev_per_u`, `extraction.grids[k].angle_deg`).

| Quantity | Units |
|---|---|
| Positions | m |
| Time | s |
| Momenta `bunch.px`, `bunch.py`, `bunch.pz`, `bunch.p` | beta gamma = p / (m0 c), dimensionless |
| Magnetic field | T |
| Electric field | V/m |
| Mass, charge | kg, C |
| Energy `bunch.families[k].w0`, `bunch.rest_energy` | J |

Deck inputs and text outputs use other units: see [Inputs](02_inputs.md) and [Outputs](03_outputs.md).

### Normalised momentum

The momentum is stored as `p = beta gamma`, with `gamma = 1 / sqrt(1 - beta^2) = sqrt(1 + p^2)`.
The velocity is `v = c p / gamma`.
The magnetic rigidity is `B rho = m0 c (beta gamma) / q` (T m).

### Angles from momenta

Plots and outputs give momenta as angles (mrad):

```
tpz = asin(pz / p) * 1000
tpr = asin(pr / p) * 1000,   pr = (x px + y py) / r
```

For small angles, `tpz ~ pz / p * 1000`.
Example: A = 100, q = 22, W = 0.27 MeV/u gives p = 0.024; a spread of 3.16e-6 in pz is 0.13 mrad.

## Phase space relative to the reference particle

Correlations, bidims and the bunch distribution use the deviations of each particle from the first particle of its family, in the curvilinear frame of that particle.
They are stored in `bunch.ps` (6 x n), computed by `lions.beam.phase_space_coords` and by `write_outputs`.

| Index | Symbol | Quantity | Units |
|---|---|---|---|
| 0 | dr | Radial deviation from the reference orbit | m |
| 1 | dpr | Radial momentum deviation | beta gamma |
| 2 | dz | Vertical deviation | m |
| 3 | dpz | Vertical momentum deviation | beta gamma |
| 4 | ds | Arc length along the reference orbit | m |
| 5 | dp | Deviation of the momentum magnitude | beta gamma |

The three planes are (dr, dpr), (dz, dpz) and (ds, dp).
Positions are in m so that the sigma matrix is homogeneous and its 2x2 determinants give emittances.

`ds` is an arc length, not an angle.
The angle seen from the centre of curvature is `ds / rc` (rad), with the radius of curvature

```
rc = -m0 c p / (q Bz)        (m, positive for Bz < 0 and q > 0)
```

The RF phase spread is `dphi = ds * omega_rf / (c p) * gamma` (rad), that is `ds / v * omega_rf`.
Example: at r = 33 mm, ds = 0.1 mm is 0.003 rad (0.17 deg) of azimuth.

`bunch.ps_cartesian` holds the same deviations in the Cartesian frame rotated by `theta - alpha` (`alpha = asin(pr / p)` of the reference particle), in the index order x 0, px 1, z 2, pz 3, y 4, py 5.

## Physical constants

`lions.constants`, values of the Fortran reference (CODATA 1986):

| Name | Quantity | Value | Units |
|---|---|---|---|
| `C_LIGHT` | Speed of light | 299792458 | m/s |
| `E_CHARGE` | Elementary charge | 1.60217733e-19 | C |
| `MEV` | 1 MeV | 1.60217733e-13 | J |
| `AMU` | Atomic mass unit | 1.6605402e-27 | kg |
| `ELECTRON_MASS_U` | Electron mass | 0.0005486 | u |
| `SQRT8` | rms to full width of a uniform 6D ball | 2.828427 | - |

## Key formulas

Rest mass of the ion, electrons removed:

```
m0 = AMU * (amat - ELECTRON_MASS_U * qn)          (kg)
```

Kinetic energy:

```
W = m0 c^2 (gamma - 1)                            (J)
W per nucleon = W / (MEV * A)                     (MeV/u)
```

Cyclotron frequency and RF:

```
f_rev = q |B| / (2 pi m0 gamma)                   (Hz)
f_rf  = harm * f_rev
```

Time step:

```
dt = harm / (fhf * npaspt)                        (s)
```

Radius of curvature: `rc = -m0 c p / (q Bz)` (m).
