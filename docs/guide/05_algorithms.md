# Algorithms

The order of the calls in one step is given in [Overview](01_overview.md#call-sequence).
Formulas below use the state names of the Fortran code where they are shorter; the [glossary](../glossary.md) maps them to `State` attributes.

## Box and regions

`lions.apertures` sets the masks used by the field routines.

| Function | Mask | Condition |
|---|---|---|
| `check_box` | `bunch.in_box` | `rmin < r < rmax - 0.02 m` (map radii), `zmin < z < zmax`, `bunch.alive` |
| `check_central` | `bunch.in_central` | `bunch.in_box`, `r < rmax1`, `abs(y) < ymax1 - hy1` |
| `check_outer` | `bunch.in_outer` | (x', y') inside the outer grid, `bunch.in_box`, not `bunch.in_central` |

`check_box` also updates `bunch.r`.
Each function stores the count of its mask (`bunch.n_in_box`, `bunch.n_central`, `bunch.n_outer`).

## Magnetic field

`lions.fields.compute_fields`, first part.
Bilinear interpolation on the (r, theta) map, with the off-plane expansion to first order in z:

```
i = int((r - rmin) / hr)          j = int(theta / htet)
a = r - rmin - i hr               b = theta - j htet
d1 = a - hr                       d2 = b - htet

F(M) = [ d1 d2 M[i,j] - a d2 M[i+1,j] + a b M[i+1,j+1] - d1 b M[i,j+1] ] / (hr htet)

Bz = F(Bz) - z^2/2 F(d2Bz)
Bx = z F(dBz/dx)
By = z F(dBz/dy)
```

`hr` (m) and `htet` (rad) are the map steps (`field.dr`, `field.dtheta`).
Particles outside the box, and particles in the box where Bz = 0, get Bz = -1e-5 T and Bx = By = 0 (avoids a division by zero in the diagnostics).
For r <= rmin, Bz = `field.bz[0, 0]` and Bx = By = 0.
If x / r falls outside [-1, 1] by rounding, the count is printed and theta is clamped.

`lions.fields.compute_bfield` does the same interpolation for B only. `distribute_bunch` calls it to get Bz at the reference particles.

## Field harmonics

`lions.fields.add_harmonics`, called by `load` before the `b_factor` scaling.
On each map node, with (x, y) the node position (m), r^2 = x^2 + y^2, phi1 and phi2 the harmonic phases (rad), sums over the coefficients `coh[k]`:

First harmonic:

```
s       = x cos(-phi1) - y sin(-phi1)
Bz     += coh[0] s + SUM_{k>=1} coh[k] r^(2k) s
dBz/dx += coh[0] cos(-phi1) + SUM_{k>=1} coh[k] r^(2k-2) (2 k x s + r^2 cos(-phi1))
dBz/dy += -coh[0] sin(-phi1) + SUM_{k>=1} coh[k] r^(2k-2) (2 k y s - r^2 sin(-phi1))
d2Bz   += SUM_{k>=1} coh[k] r^(2k-2) 4 k (k+1) s
```

Second harmonic:

```
s       = (x^2 - y^2) cos(-phi2) - 2 x y sin(-phi2)
s2      = 2 x cos(-phi2) - 2 y sin(-phi2)
s3      = -2 y cos(-phi2) - 2 x sin(-phi2)
Bz     += coh[0] s + SUM_{k>=1} coh[k] r^(2k) s
dBz/dx += coh[0] s2 + SUM_{k>=1} coh[k] r^(2k-2) (2 k x s + r^2 s2)
dBz/dy += -coh[0] s3 + SUM_{k>=1} coh[k] r^(2k-2) (2 k y s + r^2 s3)
d2Bz   += SUM_{k>=1} coh[k] r^(2k-2) 4 k (k+2) s
```

The last theta column is then set equal to the first.
The d2Bz term is the Laplacian of the added polynomial.

## Central-region electric field

For particles in `bunch.in_central`, on the Cartesian grid of `pot1` (steps hx1, hy1):

```
i = int((x - xmin1) / hx1)        j = int((y - ymin1) / hy1)
a = x - xmin1 - i hx1             b = y - ymin1 - j hy1
p1 = (a - hx1)(b - hy1)   p2 = -a (b - hy1)   p3 = a b   p4 = -(a - hx1) b
```

At the four corners, centred differences with the off-plane term, for example at (i, j):

```
Ex1 = [ -(V[i+1,j] - V[i-1,j]) + z^2/2 (d2V[i+1,j] - d2V[i-1,j]) ] / (2 hx1)
Ey1 = [ -(V[i,j+1] - V[i,j-1]) + z^2/2 (d2V[i,j+1] - d2V[i,j-1]) ] / (2 hy1)
```

RF factor, with `rf.v_right`, `rf.v_left` the stored (sign-reversed) voltages and the RF angle reduced modulo 2 pi:

```
V_rf = rf.v_right cos(omega_rf t + rf.phase_right)     x >= 0
V_rf = rf.v_left  cos(omega_rf t + rf.phase_left)      x <  0
```

Fields:

```
Ex = V_rf (Ex1 p1 + Ex2 p2 + Ex3 p3 + Ex4 p4) / (hx1 hy1)
Ey = V_rf (Ey1 p1 + Ey2 p2 + Ey3 p3 + Ey4 p4) / (hx1 hy1)
Ez = z V_rf (d2V1 p1 + d2V2 p2 + d2V3 p3 + d2V4 p4) / (hx1 hy1)
```

The indices i - 1 and j - 1 are not checked: the grid must cover r < `rmax1`.

## Outer-region electric field

For particles in `bunch.in_outer`, in the gap frame (x', y') of [Coordinates](04_coordinates.md#gap-frame-outer-region), on the `pot2` grid (steps hx2, hy2), with a, b, d1, d2, p1..p4 as above and V1..V4, D1..D4 the potential and d2pot at (i,j), (i+1,j), (i+1,j+1), (i,j+1):

```
Ex' = -V_rf [ d2 (V1 - V2 - z^2/2 (D1 - D2)) + b (V3 - V4 + z^2/2 (D3 - D4)) ] / (hx2 hy2)
Ey' = -V_rf [ d1 (V1 - V4 - z^2/2 (D1 - D4)) - a (V2 - V3 + z^2/2 (D2 - D3)) ] / (hx2 hy2)
Ez  =  z V_rf (D1 p1 + D2 p2 + D3 p3 + D4 p4) / (hx2 hy2)
```

Rotation back to the laboratory frame:

```
Ex = Ex' cos(teta1) - Ey' sin(teta1)
Ey = Ex' sin(teta1) + Ey' cos(teta1)
```

Quadrant symmetries, `heps = rf.harmonic_sign` (+1 odd harmonic, -1 even):

```
x < 0, y > 0 :  Ex = heps Ex,   Ey = -heps Ey
x > 0, y < 0 :  Ey = -Ey
x < 0, y < 0 :  Ex = heps Ex,   Ey = heps Ey
x < 0        :  Ez = -heps Ez
```

## Vertical sorter

`lions.fields.add_vertical_sorter_field`, active when `h_ver != 0`, for `rmin_ver <= r <= rmax_ver` and `teta1_ver <= theta <= teta2_ver`:

```
Ez -= V_ver (1 - dpotspot) / (2 h_ver) cos(2 omega_rf t (1 + nuz / (2 harm)) + 2 phi_ver - pi/2)
Ez -= V_ver (1 + dpotspot) / (2 h_ver) cos(2 omega_rf t (1 - nuz / (2 harm)) + 2 phi_ver - pi/2)
```

`V_ver` (V) is `pot_ver`, `h_ver` (m) the gap, `nuz` = `anuz` the vertical tune of the drive.
When particle 0 is in the sorter region, `sorter.voltage_seen = -Ez h_ver` (V), written in the z file.

## Extraction elements

`lions.fields.add_extraction_fields`, after the sorter.

### Deflectors

`check_deflector(bunch, extraction, run, i)` uses two half-planes:

```
entrance : (x - xdfi) sin(tetdfi) - (y - ydfi) cos(tetdfi) <= 0
exit     : (x - xdfo) sin(tetdfo) - (y - ydfo) cos(tetdfo) >= 0
```

First deflector:

- bunch cut (`extraction.start_step = 0`): a live particle with r > `rdf` between the two lines enters it and is flagged `bunch.extracted`;
- whole bunch: every live particle between the two lines enters it from step `extraction.start_step`.

Particles already `bunch.extracted` stay in it while between the lines.
The other deflectors act only on `bunch.extracted` particles between their lines.
Inside a deflector, the electric field is replaced by a field of magnitude `edf` (V/m) normal to the transverse momentum:

```
Ex =  edf py / p_xy
Ey = -edf px / p_xy
```

### Dipoles

`check_dipole(bunch, extraction, i)`: `bunch.extracted` and live particles between the entrance and exit lines, both through (`xdm`, `ydm`).
Inside, with `d = sqrt((x - xdm)^2 + (y - ydm)^2)` (m) and `bdm` (T/m):

```
Bz += (d - rdm) bdm
Bx += z bdm (x - xdm) / d
By += z bdm (y - ydm) / d
```

### Grids

Tested when `rex < 0`, grids exist and `run.step >= extraction.grid_start_step`.
A particle is past grid i when

```
(x - xgex) uygex - (y - ygex) uxgex <= 0
```

- One grid, no deflector, no dipole: from step `extraction.start_step`, a live particle not yet flagged and past the grid is recorded. Its flag is reset once it is back on the other side, so it is recorded at each turn.
- Otherwise: a `bunch.extracted`, live particle is recorded at grid i when it has already crossed grid i - 1 (`bunch.grid_reached`), so the grids are crossed in order.

Each crossing writes one line to the grid file (see [Outputs](03_outputs.md#extraction-grid-files-ficdebgexn)).

## Boris pusher

`lions.pusher.boris_push`. With `p` in beta gamma, B in T, E in V/m, and the per-particle constants of `load`:

```
bunch.boris_b_coef = q dt / (2 m0)        (1/T)
bunch.boris_e_coef = q dt / (m0 c)        (m/V)
run.c_dt           = c dt                 (m)
```

One step, with gamma from the previous step:

```
t   = bunch.boris_b_coef / gamma * B                 (tx, ty, tz)
s   = 1 / (1 + tx^2 + ty^2 + tz^2)
f   = bunch.boris_e_coef * E                         full-step electric kick
v   = s (f + 2 p)

         | -(ty^2+tz^2)   tx ty + tz      tx tz - ty    |
M    =   |  tx ty - tz   -(tx^2+tz^2)     ty tz + tx    |
         |  tx tz + ty    ty tz - tx     -(tx^2+ty^2)   |

p_new     = p + M v + f
gamma_new = sqrt(1 + |p_new|^2)
x_new     = x + run.c_dt / gamma_new * p_new
```

This is the standard Boris scheme: half electric kick `f/2`, rotation about B, half electric kick.
The rotation uses the gamma of the previous step, not the gamma after the first half kick.
It is second order in dt and conserves |p| exactly when E = 0.

With `traps.follow = 0`, only live particles are pushed. With `traps.follow = 1` (the default), all particles are pushed; `apply_traps` has reset the trapped ones onto particle 0.

## Traps and losses

`lions.apertures.apply_traps`, at the start of each step. It returns at once when no trap is defined.
With traps:

- a particle with z < zmin or z > zmax is lost (`bunch.alive = False`);
- up to step `traps.last_step`, a particle inside a trap circle, `(x - xp)^2 + (y - yp)^2 <= rp^2`, is lost;
- with `traps.follow = 1`, lost particles take the position, momentum and gamma of particle 0.

Without traps, particles leaving the box are not flagged: they keep moving in Bz = -1e-5 T and E = 0.

## Correlations

Computed by `write_outputs` (and `phase_space_coords` for the bidims) for each family, relative to its first particle (index i1), with `p`, `rc` of that particle:

```
rc     = -m0 c p / (q Bz)
alpha  = asin(pr / p)
(dx, dy) = (x - x_i1, y - y_i1) rotated by (theta - alpha)
ds     = rc atan(dy / (rc + dx))
v      = sin(ds / (2 rc))^2
dr     = (dx + 2 rc v) / (1 - 2 v)
```

The momentum deviations are taken relative to the reference momentum rotated by `ds / rc`, then projected on the local radial and longitudinal directions:

```
dpr = radial component
dp  = sqrt((p + dp_long)^2 + dpr^2 + dpz^2) - p
dz  = z - z_i1,   dpz = pz - pz_i1
```

Sigma matrix over the live particles of the family:

```
sigma[i,j] = SUM_k alive_k ps[i,k] ps[j,k] / n_alive
```

Output quantities, gamma of the reference particle:

```
dr    = sqrt(8 sigma_11)                     (m, written in mm)
ds    = sqrt(8 sigma_55)                     (m, written in mm)
dz    = sqrt(8 sigma_33)                     (m, written in mm)
dphi  = ds omega_rf / (c p) gamma            (rad, written in deg)
epsr  = 8 sqrt|sigma_11 sigma_22 - sigma_12^2| 1e6 / p      (mm mrad)
epsz  = 8 sqrt|sigma_33 sigma_44 - sigma_34^2| 1e6 / p      (mm mrad)
dW/W  = sqrt(8 sigma_66) p / (gamma (gamma - 1))
dW    = dW/W * W                             (MeV/u)
dwdphi = omega_rf/(c p) gamma * p/(gamma (gamma-1)) * W
         * 8 sqrt|sigma_55 sigma_66 - sigma_56^2| * 1000 (deg)  (keV/u deg)
```

Indices 1 to 6 are dr, dpr, dz, dpz, ds, dp.
`write_outputs` also computes `bunch.sigma_dc = T sigma T^T`, with T the identity except `T[0,0] = T[1,1] = -1`, `T[0,5] = r / p`, `T[4,1] = -r / p`. It is kept in the state and not written.

`lions.beam.project_on_radial_line` and `projected_phase_space` compute deviations projected on the radial line through the reference particle (Fortran `calgri`, `calgric`). Nothing calls them: they belong to the `ficmat` output, not ported.

## Runtime messages

| Situation | Behaviour |
|---|---|
| x / r outside [-1, 1] by rounding | Count and particles printed, theta clamped |
| Bz = 0 in the box, or particle outside the box | Bz = -1e-5 T, Bx = By = 0 |
| Particle in a trap | `bunch.alive = False`; it stays, or follows particle 0 |
| Particle outside the central potential grid with r < `rmax1` | No check: wrong field (index wraps in NumPy) |
