# Bunch generation

## Principle

The deck field `iflcor` (`bunch.correlations`) controls the initial bunch.

- `iflcor = 0`: all particles of a family start on the family initial conditions (deck line 13).
- `iflcor != 0`: the particles of each family are spread around its first particle with the 6x6 sigma matrix of that family, read from `ficsig`.

The first particle of each family keeps the deck initial conditions: it is the reference particle.
`lions.beam.distribute_bunch` does the work, called by `load` once the magnetic map and potentials are read.
Families with one particle are skipped.

## Deck lines

```
iflcor                 != 0 to distribute
ficsig                 sigma file
espace                 6d or 3d
ficcor ficmat ficdc    correlation outputs, no = off
```

## Sigma file

One 6x6 matrix per family, in family order, six lines of six values, one matrix row per line:

```
s(dr,dr)   s(dr,dpr)   s(dr,dz)   s(dr,dpz)   s(dr,ds)   s(dr,dp)
s(dpr,dr)  s(dpr,dpr)  ...
s(dz,dr)   ...
s(dpz,dr)  ...
s(ds,dr)   ...
s(dp,dr)   ...                                            s(dp,dp)
```

The order is that of `bunch.ps`: dr, dpr, dz, dpz, ds, dp (see [Coordinates](04_coordinates.md)).
The matrix is the covariance matrix of the deviations: diagonal terms are variances.
It must be symmetric positive definite (Cholesky factorisation).

Units of each element are the product of the units of its two coordinates:

| Coordinate | Units |
|---|---|
| dr | m |
| dpr | beta gamma |
| dz | m |
| dpz | beta gamma |
| ds | m (arc length) |
| dp | beta gamma |

So `s(dr,dr)` is in m^2, `s(dr,dpr)` in m, `s(ds,ds)` in m^2, `s(dp,dp)` dimensionless.

### Diagonal matrix

Without correlations the matrix is diagonal.
The square root of a diagonal term is the rms spread of that coordinate:

```
1.0d-6   0.0d0    0.0d0    0.0d0    0.0d0    0.0d0
0.0d0    1.0d-6   0.0d0    0.0d0    0.0d0    0.0d0
0.0d0    0.0d0    1.0d-6   0.0d0    0.0d0    0.0d0
0.0d0    0.0d0    0.0d0    1.0d-6   0.0d0    0.0d0
0.0d0    0.0d0    0.0d0    0.0d0    1.0d-6   0.0d0
0.0d0    0.0d0    0.0d0    0.0d0    0.0d0    1.0d-6
```

Here `s(dr,dr) = 1e-6 m^2`: rms radial spread 1 mm; `s(ds,ds) = 1e-6 m^2`: rms arc spread 1 mm; `s(dpr,dpr) = 1e-6`: rms dpr of 1e-3 in beta gamma.

## Algorithm

### Step 1: uniform 6D ball

`lions.beam.sample_unit_ball(6, n)`.
Particles 1 to n - 1 are drawn uniformly in the unit 6D ball by rejection:

```
repeat for the rejected particles:
    v[j] = 2 * random() - 1        j = 0..5
until SUM_j v[j]^2 <= 1 for all
```

Particle 0 is set to 0 (the reference).
The NumPy legacy generator is seeded with 42 at each call, so a run is reproducible.
The Fortran uses the gfortran generator: bunches differ particle by particle between the two codes.

### Step 2: unit covariance

The sample covariance `C = (1/n) SUM_k v_k v_k^T` is factorised, `C = R^T R` (`lions.beam.cholesky_upper`), and the sample is transformed by forward substitution, `v <- (R^T)^-1 v`.
The sample covariance is then the identity.

`cholesky_upper(m, a)` returns the upper triangular R with `a = R^T R`:

```
R[0,0] = sqrt(A[0,0])
R[0,j] = A[0,j] / R[0,0]
R[i,i] = sqrt(A[i,i] - SUM_{k<i} R[k,i]^2)
R[i,j] = (A[i,j] - SUM_{k<i} R[k,i] R[k,j]) / R[i,i]      j > i
```

### Step 3: apply the sigma matrix

```
sigma = R_s^T R_s               (cholesky_upper of the file matrix)
ps    = R_s^T v
```

so that `Cov(ps) = R_s^T I R_s = sigma`.
The distribution is a uniform ellipsoid: since a uniform 6D ball of unit covariance has radius sqrt(8), each full width is `sqrt(8)` times the rms value.
With a diagonal sigma, `ps[i] = sqrt(sigma_ii) v[i]`.

### Step 4: place the particles

With `rc = -m0 c p / (q Bz)` of the reference particle (Bz from `compute_bfield`):

1. If the family starts in the outer region (`r0 > rmax1`), dr is corrected by `f = sqrt(cos(harm ds / rc))`: `dr <- dr f - rc (1 - f)`. A message says whether the correction is applied.
2. Curvilinear to Cartesian, in the frame of the reference orbit:
   ```
   dx  = dr - 2 (rc + dr) sin(ds / (2 rc))^2
   dy  = (rc + dr) sin(ds / rc)
   dp_long = sqrt((p + dp)^2 - dpr^2 - dpz^2) - p
   dpx = dpr cos(ds/rc) - dp_long sin(ds/rc)
   dpy = dpr sin(ds/rc) + dp_long cos(ds/rc)
   ```
3. Rotate (dx, dpx, dy, dpy) by `teta0 - pr0` into the laboratory frame.
4. Add the position offsets; rotate the reference momentum (px, py) by `ds / rc` and add (dpx, dpy).
5. Add dz and dpz to z and pz.
6. In the outer region, multiply px, py, pz by `f`.

## 6d and 3d

The sigma file always holds 6x6 matrices.
In the Python code both `espace = 6d` and `espace = 3d` add the vertical deviations dz, dpz (step 5).
The Fortran adds them only for `3d`; with `6d` it leaves z and pz on the reference values.

## Loaded bunch

With `iflcor != 0` and `ficloa != no`, the positions and momenta of all particles are then replaced by those of the bunch file (format in [Inputs](02_inputs.md)).

## GUI

The **Correlations** section of the GUI form sets `iflcor`, the sigma file, `espace` and the three output names.
It also offers six diagonal spreads, in GUI units:

| Field | Coordinate | GUI units | Default |
|---|---|---|---|
| sigma dr | dr | mm | 0.1 |
| sigma ds | dpr | mrad | 1.0 |
| sigma dz | dz | mm | 0.1 |
| sigma dpz | dpz | mrad | 1.0 |
| sigma dphi | ds | mm | 0.1 |
| sigma dp | dp | mrad | 1.0 |

The GUI labels the second coordinate `ds` and the fifth `dphi`; they are dpr and ds of `bunch.ps`.
**Apply -> sigma file** writes a diagonal 6x6 matrix (one family) from these values; Run does the same when `iflcor != 0`.
Conversion to the file:

```
position:  variance = (sigma_mm * 1e-3)^2                      (m^2)
momentum:  variance = (sigma_mrad * 1e-3 * p0)^2               (beta gamma)^2
```

`p0` (beta gamma) is computed by the form from the family 0 energy as `w0mev * amat` (the engine uses `w0mev * a`).
Each variance is floored at 1e-30 so that the Cholesky factorisation succeeds.
When a deck with `iflcor != 0` is loaded, the six diagonal terms of the sigma file are read back with the inverse conversion.

### Example

Spreads 0.5 mm in dr, 0.2 mm in dz, 0.3 mm in ds, momenta 0:

```
sigma dr   = 0.5   (mm)
sigma ds   = 0     (mrad)
sigma dz   = 0.2   (mm)
sigma dpz  = 0     (mrad)
sigma dphi = 0.3   (mm)
sigma dp   = 0     (mrad)
```

gives the diagonal 2.5e-7, 1e-30, 4e-8, 1e-30, 9e-8, 1e-30 (0.5 mm squared is 2.5e-7 m^2, 0.2 mm is 4e-8 m^2, 0.3 mm is 9e-8 m^2).

For bidims with correlation plots, set more than one ion and a non-zero first bidim number (see [GUI](07_gui.md)).
