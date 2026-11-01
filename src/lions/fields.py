"""Magnetic and electric fields at the particles."""


import numpy as np

from .apertures import check_central, check_deflector, check_dipole, check_grid, check_outer


def add_harmonics(rmin, rmax, nr, ntheta, bz, d2bz, dbzdx, dbzdy, n_coef, coef, phase_deg):
    """Add first and second field harmonics to the map arrays, in place."""
    pi = np.pi
    dr = (rmax - rmin) / (nr - 1.0)
    dtheta = 360.0 / (ntheta - 1.0) * pi / 180.0
    phase1 = phase_deg[0] * pi / 180.0
    phase2 = phase_deg[1] * pi / 180.0

    if n_coef[0] != 0:
        for i in range(nr):
            r = rmin + i * dr
            for j in range(ntheta - 1):
                theta = j * dtheta
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                r2 = x**2 + y**2
                co1 = np.cos(-phase1)
                si1 = np.sin(-phase1)
                s = x * co1 - y * si1
                pb = coef[0, 0] * s
                px_val = coef[0, 0] * co1
                py_val = -coef[0, 0] * si1
                p2 = 0.0
                for k in range(1, n_coef[0]):
                    k1 = k
                    k2 = k - 1
                    r2pk1 = r2**k1
                    if k == 1 and r2 <= 0.0:
                        r2pk2 = 1.0
                    else:
                        r2pk2 = r2**k2
                    pb += coef[k, 0] * r2pk1 * s
                    px_val += coef[k, 0] * r2pk2 * (2.0 * x * k1 * s + r2 * co1)
                    py_val += coef[k, 0] * r2pk2 * (2.0 * y * k1 * s - r2 * si1)
                    p2 += coef[k, 0] * r2pk2 * 4.0 * s * k1 * (k + 1)
                bz[i, j] += pb
                dbzdx[i, j] += px_val
                dbzdy[i, j] += py_val
                d2bz[i, j] += p2

    if n_coef[1] != 0:
        for i in range(nr):
            r = rmin + i * dr
            for j in range(ntheta - 1):
                theta = j * dtheta
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                r2 = x**2 + y**2
                co2 = np.cos(-phase2)
                si2 = np.sin(-phase2)
                s = (x**2 - y**2) * co2 - 2.0 * x * y * si2
                s2 = 2.0 * x * co2 - 2.0 * y * si2
                s3 = -2.0 * y * co2 - 2.0 * x * si2
                pb = coef[0, 1] * s
                px_val = coef[0, 1] * s2
                py_val = -coef[0, 1] * s3
                p2 = 0.0
                for k in range(1, n_coef[1]):
                    k1 = k
                    k2 = k - 1
                    r2pk1 = r2**k1
                    if k == 1 and r2 <= 0.0:
                        r2pk2 = 1.0
                    else:
                        r2pk2 = r2**k2
                    pb += coef[k, 1] * r2pk1 * s
                    px_val += coef[k, 1] * r2pk2 * (2.0 * x * k1 * s + r2 * s2)
                    py_val += coef[k, 1] * r2pk2 * (2.0 * y * k1 * s + r2 * s3)
                    p2 += coef[k, 1] * r2pk2 * 4.0 * s * k1 * (k + 2)
                bz[i, j] += pb
                dbzdx[i, j] += px_val
                dbzdy[i, j] += py_val
                d2bz[i, j] += p2

    # Periodicity in theta
    bz[:, ntheta - 1] = bz[:, 0]
    dbzdx[:, ntheta - 1] = dbzdx[:, 0]
    dbzdy[:, ntheta - 1] = dbzdy[:, 0]
    d2bz[:, ntheta - 1] = d2bz[:, 0]


def prepare_hermite(field, ntheta_all):
    """Node slopes of the bicubic Hermite option from the map: dBz/dr, dBz/dtheta, d2Bz/drdtheta.

    The cross derivative is a second-order difference of dBz/dtheta in r (one-sided at the map edges).
    """
    nr = field.nr
    r = field.rmin + field.dr * np.arange(nr)[:, None]
    theta = field.dtheta * np.arange(ntheta_all)[None, :]
    gx, gy = field.dbzdx[:nr, :ntheta_all], field.dbzdy[:nr, :ntheta_all]
    field.dbzdr, field.dbzdtheta, field.d2bzdrdtheta = (np.zeros_like(field.bz) for _ in range(3))
    field.dbzdr[:nr, :ntheta_all] = np.cos(theta) * gx + np.sin(theta) * gy
    field.dbzdtheta[:nr, :ntheta_all] = r * (np.cos(theta) * gy - np.sin(theta) * gx)
    field.d2bzdrdtheta[:nr, :ntheta_all] = np.gradient(field.dbzdtheta[:nr, :ntheta_all], field.dr, axis=0,
                                                       edge_order=2)


def hermite_bz(field, ivec, jvec, veca, vecb, r, theta):
    """Bicubic Hermite Bz in the cell (ivec, jvec) at offsets (veca, vecb), and its Cartesian gradient.

    Returns Bz, dBz/dx, dBz/dy on the median plane; C1 across cells.
    """
    u, v = veca / field.dr, vecb / field.dtheta

    def basis(s):
        s2, s3 = s * s, s * s * s
        return ((2 * s3 - 3 * s2 + 1, -2 * s3 + 3 * s2), (s3 - 2 * s2 + s, s3 - s2),
                (6 * s2 - 6 * s, -6 * s2 + 6 * s), (3 * s2 - 4 * s + 1, 3 * s2 - 2 * s))

    hu, su, dhu, dsu = basis(u)
    hv, sv, dhv, dsv = basis(v)
    dr, dt = field.dr, field.dtheta
    bz = ds_r = ds_t = 0.0
    for a in (0, 1):
        for b in (0, 1):
            ik, jk = ivec + a, jvec + b
            f, fr = field.bz[ik, jk], field.dbzdr[ik, jk]
            ft, frt = field.dbzdtheta[ik, jk], field.d2bzdrdtheta[ik, jk]
            bz = bz + f * hu[a] * hv[b] + fr * dr * su[a] * hv[b] + ft * dt * hu[a] * sv[b] \
                + frt * dr * dt * su[a] * sv[b]
            ds_r = ds_r + f * dhu[a] * hv[b] / dr + fr * dsu[a] * hv[b] + ft * dt * dhu[a] * sv[b] / dr \
                + frt * dt * dsu[a] * sv[b]
            ds_t = ds_t + f * hu[a] * dhv[b] / dt + fr * dr * su[a] * dhv[b] / dt + ft * hu[a] * dsv[b] \
                + frt * dr * su[a] * dsv[b]
    c, s = np.cos(theta), np.sin(theta)
    with np.errstate(divide='ignore', invalid='ignore'):   # r = 0 only outside the map, overwritten by the callers
        return bz, c * ds_r - s * ds_t / r, s * ds_r + c * ds_t / r


def compute_fields(d):
    """Magnetic and electric fields at all particles for the current time."""
    hrhtetm1 = 1.0 / (d.field.dr * d.field.dtheta)
    hrm1 = 1.0 / d.field.dr
    htetm1 = 1.0 / d.field.dtheta
    hx1hy1m1 = 1.0 / (d.central.dx * d.central.dy)
    hx1m1 = 1.0 / d.central.dx
    hx2hy2m1 = 1.0 / (d.outer.dx * d.outer.dy)
    hx2m1 = 1.0 / d.outer.dx
    hy1m1 = 1.0 / d.central.dy
    hy2m1 = 1.0 / d.outer.dy
    ivec = np.zeros(d.bunch.n_particles, dtype=np.int32)
    jvec = np.zeros(d.bunch.n_particles, dtype=np.int32)
    vec1 = np.zeros(d.bunch.n_particles)
    vec2 = np.zeros(d.bunch.n_particles)
    vec3 = np.zeros(d.bunch.n_particles)
    vec4 = np.zeros(d.bunch.n_particles)
    veca = np.zeros(d.bunch.n_particles)
    vecb1 = np.zeros(d.bunch.n_particles)
    vecb2 = np.zeros(d.bunch.n_particles)
    vecb3 = np.zeros(d.bunch.n_particles)
    vecb4 = np.zeros(d.bunch.n_particles)
    vecb = np.zeros(d.bunch.n_particles)
    vecd1 = np.zeros(d.bunch.n_particles)
    vecd2 = np.zeros(d.bunch.n_particles)
    vecp1 = np.zeros(d.bunch.n_particles)
    vecp2 = np.zeros(d.bunch.n_particles)
    vecp3 = np.zeros(d.bunch.n_particles)
    vecp4 = np.zeros(d.bunch.n_particles)
    vecpc = np.zeros(d.bunch.n_particles)
    vecx = np.zeros(d.bunch.n_particles)
    vecy = np.zeros(d.bunch.n_particles)
    vecz = np.zeros(d.bunch.n_particles)
    vex1 = np.zeros(d.bunch.n_particles)
    vex2 = np.zeros(d.bunch.n_particles)
    vex3 = np.zeros(d.bunch.n_particles)
    vex4 = np.zeros(d.bunch.n_particles)
    vey1 = np.zeros(d.bunch.n_particles)
    vey2 = np.zeros(d.bunch.n_particles)
    vey3 = np.zeros(d.bunch.n_particles)
    vey4 = np.zeros(d.bunch.n_particles)
    n = d.bunch.n_particles

    # B from the median-plane map, first order in z
    n_bad = np.count_nonzero((d.bunch.x / d.bunch.r > 1.0) | (d.bunch.x / d.bunch.r < -1.0))
    if n_bad != 0:
        print(f'particles with |x/r| > 1: {n_bad}')
        for i in range(n):
            ratio = d.bunch.x[i] / d.bunch.r[i]
            if ratio > 1.0 or ratio < -1.0:
                print(f'{i} {d.bunch.x[i]} {d.bunch.y[i]} {d.bunch.r[i]} {ratio}')

    veca = d.bunch.x / d.bunch.r
    mask_normal = (veca > -1.0) & (veca < 1.0)
    mask_ge1 = veca >= 1.0
    mask_lem1 = veca <= -1.0
    d.bunch.theta[mask_normal] = np.arccos(veca[mask_normal])
    d.bunch.theta[mask_ge1] = np.arccos(1.0)
    d.bunch.theta[mask_lem1] = np.arccos(-1.0)
    mask_yneg = d.bunch.y < 0.0
    d.bunch.theta[mask_yneg] = 2.0 * np.pi - d.bunch.theta[mask_yneg]
    # y < 0 with x / r rounding to 1 gives exactly 2 pi, one column past the map
    d.bunch.theta[d.bunch.theta >= 2.0 * np.pi] -= 2.0 * np.pi

    ivec[:] = 0
    jvec[:] = 0
    mask = d.bunch.in_box
    ivec[mask] = ((d.bunch.r[mask] - d.field.rmin) * hrm1).astype(np.int32)
    jvec[mask] = (d.bunch.theta[mask] * htetm1).astype(np.int32)
    veca[:] = d.bunch.r - d.field.rmin - ivec * d.field.dr
    vecb[:] = d.bunch.theta - jvec * d.field.dtheta
    vecd1[:] = veca - d.field.dr
    vecd2[:] = vecb - d.field.dtheta
    vecz[:] = d.bunch.z * d.bunch.z * 0.5

    # 2x2 cell: M[a, b] = map[i + a, j + b]
    ik = ivec + np.arange(2)[:, None, None]
    jk = jvec + np.arange(2)[None, :, None]
    c1, c2, c3, c4 = vecd1 * vecd2, veca * vecd2, veca * vecb, vecd1 * vecb
    B, D = d.field.bz[ik, jk], d.field.d2bz[ik, jk]
    if d.field.interpolation == 'hermite':
        # Bz and its gradient from one C1 surface; the z^2 term stays bilinear
        bz, gx, gy = hermite_bz(d.field, ivec, jvec, veca, vecb, d.bunch.r, d.bunch.theta)
        d.bunch.bz[:] = bz - vecz * (c1 * D[0, 0] - c2 * D[1, 0] + c3 * D[1, 1] - c4 * D[0, 1]) * hrhtetm1
        d.bunch.bx[:] = gx * d.bunch.z
        d.bunch.by[:] = gy * d.bunch.z
    else:
        d.bunch.bz[:] = (c1 * (B[0, 0] - vecz * D[0, 0]) - c2 * (B[1, 0] - vecz * D[1, 0]) +
                         c3 * (B[1, 1] - vecz * D[1, 1]) - c4 * (B[0, 1] - vecz * D[0, 1])) * hrhtetm1
        for out, g in ((d.bunch.bx, d.field.dbzdx[ik, jk]), (d.bunch.by, d.field.dbzdy[ik, jk])):
            out[:] = (c1 * g[0, 0] - c2 * g[1, 0] + c3 * g[1, 1] - c4 * g[0, 1]) * d.bunch.z * hrhtetm1

    mask_bz0 = (d.bunch.in_box & (d.bunch.bz == 0.0)) | (~d.bunch.in_box)
    d.bunch.bx[mask_bz0] = 0.0
    d.bunch.by[mask_bz0] = 0.0
    d.bunch.bz[mask_bz0] = -0.00001

    mask_rle = d.bunch.r <= d.field.rmin
    d.bunch.bx[mask_rle] = 0.0
    d.bunch.by[mask_rle] = 0.0
    d.bunch.bz[mask_rle] = d.field.bz[0, 0]

    # E in the central region
    d.bunch.ex[:] = 0.0
    d.bunch.ey[:] = 0.0
    d.bunch.ez[:] = 0.0
    check_central(d.bunch, d.central)

    if d.bunch.n_central > 0:
        m_c = d.bunch.in_central
        ivec[m_c] = ((d.bunch.x[m_c] - d.central.xmin) * hx1m1).astype(np.int32)
        jvec[m_c] = ((d.bunch.y[m_c] - d.central.ymin) * hy1m1).astype(np.int32)
        veca[m_c] = d.bunch.x[m_c] - d.central.xmin - ivec[m_c] * d.central.dx
        vecb[m_c] = d.bunch.y[m_c] - d.central.ymin - jvec[m_c] * d.central.dy
        vecd1[m_c] = veca[m_c] - d.central.dx
        vecd2[m_c] = vecb[m_c] - d.central.dy
        vecp1[m_c] = vecd1[m_c] * vecd2[m_c]
        vecp2[m_c] = -veca[m_c] * vecd2[m_c]
        vecp3[m_c] = veca[m_c] * vecb[m_c]
        vecp4[m_c] = -vecd1[m_c] * vecb[m_c]
        vecz[m_c] = d.bunch.z[m_c]**2 * 0.5

        # 4x4 stencil around (i, j): P[a, b] = pot[i + a - 1, j + b - 1]
        ik = ivec[m_c] + np.arange(-1, 3)[:, None, None]
        jk = jvec[m_c] + np.arange(-1, 3)[None, :, None]
        P, D = d.central.pot[ik, jk], d.central.d2pot[ik, jk]
        vz, dx, dy = vecz[m_c], d.central.dx, d.central.dy
        vecb1[m_c], vecb2[m_c], vecb3[m_c], vecb4[m_c] = D[1, 1], D[2, 1], D[2, 2], D[1, 2]
        # Centred differences at the four corners, in the same operation order as the Fortran
        vex1[m_c] = (-(P[2, 1] - P[0, 1]) + vz * (D[2, 1] - D[0, 1])) * 0.5 / dx
        vex2[m_c] = (-(P[3, 1] - P[1, 1]) + vz * (D[3, 1] - D[1, 1])) * 0.5 / dx
        vex3[m_c] = (-(P[3, 2] - P[1, 2]) + vz * (D[3, 2] - D[1, 2])) * 0.5 / dx
        vex4[m_c] = (-(P[2, 2] - P[0, 2]) + vz * (D[2, 2] - D[0, 2])) * 0.5 / dx
        vey1[m_c] = (-(P[1, 2] - P[1, 0]) + vz * (D[1, 2] - D[1, 0])) * 0.5 / dy
        vey2[m_c] = (-(P[2, 2] - P[2, 0]) + vz * (D[2, 2] - D[2, 0])) * 0.5 / dy
        vey3[m_c] = (-(P[2, 3] - P[2, 1]) + vz * (D[2, 3] - D[2, 1])) * 0.5 / dy
        vey4[m_c] = (-(P[1, 3] - P[1, 1]) + vz * (D[1, 3] - D[1, 1])) * 0.5 / dy

        n_periods = int(d.rf.omega * d.run.t / (2.0 * np.pi))
        rf_phase = d.rf.omega * d.run.t - n_periods * 2.0 * np.pi
        pcright = d.rf.v_right * np.cos(rf_phase + d.rf.phase_right)
        pcleft = d.rf.v_left * np.cos(rf_phase + d.rf.phase_left)
        vecpc[d.bunch.x >= 0.0] = pcright
        vecpc[d.bunch.x < 0.0] = pcleft

        m_c = d.bunch.in_central
        d.bunch.ex[m_c] = vecpc[m_c] * (
            vex1[m_c] * vecp1[m_c] + vex2[m_c] * vecp2[m_c] +
            vex3[m_c] * vecp3[m_c] + vex4[m_c] * vecp4[m_c]
        ) * hx1hy1m1
        d.bunch.ey[m_c] = vecpc[m_c] * (
            vey1[m_c] * vecp1[m_c] + vey2[m_c] * vecp2[m_c] +
            vey3[m_c] * vecp3[m_c] + vey4[m_c] * vecp4[m_c]
        ) * hx1hy1m1
        d.bunch.ez[m_c] = d.bunch.z[m_c] * vecpc[m_c] * (
            vecp1[m_c] * vecb1[m_c] + vecp2[m_c] * vecb2[m_c] +
            vecp3[m_c] * vecb3[m_c] + vecp4[m_c] * vecb4[m_c]
        ) * hx1hy1m1

    # E in the outer region
    if (d.bunch.n_in_box == d.bunch.n_central) and (len(d.extraction.dipoles) == 0) and (len(d.extraction.deflectors) == 0) and (len(d.extraction.grids) == 1):
        pass  # Fortran goto 300
    elif d.bunch.n_in_box == d.bunch.n_central:
        add_vertical_sorter_field(d.bunch, d.sorter, d.rf, d.run)
        add_extraction_fields(d.bunch, d.extraction, d.outputs, d.run)
        return
    else:
        cos_dee = np.cos(d.rf.dee_half_angle)
        sin_dee = np.sin(d.rf.dee_half_angle)
        mask_nc = ~d.bunch.in_central
        vecx[mask_nc] = np.abs(d.bunch.x[mask_nc]) * cos_dee + (np.abs(d.bunch.y[mask_nc]) - d.rf.gap_offset) * sin_dee
        vecy[mask_nc] = -np.abs(d.bunch.x[mask_nc]) * sin_dee + (np.abs(d.bunch.y[mask_nc]) - d.rf.gap_offset) * cos_dee

        check_outer(d.bunch, d.outer, vecx, vecy)
        n_periods = int(d.rf.omega * d.run.t / (2.0 * np.pi))
        rf_phase = d.rf.omega * d.run.t - n_periods * 2.0 * np.pi

        if d.bunch.n_outer > 0:
            if d.bunch.n_central == 0:
                pcright = d.rf.v_right * np.cos(rf_phase + d.rf.phase_right)
                pcleft = d.rf.v_left * np.cos(rf_phase + d.rf.phase_left)
                vecpc[d.bunch.x >= 0.0] = pcright
                vecpc[d.bunch.x < 0.0] = pcleft

            m_p = d.bunch.in_outer
            ivec[m_p] = ((vecx[m_p] - d.outer.xmin) * hx2m1).astype(np.int32)
            jvec[m_p] = ((vecy[m_p] - d.outer.ymin) * hy2m1).astype(np.int32)
            veca[m_p] = vecx[m_p] - d.outer.xmin - ivec[m_p] * d.outer.dx
            vecb[m_p] = vecy[m_p] - d.outer.ymin - jvec[m_p] * d.outer.dy
            vecd1[m_p] = veca[m_p] - d.outer.dx
            vecd2[m_p] = vecb[m_p] - d.outer.dy
            vecp1[m_p] = vecd1[m_p] * vecd2[m_p]
            vecp2[m_p] = -veca[m_p] * vecd2[m_p]
            vecp3[m_p] = veca[m_p] * vecb[m_p]
            vecp4[m_p] = -vecd1[m_p] * vecb[m_p]
            vecz[m_p] = d.bunch.z[m_p]**2 * 0.5

            ik = ivec[m_p] + np.arange(2)[:, None, None]
            jk = jvec[m_p] + np.arange(2)[None, :, None]
            P, D = d.outer.pot[ik, jk], d.outer.d2pot[ik, jk]
            vec1[m_p], vec2[m_p], vec3[m_p], vec4[m_p] = P[0, 0], P[1, 0], P[1, 1], P[0, 1]
            vecb1[m_p], vecb2[m_p], vecb3[m_p], vecb4[m_p] = D[0, 0], D[1, 0], D[1, 1], D[0, 1]

            d.bunch.ex[m_p] = -vecpc[m_p] * (
                vecd2[m_p] * (vec1[m_p] - vec2[m_p] -
                                vecz[m_p] * (vecb1[m_p] - vecb2[m_p])) +
                vecb[m_p] * (vec3[m_p] - vec4[m_p] +
                                   vecz[m_p] * (vecb3[m_p] - vecb4[m_p]))
            ) * hx2hy2m1
            d.bunch.ey[m_p] = -vecpc[m_p] * (
                vecd1[m_p] * (vec1[m_p] - vec4[m_p] -
                                vecz[m_p] * (vecb1[m_p] - vecb4[m_p])) -
                veca[m_p] * (vec2[m_p] - vec3[m_p] +
                               vecz[m_p] * (vecb2[m_p] - vecb3[m_p]))
            ) * hx2hy2m1

            # Back from the gap frame
            vec_tmp = d.bunch.ex[m_p] * cos_dee - d.bunch.ey[m_p] * sin_dee
            d.bunch.ey[m_p] = d.bunch.ex[m_p] * sin_dee + d.bunch.ey[m_p] * cos_dee
            d.bunch.ex[m_p] = vec_tmp

            d.bunch.ez[m_p] = d.bunch.z[m_p] * vecpc[m_p] * (
                vecp1[m_p] * vecb1[m_p] + vecp2[m_p] * vecb2[m_p] +
                vecp3[m_p] * vecb3[m_p] + vecp4[m_p] * vecb4[m_p]
            ) * hx2hy2m1

            # Symmetries of the other quadrants
            m1 = m_p & (d.bunch.x < 0.0) & (d.bunch.y > 0.0)
            d.bunch.ex[m1] = d.rf.harmonic_sign * d.bunch.ex[m1]
            d.bunch.ey[m1] = -d.rf.harmonic_sign * d.bunch.ey[m1]

            m2 = m_p & (d.bunch.x > 0.0) & (d.bunch.y < 0.0)
            d.bunch.ey[m2] = -d.bunch.ey[m2]

            m3 = m_p & (d.bunch.x < 0.0) & (d.bunch.y < 0.0)
            d.bunch.ex[m3] = d.rf.harmonic_sign * d.bunch.ex[m3]
            d.bunch.ey[m3] = d.rf.harmonic_sign * d.bunch.ey[m3]

            m4 = m_p & (d.bunch.x < 0.0)
            d.bunch.ez[m4] = -d.rf.harmonic_sign * d.bunch.ez[m4]

    add_vertical_sorter_field(d.bunch, d.sorter, d.rf, d.run)
    add_extraction_fields(d.bunch, d.extraction, d.outputs, d.run)


def add_vertical_sorter_field(bunch, sorter, rf, run):
    """Add the vertical sorter Ez; sorter.voltage_seen is its voltage at particle 0."""
    sorter.voltage_seen = 0.0
    if sorter.gap == 0.0:
        return
    mask_ver = ((bunch.r >= sorter.rmin) & (bunch.r <= sorter.rmax) &
                (bunch.theta >= sorter.theta1) & (bunch.theta <= sorter.theta2))
    bunch.ez[mask_ver] = (bunch.ez[mask_ver] -
        sorter.voltage * (1.0 - sorter.dv) / sorter.gap / 2.0 *
        np.cos(2.0 * rf.omega * run.t * (1.0 + sorter.nuz / 2.0 / rf.harmonic) +
               2.0 * sorter.phase - np.pi / 2.0) -
        sorter.voltage * (1.0 + sorter.dv) / sorter.gap / 2.0 *
        np.cos(2.0 * rf.omega * run.t * (1.0 - sorter.nuz / 2.0 / rf.harmonic) +
               2.0 * sorter.phase - np.pi / 2.0))
    if (bunch.r[0] >= sorter.rmin and bunch.r[0] <= sorter.rmax and
            bunch.theta[0] >= sorter.theta1 and bunch.theta[0] <= sorter.theta2):
        sorter.voltage_seen = -bunch.ez[0] * sorter.gap


def add_extraction_fields(bunch, extraction, outputs, run):
    """Deflector E, dipole B and extraction grid crossings."""
    for i in range(len(extraction.deflectors)):
        check_deflector(bunch, extraction, run, i)
        if bunch.n_in_element != 0:
            mask = bunch.in_element
            vec_p = np.sqrt(bunch.px[mask]**2 + bunch.py[mask]**2)
            bunch.ex[mask] = extraction.deflectors[i].e * bunch.py[mask] / vec_p
            bunch.ey[mask] = -extraction.deflectors[i].e * bunch.px[mask] / vec_p

    for i in range(len(extraction.dipoles)):
        check_dipole(bunch, extraction, i)
        if bunch.n_in_element != 0:
            mask = bunch.in_element
            vec_dm = np.sqrt((bunch.x[mask] - extraction.dipoles[i].x)**2 + (bunch.y[mask] - extraction.dipoles[i].y)**2)
            vec1_dm = vec_dm - extraction.dipoles[i].r
            bunch.bx[mask] += bunch.z[mask] * extraction.dipoles[i].gradient * (bunch.x[mask] - extraction.dipoles[i].x) / vec_dm
            bunch.by[mask] += bunch.z[mask] * extraction.dipoles[i].gradient * (bunch.y[mask] - extraction.dipoles[i].y) / vec_dm
            bunch.bz[mask] += vec1_dm * extraction.dipoles[i].gradient

    if extraction.mode < 0.0 and len(extraction.grids) > 0:
        if run.step >= extraction.grid_start_step:
            for kg in range(len(extraction.grids)):
                check_grid(bunch, extraction, run, kg)
                bunch.grid_reached[bunch.crossed_grid] = kg + 1  # 1-based in data
                for i in range(bunch.n_particles):
                    if bunch.crossed_grid[i]:
                        if f'grid{kg}' in outputs.open:
                            outputs.open[f'grid{kg}'].write(
                                f' {i+1:5d} {run.step:6d}'
                                f' {bunch.x[i]:15.8e} {bunch.y[i]:15.8e} {bunch.z[i]:15.8e}'
                                f' {bunch.px[i]:15.8e} {bunch.py[i]:15.8e} {bunch.pz[i]:15.8e}'
                                f' {bunch.bz[i]:8.5f}\n')


def compute_bfield(bunch, field):
    """Magnetic field only, used while building the initial bunch."""
    hrhtet = field.dr * field.dtheta
    ivec = np.zeros(bunch.n_particles, dtype=np.int32)
    jvec = np.zeros(bunch.n_particles, dtype=np.int32)
    vec1 = np.zeros(bunch.n_particles)
    vec2 = np.zeros(bunch.n_particles)
    vec3 = np.zeros(bunch.n_particles)
    vec4 = np.zeros(bunch.n_particles)
    veca = np.zeros(bunch.n_particles)
    vecb1 = np.zeros(bunch.n_particles)
    vecb2 = np.zeros(bunch.n_particles)
    vecb3 = np.zeros(bunch.n_particles)
    vecb4 = np.zeros(bunch.n_particles)
    vecb = np.zeros(bunch.n_particles)
    vecd1 = np.zeros(bunch.n_particles)
    vecd2 = np.zeros(bunch.n_particles)
    vecp1 = np.zeros(bunch.n_particles)
    vecp2 = np.zeros(bunch.n_particles)
    vecp3 = np.zeros(bunch.n_particles)
    vecp4 = np.zeros(bunch.n_particles)
    vecx1 = np.zeros(bunch.n_particles)
    vecx2 = np.zeros(bunch.n_particles)
    vecx3 = np.zeros(bunch.n_particles)
    vecx4 = np.zeros(bunch.n_particles)
    vecy1 = np.zeros(bunch.n_particles)
    vecy2 = np.zeros(bunch.n_particles)
    vecy3 = np.zeros(bunch.n_particles)
    vecy4 = np.zeros(bunch.n_particles)
    vecz = np.zeros(bunch.n_particles)
    bunch.theta = np.arccos(np.clip(bunch.x / bunch.r, -1.0, 1.0))
    bunch.theta[bunch.y < 0.0] = 2.0 * np.pi - bunch.theta[bunch.y < 0.0]
    bunch.theta[bunch.theta >= 2.0 * np.pi] -= 2.0 * np.pi

    mask = bunch.in_box
    ivec[mask] = ((bunch.r[mask] - field.rmin) / field.dr).astype(np.int32)
    jvec[mask] = (bunch.theta[mask] / field.dtheta).astype(np.int32)
    veca[mask] = bunch.r[mask] - field.rmin - ivec[mask] * field.dr
    vecb[mask] = bunch.theta[mask] - jvec[mask] * field.dtheta
    vecd1[mask] = veca[mask] - field.dr
    vecd2[mask] = vecb[mask] - field.dtheta
    vecp1[mask] = vecd1[mask] * vecd2[mask]
    vecp2[mask] = -veca[mask] * vecd2[mask]
    vecp3[mask] = veca[mask] * vecb[mask]
    vecp4[mask] = -vecd1[mask] * vecb[mask]
    vecz[mask] = bunch.z[mask]**2 / 2.0

    ik = ivec[mask] + np.arange(2)[:, None, None]
    jk = jvec[mask] + np.arange(2)[None, :, None]
    B, D = field.bz[ik, jk], -vecz[mask] * field.d2bz[ik, jk]
    X, Y = field.dbzdx[ik, jk], field.dbzdy[ik, jk]
    vecb1[mask], vecb2[mask], vecb3[mask], vecb4[mask] = B[0, 0], B[1, 0], B[1, 1], B[0, 1]
    vec1[mask], vec2[mask], vec3[mask], vec4[mask] = D[0, 0], D[1, 0], D[1, 1], D[0, 1]
    vecx1[mask], vecx2[mask], vecx3[mask], vecx4[mask] = X[0, 0], X[1, 0], X[1, 1], X[0, 1]
    vecy1[mask], vecy2[mask], vecy3[mask], vecy4[mask] = Y[0, 0], Y[1, 0], Y[1, 1], Y[0, 1]

    if field.interpolation == 'hermite':
        bz, gx, gy = hermite_bz(field, ivec[mask], jvec[mask], veca[mask], vecb[mask], bunch.r[mask],
                                bunch.theta[mask])
        bunch.bz[mask] = bz + (vecp1[mask] * vec1[mask] + vecp2[mask] * vec2[mask] +
                               vecp3[mask] * vec3[mask] + vecp4[mask] * vec4[mask]) / hrhtet
        bunch.bx[mask] = bunch.z[mask] * gx
        bunch.by[mask] = bunch.z[mask] * gy
    else:
        vec1[mask] += vecb1[mask]
        vec2[mask] += vecb2[mask]
        vec3[mask] += vecb3[mask]
        vec4[mask] += vecb4[mask]
        bunch.bz[mask] = (vecp1[mask] * vec1[mask] + vecp2[mask] * vec2[mask] +
                      vecp3[mask] * vec3[mask] + vecp4[mask] * vec4[mask]) / hrhtet
        bunch.bx[mask] = bunch.z[mask] * (
            vecp1[mask] * vecx1[mask] + vecp2[mask] * vecx2[mask] +
            vecp3[mask] * vecx3[mask] + vecp4[mask] * vecx4[mask]) / hrhtet
        bunch.by[mask] = bunch.z[mask] * (
            vecp1[mask] * vecy1[mask] + vecp2[mask] * vecy2[mask] +
            vecp3[mask] * vecy3[mask] + vecp4[mask] * vecy4[mask]) / hrhtet

    mask_rle = bunch.r <= field.rmin
    bunch.bx[mask_rle] = 0.0
    bunch.by[mask_rle] = 0.0
    bunch.bz[mask_rle] = field.bz[0, 0]
