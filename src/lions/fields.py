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

    for k in range(n):
        ik = ivec[k]
        jk = jvec[k]
        d.bunch.bz[k] = (vecd1[k] * vecd2[k] * (d.field.bz[ik, jk] -
                   vecz[k] * d.field.d2bz[ik, jk]) -
                   veca[k] * vecd2[k] * (d.field.bz[ik+1, jk] -
                   vecz[k] * d.field.d2bz[ik+1, jk]) +
                   veca[k] * vecb[k] * (d.field.bz[ik+1, jk+1] -
                   vecz[k] * d.field.d2bz[ik+1, jk+1]) -
                   vecd1[k] * vecb[k] * (d.field.bz[ik, jk+1] -
                   vecz[k] * d.field.d2bz[ik, jk+1])) * hrhtetm1
        d.bunch.bx[k] = (vecd1[k] * vecd2[k] * d.field.dbzdx[ik, jk] -
                   veca[k] * vecd2[k] * d.field.dbzdx[ik+1, jk] +
                   veca[k] * vecb[k] * d.field.dbzdx[ik+1, jk+1] -
                   vecd1[k] * vecb[k] * d.field.dbzdx[ik, jk+1]
                   ) * d.bunch.z[k] * hrhtetm1
        d.bunch.by[k] = (vecd1[k] * vecd2[k] * d.field.dbzdy[ik, jk] -
                   veca[k] * vecd2[k] * d.field.dbzdy[ik+1, jk] +
                   veca[k] * vecb[k] * d.field.dbzdy[ik+1, jk+1] -
                   vecd1[k] * vecb[k] * d.field.dbzdy[ik, jk+1]
                   ) * d.bunch.z[k] * hrhtetm1

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

        for k in range(n):
            if d.bunch.in_central[k]:
                ik = ivec[k]
                jk = jvec[k]
                vecb1[k] = d.central.d2pot[ik, jk]
                vecb2[k] = d.central.d2pot[ik+1, jk]
                vecb3[k] = d.central.d2pot[ik+1, jk+1]
                vecb4[k] = d.central.d2pot[ik, jk+1]
                vz = vecz[k]
                vex1[k] = (-(d.central.pot[ik+1, jk] - d.central.pot[ik-1, jk]) +
                             vz * (d.central.d2pot[ik+1, jk] - d.central.d2pot[ik-1, jk])
                             ) * 0.5 / d.central.dx
                vex2[k] = (-(d.central.pot[ik+2, jk] - d.central.pot[ik, jk]) +
                             vz * (d.central.d2pot[ik+2, jk] - d.central.d2pot[ik, jk])
                             ) * 0.5 / d.central.dx
                vex3[k] = (-(d.central.pot[ik+2, jk+1] - d.central.pot[ik, jk+1]) +
                             vz * (d.central.d2pot[ik+2, jk+1] - d.central.d2pot[ik, jk+1])
                             ) * 0.5 / d.central.dx
                vex4[k] = (-(d.central.pot[ik+1, jk+1] - d.central.pot[ik-1, jk+1]) +
                             vz * (d.central.d2pot[ik+1, jk+1] - d.central.d2pot[ik-1, jk+1])
                             ) * 0.5 / d.central.dx
                vey1[k] = (-(d.central.pot[ik, jk+1] - d.central.pot[ik, jk-1]) +
                             vz * (d.central.d2pot[ik, jk+1] - d.central.d2pot[ik, jk-1])
                             ) * 0.5 / d.central.dy
                vey2[k] = (-(d.central.pot[ik+1, jk+1] - d.central.pot[ik+1, jk-1]) +
                             vz * (d.central.d2pot[ik+1, jk+1] - d.central.d2pot[ik+1, jk-1])
                             ) * 0.5 / d.central.dy
                vey3[k] = (-(d.central.pot[ik+1, jk+2] - d.central.pot[ik+1, jk]) +
                             vz * (d.central.d2pot[ik+1, jk+2] - d.central.d2pot[ik+1, jk])
                             ) * 0.5 / d.central.dy
                vey4[k] = (-(d.central.pot[ik, jk+2] - d.central.pot[ik, jk]) +
                             vz * (d.central.d2pot[ik, jk+2] - d.central.d2pot[ik, jk])
                             ) * 0.5 / d.central.dy

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

            for k in range(n):
                if d.bunch.in_outer[k]:
                    ik = ivec[k]
                    jk = jvec[k]
                    vec1[k] = d.outer.pot[ik, jk]
                    vec2[k] = d.outer.pot[ik+1, jk]
                    vec3[k] = d.outer.pot[ik+1, jk+1]
                    vec4[k] = d.outer.pot[ik, jk+1]
                    vecb1[k] = d.outer.d2pot[ik, jk]
                    vecb2[k] = d.outer.d2pot[ik+1, jk]
                    vecb3[k] = d.outer.d2pot[ik+1, jk+1]
                    vecb4[k] = d.outer.d2pot[ik, jk+1]

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

    for k in range(bunch.n_particles):
        if bunch.in_box[k]:
            ik = ivec[k]
            jk = jvec[k]
            vk = vecz[k]
            vecb1[k] = field.bz[ik, jk]
            vecb2[k] = field.bz[ik+1, jk]
            vecb3[k] = field.bz[ik+1, jk+1]
            vecb4[k] = field.bz[ik, jk+1]
            vec1[k] = -vk * field.d2bz[ik, jk]
            vec2[k] = -vk * field.d2bz[ik+1, jk]
            vec3[k] = -vk * field.d2bz[ik+1, jk+1]
            vec4[k] = -vk * field.d2bz[ik, jk+1]
            vecx1[k] = field.dbzdx[ik, jk]
            vecx2[k] = field.dbzdx[ik+1, jk]
            vecx3[k] = field.dbzdx[ik+1, jk+1]
            vecx4[k] = field.dbzdx[ik, jk+1]
            vecy1[k] = field.dbzdy[ik, jk]
            vecy2[k] = field.dbzdy[ik+1, jk]
            vecy3[k] = field.dbzdy[ik+1, jk+1]
            vecy4[k] = field.dbzdy[ik, jk+1]

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
