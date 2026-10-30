"""Magnetic and electric fields at the particles."""


import numpy as np

from .apertures import check_central, check_deflector, check_dipole, check_grid, check_outer


def add_harmonics(rmin, rmax, nr, ntet, tbz, d2bz, dbzx, dbzy, ncoh, coh, phihar):
    """Add first and second field harmonics to the map arrays, in place."""
    pi = np.pi
    dr = (rmax - rmin) / (nr - 1.0)
    dtheta = 360.0 / (ntet - 1.0) * pi / 180.0
    phase1 = phihar[0] * pi / 180.0
    phase2 = phihar[1] * pi / 180.0

    if ncoh[0] != 0:
        for i in range(nr):
            r = rmin + i * dr
            for j in range(ntet - 1):
                theta = j * dtheta
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                r2 = x**2 + y**2
                co1 = np.cos(-phase1)
                si1 = np.sin(-phase1)
                s = x * co1 - y * si1
                pb = coh[0, 0] * s
                px_val = coh[0, 0] * co1
                py_val = -coh[0, 0] * si1
                p2 = 0.0
                for k in range(1, ncoh[0]):
                    k1 = k
                    k2 = k - 1
                    r2pk1 = r2**k1
                    if k == 1 and r2 <= 0.0:
                        r2pk2 = 1.0
                    else:
                        r2pk2 = r2**k2
                    pb += coh[k, 0] * r2pk1 * s
                    px_val += coh[k, 0] * r2pk2 * (2.0 * x * k1 * s + r2 * co1)
                    py_val += coh[k, 0] * r2pk2 * (2.0 * y * k1 * s - r2 * si1)
                    p2 += coh[k, 0] * r2pk2 * 4.0 * s * k1 * (k + 1)
                tbz[i, j] += pb
                dbzx[i, j] += px_val
                dbzy[i, j] += py_val
                d2bz[i, j] += p2

    if ncoh[1] != 0:
        for i in range(nr):
            r = rmin + i * dr
            for j in range(ntet - 1):
                theta = j * dtheta
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                r2 = x**2 + y**2
                co2 = np.cos(-phase2)
                si2 = np.sin(-phase2)
                s = (x**2 - y**2) * co2 - 2.0 * x * y * si2
                s2 = 2.0 * x * co2 - 2.0 * y * si2
                s3 = -2.0 * y * co2 - 2.0 * x * si2
                pb = coh[0, 1] * s
                px_val = coh[0, 1] * s2
                py_val = -coh[0, 1] * s3
                p2 = 0.0
                for k in range(1, ncoh[1]):
                    k1 = k
                    k2 = k - 1
                    r2pk1 = r2**k1
                    if k == 1 and r2 <= 0.0:
                        r2pk2 = 1.0
                    else:
                        r2pk2 = r2**k2
                    pb += coh[k, 1] * r2pk1 * s
                    px_val += coh[k, 1] * r2pk2 * (2.0 * x * k1 * s + r2 * s2)
                    py_val += coh[k, 1] * r2pk2 * (2.0 * y * k1 * s + r2 * s3)
                    p2 += coh[k, 1] * r2pk2 * 4.0 * s * k1 * (k + 2)
                tbz[i, j] += pb
                dbzx[i, j] += px_val
                dbzy[i, j] += py_val
                d2bz[i, j] += p2

    # Periodicity in theta
    tbz[:, ntet - 1] = tbz[:, 0]
    dbzx[:, ntet - 1] = dbzx[:, 0]
    dbzy[:, ntet - 1] = dbzy[:, 0]
    d2bz[:, ntet - 1] = d2bz[:, 0]


def compute_fields(d):
    """Magnetic and electric fields at all particles for the current time."""
    hrhtetm1 = 1.0 / (d.map_dr * d.map_dtheta)
    hrm1 = 1.0 / d.map_dr
    htetm1 = 1.0 / d.map_dtheta
    hx1hy1m1 = 1.0 / (d.central_dx * d.central_dy)
    hx1m1 = 1.0 / d.central_dx
    hx2hy2m1 = 1.0 / (d.outer_dx * d.outer_dy)
    hx2m1 = 1.0 / d.outer_dx
    hy1m1 = 1.0 / d.central_dy
    hy2m1 = 1.0 / d.outer_dy
    ivec = np.zeros(d.n_particles, dtype=np.int32)
    jvec = np.zeros(d.n_particles, dtype=np.int32)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    veca = np.zeros(d.n_particles)
    vecb1 = np.zeros(d.n_particles)
    vecb2 = np.zeros(d.n_particles)
    vecb3 = np.zeros(d.n_particles)
    vecb4 = np.zeros(d.n_particles)
    vecb = np.zeros(d.n_particles)
    vecd1 = np.zeros(d.n_particles)
    vecd2 = np.zeros(d.n_particles)
    vecp1 = np.zeros(d.n_particles)
    vecp2 = np.zeros(d.n_particles)
    vecp3 = np.zeros(d.n_particles)
    vecp4 = np.zeros(d.n_particles)
    vecpc = np.zeros(d.n_particles)
    vecx = np.zeros(d.n_particles)
    vecy = np.zeros(d.n_particles)
    vecz = np.zeros(d.n_particles)
    vex1 = np.zeros(d.n_particles)
    vex2 = np.zeros(d.n_particles)
    vex3 = np.zeros(d.n_particles)
    vex4 = np.zeros(d.n_particles)
    vey1 = np.zeros(d.n_particles)
    vey2 = np.zeros(d.n_particles)
    vey3 = np.zeros(d.n_particles)
    vey4 = np.zeros(d.n_particles)
    n = d.n_particles

    # B from the median-plane map, first order in z
    n_bad = np.count_nonzero((d.x / d.r > 1.0) | (d.x / d.r < -1.0))
    if n_bad != 0:
        print(f'particles with |x/r| > 1: {n_bad}')
        for i in range(n):
            ratio = d.x[i] / d.r[i]
            if ratio > 1.0 or ratio < -1.0:
                print(f'{i} {d.x[i]} {d.y[i]} {d.r[i]} {ratio}')

    veca = d.x / d.r
    mask_normal = (veca > -1.0) & (veca < 1.0)
    mask_ge1 = veca >= 1.0
    mask_lem1 = veca <= -1.0
    d.theta[mask_normal] = np.arccos(veca[mask_normal])
    d.theta[mask_ge1] = np.arccos(1.0)
    d.theta[mask_lem1] = np.arccos(-1.0)
    mask_yneg = d.y < 0.0
    d.theta[mask_yneg] = 2.0 * np.pi - d.theta[mask_yneg]

    ivec[:] = 0
    jvec[:] = 0
    mask = d.in_box
    ivec[mask] = ((d.r[mask] - d.map_rmin) * hrm1).astype(np.int32)
    jvec[mask] = (d.theta[mask] * htetm1).astype(np.int32)
    veca[:] = d.r - d.map_rmin - ivec * d.map_dr
    vecb[:] = d.theta - jvec * d.map_dtheta
    vecd1[:] = veca - d.map_dr
    vecd2[:] = vecb - d.map_dtheta
    vecz[:] = d.z * d.z * 0.5

    for k in range(n):
        ik = ivec[k]
        jk = jvec[k]
        d.bz[k] = (vecd1[k] * vecd2[k] * (d.bz_map[ik, jk] -
                   vecz[k] * d.d2bz_map[ik, jk]) -
                   veca[k] * vecd2[k] * (d.bz_map[ik+1, jk] -
                   vecz[k] * d.d2bz_map[ik+1, jk]) +
                   veca[k] * vecb[k] * (d.bz_map[ik+1, jk+1] -
                   vecz[k] * d.d2bz_map[ik+1, jk+1]) -
                   vecd1[k] * vecb[k] * (d.bz_map[ik, jk+1] -
                   vecz[k] * d.d2bz_map[ik, jk+1])) * hrhtetm1
        d.bx[k] = (vecd1[k] * vecd2[k] * d.dbzdx_map[ik, jk] -
                   veca[k] * vecd2[k] * d.dbzdx_map[ik+1, jk] +
                   veca[k] * vecb[k] * d.dbzdx_map[ik+1, jk+1] -
                   vecd1[k] * vecb[k] * d.dbzdx_map[ik, jk+1]
                   ) * d.z[k] * hrhtetm1
        d.by[k] = (vecd1[k] * vecd2[k] * d.dbzdy_map[ik, jk] -
                   veca[k] * vecd2[k] * d.dbzdy_map[ik+1, jk] +
                   veca[k] * vecb[k] * d.dbzdy_map[ik+1, jk+1] -
                   vecd1[k] * vecb[k] * d.dbzdy_map[ik, jk+1]
                   ) * d.z[k] * hrhtetm1

    mask_bz0 = (d.in_box & (d.bz == 0.0)) | (~d.in_box)
    d.bx[mask_bz0] = 0.0
    d.by[mask_bz0] = 0.0
    d.bz[mask_bz0] = -0.00001

    mask_rle = d.r <= d.map_rmin
    d.bx[mask_rle] = 0.0
    d.by[mask_rle] = 0.0
    d.bz[mask_rle] = d.bz_map[0, 0]

    # E in the central region
    d.ex[:] = 0.0
    d.ey[:] = 0.0
    d.ez[:] = 0.0
    check_central(d)

    if d.n_central > 0:
        m_c = d.in_central
        ivec[m_c] = ((d.x[m_c] - d.central_xmin) * hx1m1).astype(np.int32)
        jvec[m_c] = ((d.y[m_c] - d.central_ymin) * hy1m1).astype(np.int32)
        veca[m_c] = d.x[m_c] - d.central_xmin - ivec[m_c] * d.central_dx
        vecb[m_c] = d.y[m_c] - d.central_ymin - jvec[m_c] * d.central_dy
        vecd1[m_c] = veca[m_c] - d.central_dx
        vecd2[m_c] = vecb[m_c] - d.central_dy
        vecp1[m_c] = vecd1[m_c] * vecd2[m_c]
        vecp2[m_c] = -veca[m_c] * vecd2[m_c]
        vecp3[m_c] = veca[m_c] * vecb[m_c]
        vecp4[m_c] = -vecd1[m_c] * vecb[m_c]
        vecz[m_c] = d.z[m_c]**2 * 0.5

        for k in range(n):
            if d.in_central[k]:
                ik = ivec[k]
                jk = jvec[k]
                vecb1[k] = d.central_d2pot[ik, jk]
                vecb2[k] = d.central_d2pot[ik+1, jk]
                vecb3[k] = d.central_d2pot[ik+1, jk+1]
                vecb4[k] = d.central_d2pot[ik, jk+1]
                vz = vecz[k]
                vex1[k] = (-(d.central_pot[ik+1, jk] - d.central_pot[ik-1, jk]) +
                             vz * (d.central_d2pot[ik+1, jk] - d.central_d2pot[ik-1, jk])
                             ) * 0.5 / d.central_dx
                vex2[k] = (-(d.central_pot[ik+2, jk] - d.central_pot[ik, jk]) +
                             vz * (d.central_d2pot[ik+2, jk] - d.central_d2pot[ik, jk])
                             ) * 0.5 / d.central_dx
                vex3[k] = (-(d.central_pot[ik+2, jk+1] - d.central_pot[ik, jk+1]) +
                             vz * (d.central_d2pot[ik+2, jk+1] - d.central_d2pot[ik, jk+1])
                             ) * 0.5 / d.central_dx
                vex4[k] = (-(d.central_pot[ik+1, jk+1] - d.central_pot[ik-1, jk+1]) +
                             vz * (d.central_d2pot[ik+1, jk+1] - d.central_d2pot[ik-1, jk+1])
                             ) * 0.5 / d.central_dx
                vey1[k] = (-(d.central_pot[ik, jk+1] - d.central_pot[ik, jk-1]) +
                             vz * (d.central_d2pot[ik, jk+1] - d.central_d2pot[ik, jk-1])
                             ) * 0.5 / d.central_dy
                vey2[k] = (-(d.central_pot[ik+1, jk+1] - d.central_pot[ik+1, jk-1]) +
                             vz * (d.central_d2pot[ik+1, jk+1] - d.central_d2pot[ik+1, jk-1])
                             ) * 0.5 / d.central_dy
                vey3[k] = (-(d.central_pot[ik+1, jk+2] - d.central_pot[ik+1, jk]) +
                             vz * (d.central_d2pot[ik+1, jk+2] - d.central_d2pot[ik+1, jk])
                             ) * 0.5 / d.central_dy
                vey4[k] = (-(d.central_pot[ik, jk+2] - d.central_pot[ik, jk]) +
                             vz * (d.central_d2pot[ik, jk+2] - d.central_d2pot[ik, jk])
                             ) * 0.5 / d.central_dy

        n_periods = int(d.omega_rf * d.t / (2.0 * np.pi))
        rf_phase = d.omega_rf * d.t - n_periods * 2.0 * np.pi
        pcright = d.v_dee_right * np.cos(rf_phase + d.phase_right)
        pcleft = d.v_dee_left * np.cos(rf_phase + d.phase_left)
        vecpc[d.x >= 0.0] = pcright
        vecpc[d.x < 0.0] = pcleft

        m_c = d.in_central
        d.ex[m_c] = vecpc[m_c] * (
            vex1[m_c] * vecp1[m_c] + vex2[m_c] * vecp2[m_c] +
            vex3[m_c] * vecp3[m_c] + vex4[m_c] * vecp4[m_c]
        ) * hx1hy1m1
        d.ey[m_c] = vecpc[m_c] * (
            vey1[m_c] * vecp1[m_c] + vey2[m_c] * vecp2[m_c] +
            vey3[m_c] * vecp3[m_c] + vey4[m_c] * vecp4[m_c]
        ) * hx1hy1m1
        d.ez[m_c] = d.z[m_c] * vecpc[m_c] * (
            vecp1[m_c] * vecb1[m_c] + vecp2[m_c] * vecb2[m_c] +
            vecp3[m_c] * vecb3[m_c] + vecp4[m_c] * vecb4[m_c]
        ) * hx1hy1m1

    # E in the outer region
    if (d.n_in_box == d.n_central) and (d.n_dipoles == 0) and (d.n_deflectors == 0) and (d.n_grids == 1):
        pass  # Fortran goto 300
    elif d.n_in_box == d.n_central:
        add_vertical_sorter_field(d)
        add_extraction_fields(d)
        return
    else:
        cos_dee = np.cos(d.dee_half_angle)
        sin_dee = np.sin(d.dee_half_angle)
        mask_nc = ~d.in_central
        vecx[mask_nc] = np.abs(d.x[mask_nc]) * cos_dee + (np.abs(d.y[mask_nc]) - d.gap_offset) * sin_dee
        vecy[mask_nc] = -np.abs(d.x[mask_nc]) * sin_dee + (np.abs(d.y[mask_nc]) - d.gap_offset) * cos_dee

        check_outer(d, vecx, vecy)
        n_periods = int(d.omega_rf * d.t / (2.0 * np.pi))
        rf_phase = d.omega_rf * d.t - n_periods * 2.0 * np.pi

        if d.n_outer > 0:
            if d.n_central == 0:
                pcright = d.v_dee_right * np.cos(rf_phase + d.phase_right)
                pcleft = d.v_dee_left * np.cos(rf_phase + d.phase_left)
                vecpc[d.x >= 0.0] = pcright
                vecpc[d.x < 0.0] = pcleft

            m_p = d.in_outer
            ivec[m_p] = ((vecx[m_p] - d.outer_xmin) * hx2m1).astype(np.int32)
            jvec[m_p] = ((vecy[m_p] - d.outer_ymin) * hy2m1).astype(np.int32)
            veca[m_p] = vecx[m_p] - d.outer_xmin - ivec[m_p] * d.outer_dx
            vecb[m_p] = vecy[m_p] - d.outer_ymin - jvec[m_p] * d.outer_dy
            vecd1[m_p] = veca[m_p] - d.outer_dx
            vecd2[m_p] = vecb[m_p] - d.outer_dy
            vecp1[m_p] = vecd1[m_p] * vecd2[m_p]
            vecp2[m_p] = -veca[m_p] * vecd2[m_p]
            vecp3[m_p] = veca[m_p] * vecb[m_p]
            vecp4[m_p] = -vecd1[m_p] * vecb[m_p]
            vecz[m_p] = d.z[m_p]**2 * 0.5

            for k in range(n):
                if d.in_outer[k]:
                    ik = ivec[k]
                    jk = jvec[k]
                    vec1[k] = d.outer_pot[ik, jk]
                    vec2[k] = d.outer_pot[ik+1, jk]
                    vec3[k] = d.outer_pot[ik+1, jk+1]
                    vec4[k] = d.outer_pot[ik, jk+1]
                    vecb1[k] = d.outer_d2pot[ik, jk]
                    vecb2[k] = d.outer_d2pot[ik+1, jk]
                    vecb3[k] = d.outer_d2pot[ik+1, jk+1]
                    vecb4[k] = d.outer_d2pot[ik, jk+1]

            d.ex[m_p] = -vecpc[m_p] * (
                vecd2[m_p] * (vec1[m_p] - vec2[m_p] -
                                vecz[m_p] * (vecb1[m_p] - vecb2[m_p])) +
                vecb[m_p] * (vec3[m_p] - vec4[m_p] +
                                   vecz[m_p] * (vecb3[m_p] - vecb4[m_p]))
            ) * hx2hy2m1
            d.ey[m_p] = -vecpc[m_p] * (
                vecd1[m_p] * (vec1[m_p] - vec4[m_p] -
                                vecz[m_p] * (vecb1[m_p] - vecb4[m_p])) -
                veca[m_p] * (vec2[m_p] - vec3[m_p] +
                               vecz[m_p] * (vecb2[m_p] - vecb3[m_p]))
            ) * hx2hy2m1

            # Back from the gap frame
            vec_tmp = d.ex[m_p] * cos_dee - d.ey[m_p] * sin_dee
            d.ey[m_p] = d.ex[m_p] * sin_dee + d.ey[m_p] * cos_dee
            d.ex[m_p] = vec_tmp

            d.ez[m_p] = d.z[m_p] * vecpc[m_p] * (
                vecp1[m_p] * vecb1[m_p] + vecp2[m_p] * vecb2[m_p] +
                vecp3[m_p] * vecb3[m_p] + vecp4[m_p] * vecb4[m_p]
            ) * hx2hy2m1

            # Symmetries of the other quadrants
            m1 = m_p & (d.x < 0.0) & (d.y > 0.0)
            d.ex[m1] = d.harmonic_sign * d.ex[m1]
            d.ey[m1] = -d.harmonic_sign * d.ey[m1]

            m2 = m_p & (d.x > 0.0) & (d.y < 0.0)
            d.ey[m2] = -d.ey[m2]

            m3 = m_p & (d.x < 0.0) & (d.y < 0.0)
            d.ex[m3] = d.harmonic_sign * d.ex[m3]
            d.ey[m3] = d.harmonic_sign * d.ey[m3]

            m4 = m_p & (d.x < 0.0)
            d.ez[m4] = -d.harmonic_sign * d.ez[m4]

    add_vertical_sorter_field(d)
    add_extraction_fields(d)


def add_vertical_sorter_field(d):
    """Add the vertical sorter Ez; sorter_voltage_seen is its voltage at particle 0."""
    d.sorter_voltage_seen = 0.0
    if d.sorter_gap == 0.0:
        return
    mask_ver = ((d.r >= d.sorter_rmin) & (d.r <= d.sorter_rmax) &
                (d.theta >= d.sorter_theta1) & (d.theta <= d.sorter_theta2))
    d.ez[mask_ver] = (d.ez[mask_ver] -
        d.sorter_voltage * (1.0 - d.sorter_dv) / d.sorter_gap / 2.0 *
        np.cos(2.0 * d.omega_rf * d.t * (1.0 + d.sorter_nuz / 2.0 / d.harmonic) +
               2.0 * d.sorter_phase - np.pi / 2.0) -
        d.sorter_voltage * (1.0 + d.sorter_dv) / d.sorter_gap / 2.0 *
        np.cos(2.0 * d.omega_rf * d.t * (1.0 - d.sorter_nuz / 2.0 / d.harmonic) +
               2.0 * d.sorter_phase - np.pi / 2.0))
    if (d.r[0] >= d.sorter_rmin and d.r[0] <= d.sorter_rmax and
            d.theta[0] >= d.sorter_theta1 and d.theta[0] <= d.sorter_theta2):
        d.sorter_voltage_seen = -d.ez[0] * d.sorter_gap


def add_extraction_fields(d):
    """Deflector E, dipole B and extraction grid crossings."""
    for i in range(d.n_deflectors):
        check_deflector(d, i)
        if d.n_in_element != 0:
            mask = d.in_element
            vec_p = np.sqrt(d.px[mask]**2 + d.py[mask]**2)
            d.ex[mask] = d.deflector_e[i] * d.py[mask] / vec_p
            d.ey[mask] = -d.deflector_e[i] * d.px[mask] / vec_p

    for i in range(d.n_dipoles):
        check_dipole(d, i)
        if d.n_in_element != 0:
            mask = d.in_element
            vec_dm = np.sqrt((d.x[mask] - d.dipole_x[i])**2 + (d.y[mask] - d.dipole_y[i])**2)
            vec1_dm = vec_dm - d.dipole_r[i]
            d.bx[mask] += d.z[mask] * d.dipole_gradient[i] * (d.x[mask] - d.dipole_x[i]) / vec_dm
            d.by[mask] += d.z[mask] * d.dipole_gradient[i] * (d.y[mask] - d.dipole_y[i]) / vec_dm
            d.bz[mask] += vec1_dm * d.dipole_gradient[i]

    if d.extraction < 0.0 and d.n_grids > 0:
        if d.step >= d.grid_start_step:
            for kg in range(d.n_grids):
                check_grid(d, kg)
                d.grid_reached[d.crossed_grid] = kg + 1  # 1-based in data
                for i in range(d.n_particles):
                    if d.crossed_grid[i]:
                        if f'grid{kg}' in d.outputs:
                            d.outputs[f'grid{kg}'].write(
                                f' {i+1:5d} {d.step:6d}'
                                f' {d.x[i]:15.8e} {d.y[i]:15.8e} {d.z[i]:15.8e}'
                                f' {d.px[i]:15.8e} {d.py[i]:15.8e} {d.pz[i]:15.8e}'
                                f' {d.bz[i]:8.5f}\n')


def compute_bfield(d):
    """Magnetic field only, used while building the initial bunch."""
    hrhtet = d.map_dr * d.map_dtheta
    ivec = np.zeros(d.n_particles, dtype=np.int32)
    jvec = np.zeros(d.n_particles, dtype=np.int32)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    veca = np.zeros(d.n_particles)
    vecb1 = np.zeros(d.n_particles)
    vecb2 = np.zeros(d.n_particles)
    vecb3 = np.zeros(d.n_particles)
    vecb4 = np.zeros(d.n_particles)
    vecb = np.zeros(d.n_particles)
    vecd1 = np.zeros(d.n_particles)
    vecd2 = np.zeros(d.n_particles)
    vecp1 = np.zeros(d.n_particles)
    vecp2 = np.zeros(d.n_particles)
    vecp3 = np.zeros(d.n_particles)
    vecp4 = np.zeros(d.n_particles)
    vecx1 = np.zeros(d.n_particles)
    vecx2 = np.zeros(d.n_particles)
    vecx3 = np.zeros(d.n_particles)
    vecx4 = np.zeros(d.n_particles)
    vecy1 = np.zeros(d.n_particles)
    vecy2 = np.zeros(d.n_particles)
    vecy3 = np.zeros(d.n_particles)
    vecy4 = np.zeros(d.n_particles)
    vecz = np.zeros(d.n_particles)
    d.theta = np.arccos(np.clip(d.x / d.r, -1.0, 1.0))
    d.theta[d.y < 0.0] = 2.0 * np.pi - d.theta[d.y < 0.0]

    mask = d.in_box
    ivec[mask] = ((d.r[mask] - d.map_rmin) / d.map_dr).astype(np.int32)
    jvec[mask] = (d.theta[mask] / d.map_dtheta).astype(np.int32)
    veca[mask] = d.r[mask] - d.map_rmin - ivec[mask] * d.map_dr
    vecb[mask] = d.theta[mask] - jvec[mask] * d.map_dtheta
    vecd1[mask] = veca[mask] - d.map_dr
    vecd2[mask] = vecb[mask] - d.map_dtheta
    vecp1[mask] = vecd1[mask] * vecd2[mask]
    vecp2[mask] = -veca[mask] * vecd2[mask]
    vecp3[mask] = veca[mask] * vecb[mask]
    vecp4[mask] = -vecd1[mask] * vecb[mask]
    vecz[mask] = d.z[mask]**2 / 2.0

    for k in range(d.n_particles):
        if d.in_box[k]:
            ik = ivec[k]
            jk = jvec[k]
            vk = vecz[k]
            vecb1[k] = d.bz_map[ik, jk]
            vecb2[k] = d.bz_map[ik+1, jk]
            vecb3[k] = d.bz_map[ik+1, jk+1]
            vecb4[k] = d.bz_map[ik, jk+1]
            vec1[k] = -vk * d.d2bz_map[ik, jk]
            vec2[k] = -vk * d.d2bz_map[ik+1, jk]
            vec3[k] = -vk * d.d2bz_map[ik+1, jk+1]
            vec4[k] = -vk * d.d2bz_map[ik, jk+1]
            vecx1[k] = d.dbzdx_map[ik, jk]
            vecx2[k] = d.dbzdx_map[ik+1, jk]
            vecx3[k] = d.dbzdx_map[ik+1, jk+1]
            vecx4[k] = d.dbzdx_map[ik, jk+1]
            vecy1[k] = d.dbzdy_map[ik, jk]
            vecy2[k] = d.dbzdy_map[ik+1, jk]
            vecy3[k] = d.dbzdy_map[ik+1, jk+1]
            vecy4[k] = d.dbzdy_map[ik, jk+1]

    vec1[mask] += vecb1[mask]
    vec2[mask] += vecb2[mask]
    vec3[mask] += vecb3[mask]
    vec4[mask] += vecb4[mask]
    d.bz[mask] = (vecp1[mask] * vec1[mask] + vecp2[mask] * vec2[mask] +
                  vecp3[mask] * vec3[mask] + vecp4[mask] * vec4[mask]) / hrhtet
    d.bx[mask] = d.z[mask] * (
        vecp1[mask] * vecx1[mask] + vecp2[mask] * vecx2[mask] +
        vecp3[mask] * vecx3[mask] + vecp4[mask] * vecx4[mask]) / hrhtet
    d.by[mask] = d.z[mask] * (
        vecp1[mask] * vecy1[mask] + vecp2[mask] * vecy2[mask] +
        vecp3[mask] * vecy3[mask] + vecp4[mask] * vecy4[mask]) / hrhtet

    mask_rle = d.r <= d.map_rmin
    d.bx[mask_rle] = 0.0
    d.by[mask_rle] = 0.0
    d.bz[mask_rle] = d.bz_map[0, 0]
