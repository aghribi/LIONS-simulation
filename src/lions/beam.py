"""Initial 6D distribution and phase-space coordinates relative to particle 0."""


import numpy as np

from .apertures import check_box
from .fields import compute_bfield


def cholesky_upper(m, a):
    """Décomposition de Cholesky. Retourne R tel que A = R^T * R."""
    r = np.zeros((m, m), dtype=np.float64)
    r[0, 0] = np.sqrt(a[0, 0])
    r[0, 1:m] = a[0, 1:m] / r[0, 0]
    for i in range(1, m - 1):
        r[i, i] = np.sqrt(a[i, i] - np.sum(r[:i, i]**2))
        vec_im = np.sum(r[:i, i].reshape(-1, 1) * r[:i, i+1:m], axis=0)
        r[i, i+1:m] = (a[i, i+1:m] - vec_im) / r[i, i]
    r[m-1, m-1] = np.sqrt(a[m-1, m-1] - np.sum(r[:m-1, m-1]**2))
    return r


def sample_unit_ball(m, n_part):
    """Génère une distribution aléatoire normalisée dans l'hyper-sphère."""
    np.random.seed(42)
    ps = np.zeros((m, n_part), dtype=np.float64)

    vtab = 10.0 * np.ones((m, n_part - 1), dtype=np.float64)
    while True:
        vec_sum = np.sum(vtab * vtab, axis=0)
        vmax = np.max(vec_sum)
        if vmax <= 1.0:
            break
        vmask = vec_sum > 1.0
        for j in range(m):
            rand_vals = np.random.random(n_part - 1)
            vtab[j, vmask] = 2.0 * rand_vals[vmask] - 1.0

    ps[:, 1:] = vtab
    ps[:, 0] = 0.0

    # Calcul de la matrice de corrélation
    c = np.zeros((m, m), dtype=np.float64)
    for i in range(m):
        for j in range(m):
            c[i, j] = np.sum(ps[i, :] * ps[j, :]) / n_part

    r = cholesky_upper(m, c)
    tr = r.T

    ps[0, :] = ps[0, :] / tr[0, 0]
    for i in range(1, m):
        ps[i, :] = (ps[i, :] - np.sum(
            tr[i, :i].reshape(-1, 1) * ps[:i, :], axis=0
        )) / tr[i, i]

    return ps


def phase_space_coords(d):
    """Calcul des coordonnées dans l'espace de phase."""
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    for k in range(d.n_families):
        nbr = d.family_size[k]
        i1 = d.family_start[k]
        i2 = i1 + nbr
        d.p[i1] = np.sqrt(d.px[i1]**2 + d.py[i1]**2 + d.pz[i1]**2)
        rc = -d.mass_c[k] * d.p[i1] / d.charge[i1] / d.bz[i1]
        d.ps_cartesian[0, :nbr] = d.x[i1:i2] - d.x[i1]
        d.ps_cartesian[4, :nbr] = d.y[i1:i2] - d.y[i1]
        ri1 = np.sqrt(d.x[i1]**2 + d.y[i1]**2)
        tetai1 = np.arccos(d.x[i1] / ri1)
        if d.y[i1] < 0.0:
            tetai1 = 2.0 * np.pi - tetai1
        pri1 = (d.x[i1] * d.px[i1] + d.y[i1] * d.py[i1]) / ri1
        tpri1 = np.arcsin(pri1 / d.p[i1])
        co = np.cos(tetai1 - tpri1)
        si = np.sin(tetai1 - tpri1)
        vec[:nbr] = d.ps_cartesian[0, :nbr] * co + d.ps_cartesian[4, :nbr] * si
        d.ps_cartesian[4, :nbr] = -d.ps_cartesian[0, :nbr] * si + d.ps_cartesian[4, :nbr] * co
        d.ps_cartesian[0, :nbr] = vec[:nbr]
        for ii in range(nbr):
            d.ps[4, ii] = rc * np.arctan(d.ps_cartesian[4, ii] / (rc + d.ps_cartesian[0, ii]))
        vec[:nbr] = np.sin(d.ps[4, :nbr] / 2.0 / rc)**2
        d.ps[0, :nbr] = (d.ps_cartesian[0, :nbr] + 2.0 * rc * vec[:nbr]) / (1.0 - 2.0 * vec[:nbr])
        vec1[:nbr] = np.cos(d.ps[4, :nbr] / rc)
        vec2[:nbr] = np.sin(d.ps[4, :nbr] / rc)
        d.ps_cartesian[1, :nbr] = d.px[i1:i2] - d.px[i1] * vec1[:nbr] + d.py[i1] * vec2[:nbr]
        d.ps_cartesian[5, :nbr] = d.py[i1:i2] - d.px[i1] * vec2[:nbr] - d.py[i1] * vec1[:nbr]
        vec[:nbr] = d.ps_cartesian[1, :nbr] * co + d.ps_cartesian[5, :nbr] * si
        d.ps_cartesian[5, :nbr] = -d.ps_cartesian[1, :nbr] * si + d.ps_cartesian[5, :nbr] * co
        d.ps_cartesian[1, :nbr] = vec[:nbr]
        d.ps[1, :nbr] = d.ps_cartesian[1, :nbr] * vec1[:nbr] + d.ps_cartesian[5, :nbr] * vec2[:nbr]
        d.ps[5, :nbr] = -d.ps_cartesian[1, :nbr] * vec2[:nbr] + d.ps_cartesian[5, :nbr] * vec1[:nbr]
        d.ps[2, :nbr] = d.z[i1:i2] - d.z[i1]
        d.ps[3, :nbr] = d.pz[i1:i2] - d.pz[i1]
        d.ps[5, :nbr] = np.sqrt(
            (d.p[i1] + d.ps[5, :nbr])**2 +
            d.ps[1, :nbr]**2 + d.ps[3, :nbr]**2) - d.p[i1]

        for i in range(6):
            for j in range(6):
                d.sigma[i, j] = np.sum(d.ps[i, :nbr] * d.ps[j, :nbr]) / nbr

        tetai1 = tetai1 * 180.0 / np.pi
        d.dc_transform[:] = 0.0
        d.dc_transform[0, 0] = -1.0
        d.dc_transform[1, 1] = -1.0
        for i in range(2, 6):
            d.dc_transform[i, i] = 1.0
        d.dc_transform[0, 5] = ri1 / d.p[i1]
        d.dc_transform[4, 1] = -d.dc_transform[0, 5]
        d.sigma_dc = d.dc_transform @ d.sigma @ d.dc_transform.T


def project_on_radial_line(d):
    """Calcul des coordonnées sur grille."""
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    vec5 = np.zeros(d.n_particles)
    vec6 = np.zeros(d.n_particles)
    for k in range(d.n_families):
        nbr = d.family_size[k]
        i1 = d.family_start[k]
        i2 = i1 + nbr
        d.p[i1:i2] = np.sqrt(d.px[i1:i2]**2 + d.py[i1:i2]**2 + d.pz[i1:i2]**2)
        vec[:nbr] = -d.mass_c[k] * d.p[i1:i2] / d.charge[i1:i2] / d.bz[i1]
        vec1[i1:i2] = d.x[i1:i2] + d.py[i1:i2] * d.dt / d.boris_e_coef[i1] / d.bz[i1]
        vec2[i1:i2] = d.y[i1:i2] - d.px[i1:i2] * d.dt / d.boris_e_coef[i1] / d.bz[i1]
        ux = (d.x[i1] - vec1[i1]) / vec[0]
        uy = (d.y[i1] - vec2[i1]) / vec[0]
        vec3[i1:i2] = ((d.x[i1] - vec1[i1:i2]) * ux +
                         (d.y[i1] - vec2[i1:i2]) * uy)
        vec4[i1:i2] = ((d.x[i1] - vec1[i1:i2])**2 +
                         (d.y[i1] - vec2[i1:i2])**2 - vec[:nbr]**2)
        d.proj_offset[i1:i2] = -vec3[i1:i2] + np.sqrt(np.abs(vec3[i1:i2]**2 - vec4[i1:i2]))
        d.proj_x[i1:i2] = d.x[i1] + d.proj_offset[i1:i2] * ux
        d.proj_y[i1:i2] = d.y[i1] + d.proj_offset[i1:i2] * uy
        vec3[i1:i2] = np.sqrt((d.proj_x[i1:i2] - vec1[i1:i2])**2 +
                                (d.proj_y[i1:i2] - vec2[i1:i2])**2)
        vec4[i1:i2] = np.sqrt((d.proj_x[i1:i2] - vec1[i1])**2 +
                                (d.proj_y[i1:i2] - vec2[i1])**2)
        vec5[i1:i2] = (
            (d.proj_x[i1:i2] - vec1[i1:i2]) * (d.proj_y[i1:i2] - vec2[i1]) -
            (d.proj_y[i1:i2] - vec2[i1:i2]) * (d.proj_x[i1:i2] - vec1[i1])
        ) / vec3[i1:i2] / vec4[i1:i2]
        d.proj_angle[i1:i2] = np.arcsin(vec5[i1:i2])
        vec4[i1:i2] = np.sqrt((d.x[i1:i2] - vec1[i1:i2])**2 +
                                (d.y[i1:i2] - vec2[i1:i2])**2)
        vec6[i1:i2] = (
            (d.proj_x[i1:i2] - vec1[i1:i2]) * (d.y[i1:i2] - vec2[i1:i2]) -
            (d.proj_y[i1:i2] - vec2[i1:i2]) * (d.x[i1:i2] - vec1[i1:i2])
        ) / vec3[i1:i2] / vec4[i1:i2]
        d.proj_phase[i1:i2] = d.harmonic * np.arcsin(vec6[i1:i2])
    return vec


def projected_phase_space(d):
    """Calcul sur grille pour corrélations."""
    nbr = d.family_size[0]
    vec = project_on_radial_line(d)
    d.ps_projected[0, :nbr] = d.proj_offset[:nbr]
    d.ps_projected[1, :nbr] = d.proj_angle[:nbr] * d.p[0]
    d.ps_projected[2, :nbr] = d.ps[2, :nbr]
    d.ps_projected[3, :nbr] = d.ps[3, :nbr]
    d.ps_projected[4, :nbr] = d.proj_phase[:nbr] * vec[:nbr] / d.harmonic
    d.ps_projected[5, :nbr] = d.ps[5, :nbr]


def distribute_bunch(d, sigma_matrices):
    """Bunch from the sigma matrices (after the field map is read), then optional bunch load."""
    dxtmp = np.zeros((6, d.n_particles))
    rot = np.zeros((6, 6))
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    if d.correlations != 0:
        m = 6
        check_box(d)
        compute_bfield(d)
        for k in range(d.n_families):
            nbr = d.family_size[k]
            if nbr <= 1:
                continue  # une seule particule, pas de distribution a appliquer
            ps = sample_unit_ball(m, nbr)
            d.ps[:, :nbr] = ps
            rcor = cholesky_upper(m, sigma_matrices[k])
            trcor = rcor.T
            dxtmp[:6, :nbr] = 0.0
            for i in range(6):
                for j in range(6):
                    dxtmp[i, :nbr] += trcor[i, j] * d.ps[j, :nbr]
            d.ps[:6, :nbr] = dxtmp[:6, :nbr]

            i1 = d.family_start[k]
            i2 = i1 + d.family_size[k]
            rc = -d.mass_c[k] * d.p[i1] / d.charge[i1] / d.bz[i1]

            if d.r0[k] > d.central_rmax:
                print('prise en compte de sqrt(cos(phi))...')
                vec4[:nbr] = np.sqrt(np.cos(d.harmonic * d.ps[4, :nbr] / rc))
                d.ps[0, :nbr] = d.ps[0, :nbr] * vec4[:nbr] - rc * (1.0 - vec4[:nbr])
            else:
                print('pas de sqrt(cos(phi))...')

            vec1[:nbr] = np.cos(d.ps[4, :nbr] / rc)
            vec2[:nbr] = np.sin(d.ps[4, :nbr] / rc)
            vec3[:nbr] = np.sin(d.ps[4, :nbr] / rc / 2.0)**2
            d.ps_cartesian[0, :nbr] = d.ps[0, :nbr] - 2.0 * (rc + d.ps[0, :nbr]) * vec3[:nbr]
            d.ps[5, :nbr] = np.sqrt(
                (d.p[i1] + d.ps[5, :nbr])**2 -
                d.ps[1, :nbr]**2 - d.ps[3, :nbr]**2) - d.p[i1]
            d.ps_cartesian[4, :nbr] = (rc + d.ps[0, :nbr]) * vec2[:nbr]
            d.ps_cartesian[1, :nbr] = d.ps[1, :nbr] * vec1[:nbr] - d.ps[5, :nbr] * vec2[:nbr]
            d.ps_cartesian[5, :nbr] = d.ps[1, :nbr] * vec2[:nbr] + d.ps[5, :nbr] * vec1[:nbr]
            d.ps_cartesian[2:4, :nbr] = d.ps[2:4, :nbr]

            co = np.cos(d.theta0[k] - d.pr0[k])
            si = np.sin(d.theta0[k] - d.pr0[k])
            rot[:] = 0.0
            rot[0, 0] = co
            rot[0, 4] = -si
            rot[1, 1] = co
            rot[1, 5] = -si
            rot[2, 2] = 1.0
            rot[3, 3] = 1.0
            rot[4, 0] = si
            rot[4, 4] = co
            rot[5, 1] = si
            rot[5, 5] = co
            dxtmp[:6, :nbr] = 0.0
            for i in range(6):
                for j in range(6):
                    dxtmp[i, :nbr] += rot[i, j] * d.ps_cartesian[j, :nbr]
            d.ps_cartesian[:6, :nbr] = dxtmp[:6, :nbr]

            d.x[i1:i2] += d.ps_cartesian[0, :nbr]
            d.y[i1:i2] += d.ps_cartesian[4, :nbr]
            vec[:nbr] = d.px[i1:i2] * vec1[:nbr] - d.py[i1:i2] * vec2[:nbr]
            d.py[i1:i2] = d.px[i1:i2] * vec2[:nbr] + d.py[i1:i2] * vec1[:nbr] + d.ps_cartesian[5, :nbr]
            d.px[i1:i2] = vec[:nbr] + d.ps_cartesian[1, :nbr]
            if d.phase_space[:2] in ('3d', '6d'):
                d.z[i1:i2] += d.ps_cartesian[2, :nbr]
                d.pz[i1:i2] += d.ps_cartesian[3, :nbr]
            if d.r0[k] > d.central_rmax:
                d.px[i1:i2] *= vec4[:nbr]
                d.py[i1:i2] *= vec4[:nbr]
                d.pz[i1:i2] *= vec4[:nbr]

        if d.bunch_load_file != 'no':
            print('attention : chargement paquet d ions...')
            with open(d.bunch_load_file, 'rb') as f26:
                for i in range(d.n_particles):
                    data = np.fromfile(f26, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.x[i], d.y[i], d.z[i] = data
                for i in range(d.n_particles):
                    data = np.fromfile(f26, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.px[i], d.py[i], d.pz[i] = data
