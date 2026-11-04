"""Initial 6D distribution and phase-space coordinates relative to
particle 0.
"""

import numpy as np

from .apertures import check_box
from .fields import compute_bfield


def cholesky_upper(m, a):
    """Cholesky factor: upper-triangular R with A = R^T R (m x m)."""
    r = np.zeros((m, m), dtype=np.float64)
    r[0, 0] = np.sqrt(a[0, 0])
    r[0, 1:m] = a[0, 1:m] / r[0, 0]
    for i in range(1, m - 1):
        r[i, i] = np.sqrt(a[i, i] - np.sum(r[:i, i] ** 2))
        vec_im = np.sum(r[:i, i].reshape(-1, 1) * r[:i, i + 1 : m], axis=0)
        r[i, i + 1 : m] = (a[i, i + 1 : m] - vec_im) / r[i, i]
    r[m - 1, m - 1] = np.sqrt(a[m - 1, m - 1] - np.sum(r[: m - 1, m - 1] ** 2))
    return r


def sample_unit_ball(m, n_part):
    """n_part points uniform in the unit m-ball, point 0 at the origin,
    whitened to unit second moments. Fixed seed 42.
    """
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

    c = np.zeros((m, m), dtype=np.float64)
    for i in range(m):
        for j in range(m):
            c[i, j] = np.sum(ps[i, :] * ps[j, :]) / n_part

    r = cholesky_upper(m, c)
    tr = r.T

    ps[0, :] = ps[0, :] / tr[0, 0]
    for i in range(1, m):
        ps[i, :] = (
            ps[i, :] - np.sum(tr[i, :i].reshape(-1, 1) * ps[:i, :], axis=0)
        ) / tr[i, i]

    return ps


def phase_space_coords(bunch):
    """Phase-space coordinates of each family relative to its first
    particle, and sigma matrices.
    """
    vec = np.zeros(bunch.n_particles)
    vec1 = np.zeros(bunch.n_particles)
    vec2 = np.zeros(bunch.n_particles)
    for k in range(len(bunch.families)):
        nbr = bunch.families[k].size
        i1 = bunch.families[k].start
        i2 = i1 + nbr
        bunch.p[i1] = np.sqrt(
            bunch.px[i1] ** 2 + bunch.py[i1] ** 2 + bunch.pz[i1] ** 2
        )
        rc = -bunch.mass_c[k] * bunch.p[i1] / bunch.charge[i1] / bunch.bz[i1]
        bunch.ps_cartesian[0, :nbr] = bunch.x[i1:i2] - bunch.x[i1]
        bunch.ps_cartesian[4, :nbr] = bunch.y[i1:i2] - bunch.y[i1]
        r_ref = np.sqrt(bunch.x[i1] ** 2 + bunch.y[i1] ** 2)
        theta_ref = np.arccos(bunch.x[i1] / r_ref)
        if bunch.y[i1] < 0.0:
            theta_ref = 2.0 * np.pi - theta_ref
        pr_ref = (
            bunch.x[i1] * bunch.px[i1] + bunch.y[i1] * bunch.py[i1]
        ) / r_ref
        pr_angle_ref = np.arcsin(pr_ref / bunch.p[i1])
        co = np.cos(theta_ref - pr_angle_ref)
        si = np.sin(theta_ref - pr_angle_ref)
        vec[:nbr] = (
            bunch.ps_cartesian[0, :nbr] * co + bunch.ps_cartesian[4, :nbr] * si
        )
        bunch.ps_cartesian[4, :nbr] = (
            -bunch.ps_cartesian[0, :nbr] * si
            + bunch.ps_cartesian[4, :nbr] * co
        )
        bunch.ps_cartesian[0, :nbr] = vec[:nbr]
        for ii in range(nbr):
            bunch.ps[4, ii] = rc * np.arctan(
                bunch.ps_cartesian[4, ii] / (rc + bunch.ps_cartesian[0, ii])
            )
        vec[:nbr] = np.sin(bunch.ps[4, :nbr] / 2.0 / rc) ** 2
        bunch.ps[0, :nbr] = (
            bunch.ps_cartesian[0, :nbr] + 2.0 * rc * vec[:nbr]
        ) / (1.0 - 2.0 * vec[:nbr])
        vec1[:nbr] = np.cos(bunch.ps[4, :nbr] / rc)
        vec2[:nbr] = np.sin(bunch.ps[4, :nbr] / rc)
        bunch.ps_cartesian[1, :nbr] = (
            bunch.px[i1:i2]
            - bunch.px[i1] * vec1[:nbr]
            + bunch.py[i1] * vec2[:nbr]
        )
        bunch.ps_cartesian[5, :nbr] = (
            bunch.py[i1:i2]
            - bunch.px[i1] * vec2[:nbr]
            - bunch.py[i1] * vec1[:nbr]
        )
        vec[:nbr] = (
            bunch.ps_cartesian[1, :nbr] * co + bunch.ps_cartesian[5, :nbr] * si
        )
        bunch.ps_cartesian[5, :nbr] = (
            -bunch.ps_cartesian[1, :nbr] * si
            + bunch.ps_cartesian[5, :nbr] * co
        )
        bunch.ps_cartesian[1, :nbr] = vec[:nbr]
        bunch.ps[1, :nbr] = (
            bunch.ps_cartesian[1, :nbr] * vec1[:nbr]
            + bunch.ps_cartesian[5, :nbr] * vec2[:nbr]
        )
        bunch.ps[5, :nbr] = (
            -bunch.ps_cartesian[1, :nbr] * vec2[:nbr]
            + bunch.ps_cartesian[5, :nbr] * vec1[:nbr]
        )
        bunch.ps[2, :nbr] = bunch.z[i1:i2] - bunch.z[i1]
        bunch.ps[3, :nbr] = bunch.pz[i1:i2] - bunch.pz[i1]
        bunch.ps[5, :nbr] = (
            np.sqrt(
                (bunch.p[i1] + bunch.ps[5, :nbr]) ** 2
                + bunch.ps[1, :nbr] ** 2
                + bunch.ps[3, :nbr] ** 2
            )
            - bunch.p[i1]
        )

        for i in range(6):
            for j in range(6):
                bunch.sigma[i, j] = (
                    np.sum(bunch.ps[i, :nbr] * bunch.ps[j, :nbr]) / nbr
                )

        bunch.dc_transform[:] = 0.0
        bunch.dc_transform[0, 0] = -1.0
        bunch.dc_transform[1, 1] = -1.0
        for i in range(2, 6):
            bunch.dc_transform[i, i] = 1.0
        bunch.dc_transform[0, 5] = r_ref / bunch.p[i1]
        bunch.dc_transform[4, 1] = -bunch.dc_transform[0, 5]
        bunch.sigma_dc = (
            bunch.dc_transform @ bunch.sigma @ bunch.dc_transform.T
        )


def project_on_radial_line(bunch, rf, run):
    """Project each particle on the radial line through particle 0 of
    its family. Returns the curvature radii of the last family.
    """
    vec = np.zeros(bunch.n_particles)
    vec1 = np.zeros(bunch.n_particles)
    vec2 = np.zeros(bunch.n_particles)
    vec3 = np.zeros(bunch.n_particles)
    vec4 = np.zeros(bunch.n_particles)
    vec5 = np.zeros(bunch.n_particles)
    vec6 = np.zeros(bunch.n_particles)
    for k in range(len(bunch.families)):
        nbr = bunch.families[k].size
        i1 = bunch.families[k].start
        i2 = i1 + nbr
        bunch.p[i1:i2] = np.sqrt(
            bunch.px[i1:i2] ** 2 + bunch.py[i1:i2] ** 2 + bunch.pz[i1:i2] ** 2
        )
        vec[:nbr] = (
            -bunch.mass_c[k]
            * bunch.p[i1:i2]
            / bunch.charge[i1:i2]
            / bunch.bz[i1]
        )
        vec1[i1:i2] = (
            bunch.x[i1:i2]
            + bunch.py[i1:i2] * run.dt / bunch.boris_e_coef[i1] / bunch.bz[i1]
        )
        vec2[i1:i2] = (
            bunch.y[i1:i2]
            - bunch.px[i1:i2] * run.dt / bunch.boris_e_coef[i1] / bunch.bz[i1]
        )
        ux = (bunch.x[i1] - vec1[i1]) / vec[0]
        uy = (bunch.y[i1] - vec2[i1]) / vec[0]
        vec3[i1:i2] = (bunch.x[i1] - vec1[i1:i2]) * ux + (
            bunch.y[i1] - vec2[i1:i2]
        ) * uy
        vec4[i1:i2] = (
            (bunch.x[i1] - vec1[i1:i2]) ** 2
            + (bunch.y[i1] - vec2[i1:i2]) ** 2
            - vec[:nbr] ** 2
        )
        bunch.proj_offset[i1:i2] = -vec3[i1:i2] + np.sqrt(
            np.abs(vec3[i1:i2] ** 2 - vec4[i1:i2])
        )
        bunch.proj_x[i1:i2] = bunch.x[i1] + bunch.proj_offset[i1:i2] * ux
        bunch.proj_y[i1:i2] = bunch.y[i1] + bunch.proj_offset[i1:i2] * uy
        vec3[i1:i2] = np.sqrt(
            (bunch.proj_x[i1:i2] - vec1[i1:i2]) ** 2
            + (bunch.proj_y[i1:i2] - vec2[i1:i2]) ** 2
        )
        vec4[i1:i2] = np.sqrt(
            (bunch.proj_x[i1:i2] - vec1[i1]) ** 2
            + (bunch.proj_y[i1:i2] - vec2[i1]) ** 2
        )
        vec5[i1:i2] = (
            (
                (bunch.proj_x[i1:i2] - vec1[i1:i2])
                * (bunch.proj_y[i1:i2] - vec2[i1])
                - (bunch.proj_y[i1:i2] - vec2[i1:i2])
                * (bunch.proj_x[i1:i2] - vec1[i1])
            )
            / vec3[i1:i2]
            / vec4[i1:i2]
        )
        bunch.proj_angle[i1:i2] = np.arcsin(vec5[i1:i2])
        vec4[i1:i2] = np.sqrt(
            (bunch.x[i1:i2] - vec1[i1:i2]) ** 2
            + (bunch.y[i1:i2] - vec2[i1:i2]) ** 2
        )
        vec6[i1:i2] = (
            (
                (bunch.proj_x[i1:i2] - vec1[i1:i2])
                * (bunch.y[i1:i2] - vec2[i1:i2])
                - (bunch.proj_y[i1:i2] - vec2[i1:i2])
                * (bunch.x[i1:i2] - vec1[i1:i2])
            )
            / vec3[i1:i2]
            / vec4[i1:i2]
        )
        bunch.proj_phase[i1:i2] = rf.harmonic * np.arcsin(vec6[i1:i2])
    return vec


def projected_phase_space(bunch, rf, run):
    """Phase-space coordinates on the radial line through particle 0,
    for the correlations.
    """
    nbr = bunch.families[0].size
    vec = project_on_radial_line(bunch, rf, run)
    bunch.ps_projected[0, :nbr] = bunch.proj_offset[:nbr]
    bunch.ps_projected[1, :nbr] = bunch.proj_angle[:nbr] * bunch.p[0]
    bunch.ps_projected[2, :nbr] = bunch.ps[2, :nbr]
    bunch.ps_projected[3, :nbr] = bunch.ps[3, :nbr]
    bunch.ps_projected[4, :nbr] = (
        bunch.proj_phase[:nbr] * vec[:nbr] / rf.harmonic
    )
    bunch.ps_projected[5, :nbr] = bunch.ps[5, :nbr]


def distribute_bunch(d, sigma_matrices):
    """Bunch from the sigma matrices (after the field map is read), then
    optional bunch load.
    """
    dxtmp = np.zeros((6, d.bunch.n_particles))
    rot = np.zeros((6, 6))
    vec = np.zeros(d.bunch.n_particles)
    vec1 = np.zeros(d.bunch.n_particles)
    vec2 = np.zeros(d.bunch.n_particles)
    vec3 = np.zeros(d.bunch.n_particles)
    vec4 = np.zeros(d.bunch.n_particles)
    if d.bunch.correlations != 0:
        m = 6
        check_box(d.bunch, d.field, d.box)
        compute_bfield(d.bunch, d.field)
        for k in range(len(d.bunch.families)):
            nbr = d.bunch.families[k].size
            if nbr <= 1:
                continue  # single particle, nothing to distribute
            ps = sample_unit_ball(m, nbr)
            d.bunch.ps[:, :nbr] = ps
            rcor = cholesky_upper(m, sigma_matrices[k])
            trcor = rcor.T
            dxtmp[:6, :nbr] = 0.0
            for i in range(6):
                for j in range(6):
                    dxtmp[i, :nbr] += trcor[i, j] * d.bunch.ps[j, :nbr]
            d.bunch.ps[:6, :nbr] = dxtmp[:6, :nbr]

            i1 = d.bunch.families[k].start
            i2 = i1 + d.bunch.families[k].size
            rc = (
                -d.bunch.mass_c[k]
                * d.bunch.p[i1]
                / d.bunch.charge[i1]
                / d.bunch.bz[i1]
            )

            if d.bunch.families[k].r0 > d.central.rmax:
                print('sqrt(cos(phi)) correction applied')
                vec4[:nbr] = np.sqrt(
                    np.cos(d.rf.harmonic * d.bunch.ps[4, :nbr] / rc)
                )
                d.bunch.ps[0, :nbr] = d.bunch.ps[0, :nbr] * vec4[:nbr] - rc * (
                    1.0 - vec4[:nbr]
                )
            else:
                print('no sqrt(cos(phi)) correction')

            vec1[:nbr] = np.cos(d.bunch.ps[4, :nbr] / rc)
            vec2[:nbr] = np.sin(d.bunch.ps[4, :nbr] / rc)
            vec3[:nbr] = np.sin(d.bunch.ps[4, :nbr] / rc / 2.0) ** 2
            d.bunch.ps_cartesian[0, :nbr] = (
                d.bunch.ps[0, :nbr]
                - 2.0 * (rc + d.bunch.ps[0, :nbr]) * vec3[:nbr]
            )
            d.bunch.ps[5, :nbr] = (
                np.sqrt(
                    (d.bunch.p[i1] + d.bunch.ps[5, :nbr]) ** 2
                    - d.bunch.ps[1, :nbr] ** 2
                    - d.bunch.ps[3, :nbr] ** 2
                )
                - d.bunch.p[i1]
            )
            d.bunch.ps_cartesian[4, :nbr] = (rc + d.bunch.ps[0, :nbr]) * vec2[
                :nbr
            ]
            d.bunch.ps_cartesian[1, :nbr] = (
                d.bunch.ps[1, :nbr] * vec1[:nbr]
                - d.bunch.ps[5, :nbr] * vec2[:nbr]
            )
            d.bunch.ps_cartesian[5, :nbr] = (
                d.bunch.ps[1, :nbr] * vec2[:nbr]
                + d.bunch.ps[5, :nbr] * vec1[:nbr]
            )
            d.bunch.ps_cartesian[2:4, :nbr] = d.bunch.ps[2:4, :nbr]

            co = np.cos(d.bunch.families[k].theta0 - d.bunch.families[k].pr0)
            si = np.sin(d.bunch.families[k].theta0 - d.bunch.families[k].pr0)
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
                    dxtmp[i, :nbr] += rot[i, j] * d.bunch.ps_cartesian[j, :nbr]
            d.bunch.ps_cartesian[:6, :nbr] = dxtmp[:6, :nbr]

            d.bunch.x[i1:i2] += d.bunch.ps_cartesian[0, :nbr]
            d.bunch.y[i1:i2] += d.bunch.ps_cartesian[4, :nbr]
            vec[:nbr] = (
                d.bunch.px[i1:i2] * vec1[:nbr] - d.bunch.py[i1:i2] * vec2[:nbr]
            )
            d.bunch.py[i1:i2] = (
                d.bunch.px[i1:i2] * vec2[:nbr]
                + d.bunch.py[i1:i2] * vec1[:nbr]
                + d.bunch.ps_cartesian[5, :nbr]
            )
            d.bunch.px[i1:i2] = vec[:nbr] + d.bunch.ps_cartesian[1, :nbr]
            if d.bunch.phase_space[:2] in ('3d', '6d'):
                d.bunch.z[i1:i2] += d.bunch.ps_cartesian[2, :nbr]
                d.bunch.pz[i1:i2] += d.bunch.ps_cartesian[3, :nbr]
            if d.bunch.families[k].r0 > d.central.rmax:
                d.bunch.px[i1:i2] *= vec4[:nbr]
                d.bunch.py[i1:i2] *= vec4[:nbr]
                d.bunch.pz[i1:i2] *= vec4[:nbr]

        if d.files.bunch_load != 'no':
            print('warning: loading the ion bunch')
            with open(d.files.bunch_load, 'rb') as fbunch:
                for i in range(d.bunch.n_particles):
                    data = np.fromfile(fbunch, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.bunch.x[i], d.bunch.y[i], d.bunch.z[i] = data
                for i in range(d.bunch.n_particles):
                    data = np.fromfile(fbunch, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.bunch.px[i], d.bunch.py[i], d.bunch.pz[i] = data
