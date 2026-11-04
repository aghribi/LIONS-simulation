"""Text outputs: z and xy tracks of particle 0, correlations, end-of-run
summary.
"""

import struct

import numpy as np

from .constants import C_LIGHT, MEV, SQRT8


def write_outputs(d):
    """Write the z, xy and correlation outputs for the current step, and
    the end-of-run summary.
    """
    vec = np.zeros(d.bunch.n_particles)
    vec1 = np.zeros(d.bunch.n_particles)
    vec2 = np.zeros(d.bunch.n_particles)

    if not d.outputs.headers_written:
        d.outputs.headers_written = True

        if d.files.z_output != 'no' and 'z' in d.outputs.open:
            fz = d.outputs.open['z']
            fz.write('LIONS program\n')
            fz.write('vertical (z) motion\n')

        if d.files.xy_output != 'no' and 'xy' in d.outputs.open:
            fxy = d.outputs.open['xy']
            fxy.write('LIONS program\n')
            fxy.write(f'{d.run.title}\n')
            fxy.write(
                f' atomic mass           :    {d.bunch.families[0].atomic_mass:11.7f}\n'
            )
            fxy.write(
                f' A                     :          {d.bunch.families[0].mass_number:5.1f}\n'
            )
            fxy.write(
                f' q                     :          {d.bunch.families[0].charge_state:5.1f}\n'
            )
            r0mm = d.bunch.families[0].r0 * 1000.0
            theta0_deg = d.bunch.families[0].theta0 * 180.0 / np.pi
            z0mm = d.bunch.families[0].z0 * 1000.0
            pr0m = d.bunch.families[0].pr0 * 1000.0
            pz0m = d.bunch.families[0].pz0 * 1000.0
            fxy.write(f' r0     : {r0mm:6.1f}\n')
            fxy.write(f' theta0 : {theta0_deg:6.1f}\n')
            fxy.write(f' z0     : {z0mm:6.1f}\n')
            fxy.write(f' pr0    : {pr0m:8.2f}\n')
            fxy.write(f' pz0    : {pz0m:8.2f}\n')
            fxy.write(f' w0     : {d.bunch.w0_mev_per_u:14.10f}\n')
            fxy.write(f' phi0   : {d.bunch.families[0].phi0_deg:11.2f}\n')
            fxy.write(f' fhf    : {d.rf.frequency:9.0f}\n')
            fxy.write(f' harm   : {d.rf.harmonic:9.0f}\n')
            fxy.write(f' v_right: {-d.rf.v_right:9.1f}\n')
            fxy.write(f' v_left : {-d.rf.v_left:9.1f}\n')
            fxy.write(f' db/b   : {d.field.db_over_b:9.6f}\n')
            if d.field.b_factor is not None:
                fxy.write(f' B fact : {d.field.b_factor:9.6f} T\n')
            dtnano = d.run.dt * 1e9
            fxy.write(f' n part : {d.bunch.n_particles:14d}\n')
            fxy.write(f' n steps: {d.run.n_steps:14d}\n')
            fxy.write(f' n/turn : {d.run.steps_per_turn:14d}\n')
            fxy.write(f' dt(ns) : {dtnano:9.6f}\n')
            fxy.write(
                '   step  turn  theta     x       y      xc'
                '      yc        bz '
                '    es      er           r      tpr'
                '       w           dphi\n'
            )

    d.run.theta_prev_deg = d.run.theta_deg
    r = np.sqrt(d.bunch.x[0] ** 2 + d.bunch.y[0] ** 2)
    d.run.theta_deg = np.arccos(np.clip(d.bunch.x[0] / r, -1.0, 1.0))
    if d.bunch.y[0] < 0.0:
        d.run.theta_deg = 2.0 * np.pi - d.run.theta_deg
    d.run.theta_deg = 180.0 / np.pi * d.run.theta_deg

    if d.run.steps_per_turn > 0:
        if d.run.theta_deg < d.run.theta_prev_deg:
            d.run.turn += 1
    else:
        if d.run.theta_deg > d.run.theta_prev_deg:
            d.run.turn -= 1

    d.run.theta_total_deg = d.run.theta_deg + 360.0 * (d.run.turn - 1)

    if (
        d.run.step % d.run.output_every != 0
        and d.run.step != 1
        and d.run.step != d.run.n_steps
        and d.run.step != (d.run.n_steps - 1)
    ):
        return

    rf_phase_deg = 180.0 / np.pi * d.rf.omega * d.run.t
    phase_lag_deg = rf_phase_deg - d.rf.harmonic * d.run.theta_total_deg - 90.0
    w_mev_per_u = (
        d.bunch.rest_energy[0]
        * (d.bunch.gamma[0] - 1.0)
        / MEV
        / d.bunch.families[0].mass_number
    )
    xmm = d.bunch.x[0] * 1000.0
    ymm = d.bunch.y[0] * 1000.0
    zmm = d.bunch.z[0] * 1000.0
    rmm = r * 1000.0
    pp = np.sqrt(d.bunch.px[0] ** 2 + d.bunch.py[0] ** 2)
    ppp = np.sqrt(d.bunch.px[0] ** 2 + d.bunch.py[0] ** 2 + d.bunch.pz[0] ** 2)
    pr_val = (xmm * d.bunch.px[0] + ymm * d.bunch.py[0]) / rmm
    tpr = np.arcsin(pr_val / ppp) * 1000.0
    xc1mm = (
        d.bunch.x[0]
        + d.bunch.py[0] * d.run.dt / d.bunch.boris_e_coef[0] / d.bunch.bz[0]
    ) * 1000.0
    yc1mm = (
        d.bunch.y[0]
        - d.bunch.px[0] * d.run.dt / d.bunch.boris_e_coef[0] / d.bunch.bz[0]
    ) * 1000.0
    eskvcm = (
        (d.bunch.ex[0] * d.bunch.px[0] + d.bunch.ey[0] * d.bunch.py[0])
        / pp
        / 100000.0
    )
    erkvcm = (
        (d.bunch.ex[0] * d.bunch.py[0] - d.bunch.ey[0] * d.bunch.px[0])
        / pp
        / 100000.0
    )
    ezkvcm = d.bunch.ez[0] / 100000.0
    tpz = np.arcsin(d.bunch.pz[0] / ppp) * 1000.0

    if d.run.step % d.run.output_every != 0 and d.run.step != 1:
        pass  # skip normal output, go to end-of-run output
    else:
        if d.files.z_output != 'no' and 'z' in d.outputs.open:
            if abs(d.bunch.z[0]) > d.sorter.nuz_z_cut:
                nuz2 = (
                    -1.0
                    / d.bunch.mass[0]
                    / d.bunch.gamma[0]
                    * (d.rf.harmonic / d.rf.omega) ** 2
                    / d.bunch.z[0]
                    * d.bunch.charge[0]
                    * (
                        d.bunch.ez[0]
                        + C_LIGHT
                        / d.bunch.gamma[0]
                        * (
                            d.bunch.px[0] * d.bunch.by[0]
                            - d.bunch.py[0] * d.bunch.bx[0]
                        )
                    )
                )
                d.outputs.nuz = np.sqrt(abs(nuz2))
            else:
                d.outputs.nuz = 0.0
            d.outputs.open['z'].write(
                f' {d.run.step:6d}  {d.run.turn:3d}  {d.run.theta_deg:5.1f}'
                f' {xmm:8.2f} {ymm:8.2f} {zmm:8.2f} {rmm:8.2f}'
                f'   {d.bunch.bx[0]:7.4f} {d.bunch.by[0]:7.4f} {d.bunch.bz[0]:7.4f}'
                f' {ezkvcm:7.3f}     {tpz:6.1f} {d.outputs.nuz:12.5f}'
                f' {d.sorter.voltage_seen:12.5f}\n'
            )

        if d.files.xy_output != 'no' and 'xy' in d.outputs.open:
            d.outputs.open['xy'].write(
                f' {d.run.step:7d}  {d.run.turn:3d}  {d.run.theta_deg:5.1f}'
                f' {xmm:7.1f} {ymm:7.1f} {xc1mm:7.1f} {yc1mm:7.1f}'
                f'   {d.bunch.bz[0]:7.4f} {eskvcm:7.3f} {erkvcm:7.3f}'
                f'     {rmm:7.1f}  {tpr:6.1f} {w_mev_per_u:9.5f} {phase_lag_deg:15.2f}\n'
            )

    # Correlations, skipped together with the z/xy rows (Fortran goto
    # 800)
    if (
        (d.run.step % d.run.output_every == 0 or d.run.step == 1)
        and d.bunch.correlations != 0
        and not (
            d.files.correlation == 'no'
            and d.files.matrix == 'no'
            and d.files.dc_matrix == 'no'
        )
    ):
        for k in range(len(d.bunch.families)):
            nbr = d.bunch.families[k].size
            i1 = d.bunch.families[k].start
            i2 = i1 + nbr
            d.bunch.p[i1] = np.sqrt(
                d.bunch.px[i1] ** 2 + d.bunch.py[i1] ** 2 + d.bunch.pz[i1] ** 2
            )
            rc = (
                -d.bunch.mass_c[k]
                * d.bunch.p[i1]
                / d.bunch.charge[i1]
                / d.bunch.bz[i1]
            )
            # Same phase-space transform as beam.phase_space_coords
            # (Fortran caldxs)
            d.bunch.ps_cartesian[0, :nbr] = d.bunch.x[i1:i2] - d.bunch.x[i1]
            d.bunch.ps_cartesian[4, :nbr] = d.bunch.y[i1:i2] - d.bunch.y[i1]
            r_ref = np.sqrt(d.bunch.x[i1] ** 2 + d.bunch.y[i1] ** 2)
            theta_ref = np.arccos(np.clip(d.bunch.x[i1] / r_ref, -1.0, 1.0))
            if d.bunch.y[i1] < 0.0:
                theta_ref = 2.0 * np.pi - theta_ref
            pr_ref = (
                d.bunch.x[i1] * d.bunch.px[i1] + d.bunch.y[i1] * d.bunch.py[i1]
            ) / r_ref
            pr_angle_ref = np.arcsin(pr_ref / d.bunch.p[i1])
            co = np.cos(theta_ref - pr_angle_ref)
            si = np.sin(theta_ref - pr_angle_ref)
            vec[:nbr] = (
                d.bunch.ps_cartesian[0, :nbr] * co
                + d.bunch.ps_cartesian[4, :nbr] * si
            )
            d.bunch.ps_cartesian[4, :nbr] = (
                -d.bunch.ps_cartesian[0, :nbr] * si
                + d.bunch.ps_cartesian[4, :nbr] * co
            )
            d.bunch.ps_cartesian[0, :nbr] = vec[:nbr]
            for ii in range(nbr):
                d.bunch.ps[4, ii] = rc * np.arctan(
                    d.bunch.ps_cartesian[4, ii]
                    / (rc + d.bunch.ps_cartesian[0, ii])
                )
            vec[:nbr] = np.sin(d.bunch.ps[4, :nbr] / 2.0 / rc) ** 2
            d.bunch.ps[0, :nbr] = (
                d.bunch.ps_cartesian[0, :nbr] + 2.0 * rc * vec[:nbr]
            ) / (1.0 - 2.0 * vec[:nbr])
            vec1[:nbr] = np.cos(d.bunch.ps[4, :nbr] / rc)
            vec2[:nbr] = np.sin(d.bunch.ps[4, :nbr] / rc)
            d.bunch.ps_cartesian[1, :nbr] = (
                d.bunch.px[i1:i2]
                - d.bunch.px[i1] * vec1[:nbr]
                + d.bunch.py[i1] * vec2[:nbr]
            )
            d.bunch.ps_cartesian[5, :nbr] = (
                d.bunch.py[i1:i2]
                - d.bunch.px[i1] * vec2[:nbr]
                - d.bunch.py[i1] * vec1[:nbr]
            )
            vec[:nbr] = (
                d.bunch.ps_cartesian[1, :nbr] * co
                + d.bunch.ps_cartesian[5, :nbr] * si
            )
            d.bunch.ps_cartesian[5, :nbr] = (
                -d.bunch.ps_cartesian[1, :nbr] * si
                + d.bunch.ps_cartesian[5, :nbr] * co
            )
            d.bunch.ps_cartesian[1, :nbr] = vec[:nbr]
            d.bunch.ps[1, :nbr] = (
                d.bunch.ps_cartesian[1, :nbr] * vec1[:nbr]
                + d.bunch.ps_cartesian[5, :nbr] * vec2[:nbr]
            )
            d.bunch.ps[5, :nbr] = (
                -d.bunch.ps_cartesian[1, :nbr] * vec2[:nbr]
                + d.bunch.ps_cartesian[5, :nbr] * vec1[:nbr]
            )
            d.bunch.ps[2, :nbr] = d.bunch.z[i1:i2] - d.bunch.z[i1]
            d.bunch.ps[3, :nbr] = d.bunch.pz[i1:i2] - d.bunch.pz[i1]
            d.bunch.ps[5, :nbr] = (
                np.sqrt(
                    (d.bunch.p[i1] + d.bunch.ps[5, :nbr]) ** 2
                    + d.bunch.ps[1, :nbr] ** 2
                    + d.bunch.ps[3, :nbr] ** 2
                )
                - d.bunch.p[i1]
            )

            d.bunch.sigma[:] = 0.0
            for i in range(6):
                for j in range(6):
                    d.bunch.sigma[i, j] = (
                        np.sum(
                            d.bunch.ps[i, :nbr]
                            * d.bunch.ps[j, :nbr]
                            * d.bunch.alive[:nbr].astype(float)
                        )
                        / d.bunch.n_alive
                    )

            theta_ref_deg = theta_ref * 180.0 / np.pi
            r_ref_mm = r_ref * 1000.0
            d.bunch.dc_transform[:] = 0.0
            d.bunch.dc_transform[0, 0] = -1.0
            d.bunch.dc_transform[1, 1] = -1.0
            for i in range(2, 6):
                d.bunch.dc_transform[i, i] = 1.0
            d.bunch.dc_transform[0, 5] = r_ref / d.bunch.p[i1]
            d.bunch.dc_transform[4, 1] = -d.bunch.dc_transform[0, 5]
            d.bunch.sigma_dc = (
                d.bunch.dc_transform @ d.bunch.sigma @ d.bunch.dc_transform.T
            )

            if d.files.correlation != 'no' and 'correlation' in d.outputs.open:
                gamma_ref = np.sqrt(1.0 + d.bunch.p[i1] ** 2)
                dr = np.sqrt(d.bunch.sigma[0, 0]) * SQRT8
                drmm = dr * 1000.0
                ds = np.sqrt(d.bunch.sigma[4, 4]) * SQRT8
                dsmm = ds * 1000.0
                dz_val = np.sqrt(d.bunch.sigma[2, 2]) * SQRT8
                dzmm = dz_val * 1000.0
                dphi = (
                    ds
                    * d.rf.omega
                    / C_LIGHT
                    / d.bunch.p[i1]
                    * 180.0
                    / np.pi
                    * gamma_ref
                )
                drdpr = 8.0 * np.sqrt(
                    abs(
                        d.bunch.sigma[0, 0] * d.bunch.sigma[1, 1]
                        - d.bunch.sigma[0, 1] ** 2
                    )
                )
                epsr = drdpr * 1e6 / d.bunch.p[i1]
                dzdpz = 8.0 * np.sqrt(
                    abs(
                        d.bunch.sigma[2, 2] * d.bunch.sigma[3, 3]
                        - d.bunch.sigma[2, 3] ** 2
                    )
                )
                epsz = dzdpz * 1e6 / d.bunch.p[i1]
                dp = np.sqrt(d.bunch.sigma[5, 5]) * SQRT8
                dsdp = 8.0 * np.sqrt(
                    abs(
                        d.bunch.sigma[4, 4] * d.bunch.sigma[5, 5]
                        - d.bunch.sigma[4, 5] ** 2
                    )
                )
                dwsw = dp * d.bunch.p[i1] / gamma_ref / (gamma_ref - 1.0)
                wmpa = (
                    d.bunch.rest_energy[k]
                    * (gamma_ref - 1.0)
                    / MEV
                    / d.bunch.families[k].mass_number
                )
                dwmpa = dwsw * wmpa
                dwdphi = (
                    d.rf.omega
                    / C_LIGHT
                    / d.bunch.p[i1]
                    * 180.0
                    / np.pi
                    * gamma_ref
                    * d.bunch.p[i1]
                    / gamma_ref
                    / (gamma_ref - 1.0)
                    * wmpa
                    * dsdp
                    * 1000.0
                )
                d.outputs.open['correlation'].write(
                    f' {d.run.turn:3d}  {r_ref_mm:7.2f}  {theta_ref_deg:5.1f}'
                    f' {drmm:6.2f} {dsmm:6.2f} {dphi:7.2f}'
                    f' {epsr:6.2f} {dwsw:12.7e} {wmpa:6.3f}'
                    f' {dwmpa:12.7e} {dwdphi:12.7e}'
                    f'  {dzmm:6.2f} {epsz:6.2f}\n'
                )

    # One step before the end: save the bunch positions and the
    # reverse-tracking start point
    if d.run.step == (d.run.n_steps - 1):
        d.outputs.end_x = d.bunch.x[0]
        d.outputs.end_y = d.bunch.y[0]
        d.outputs.end_z = d.bunch.z[0]
        d.outputs.end_r = r
        d.outputs.end_theta_deg = d.run.theta_deg
        d.outputs.end_phase_deg = rf_phase_deg
        if d.files.bunch_save != 'no':
            print('warning: saving the ion bunch')
            with open(d.files.bunch_save, 'wb') as fbunch:
                for i in range(d.bunch.n_particles):
                    fbunch.write(
                        struct.pack(
                            'ddd', d.bunch.x[i], d.bunch.y[i], d.bunch.z[i]
                        )
                    )

    if d.run.step == d.run.n_steps:
        pxf = d.bunch.px[0]
        pyf = d.bunch.py[0]
        pzf = d.bunch.pz[0]
        pf = np.sqrt(pxf**2 + pyf**2 + pzf**2)
        copzf = np.cos(pzf)
        alphaf = np.arccos(pxf / pf / copzf)
        if pyf < 0.0:
            alphaf = 2.0 * np.pi - alphaf
        tprf = (
            (np.pi / 2.0 + d.outputs.end_theta_deg * np.pi / 180.0) - alphaf
        ) * 1000.0
        wf = w_mev_per_u
        tpzf = np.arcsin(pzf / pf) * 1000.0
        if d.files.bunch_save != 'no':
            with open(d.files.bunch_save, 'ab') as fbunch:
                for i in range(d.bunch.n_particles):
                    fbunch.write(
                        struct.pack(
                            'ddd', d.bunch.px[i], d.bunch.py[i], d.bunch.pz[i]
                        )
                    )

        if d.files.xy_output != 'no' and 'xy' in d.outputs.open:
            fxy = d.outputs.open['xy']
            fxy.write('---------------------------------\n')
            fxy.write('start point for reverse tracking:\n')
            fxy.write('r0 theta0 z0 pr0 pz0 w0mev phi0 :\n')
            fxy.write(
                f'{d.outputs.end_r} {d.outputs.end_theta_deg} {d.outputs.end_z} {tprf}\n'
            )
            fxy.write(f'{tpzf} {wf} {d.outputs.end_phase_deg}\n')
            fxy.write('start point to continue tracking:\n')
            fxy.write(f'{r} {d.run.theta_deg} {d.bunch.z[0]} {tprf}\n')
            fxy.write(f'{tpzf} {wf} {rf_phase_deg}\n')
            fxy.write('input for DEDUCE_SIGMA:\n')
            fxy.write('amat a qn w0mev : \n')
            fxy.write(
                f'{d.bunch.families[0].atomic_mass} {d.bunch.families[0].mass_number} {d.bunch.families[0].charge_state} {wf}\n'
            )
            fxy.write(f'magnetic field (absolute value):\n{-d.bunch.bz[0]}\n')
            fxy.write(f'fhf harmonic:\n{d.rf.frequency} {d.rf.harmonic}\n')
