"""Text outputs: z and xy tracks of particle 0, correlations, end-of-run summary."""

import struct

import numpy as np

from .constants import C_LIGHT, MEV, SQRT8


def write_outputs(d):
    """Affiche les informations de suivi."""
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)

    if not d.headers_written:
        d.headers_written = True

        if d.z_output_file != 'no' and 'z' in d.outputs:
            f15 = d.outputs['z']
            f15.write('Programme LIONS\n')
            f15.write('informations selon z\n')

        if d.xy_output_file != 'no' and 'xy' in d.outputs:
            f16 = d.outputs['xy']
            f16.write('Programme LIONS\n')
            f16.write(f'{d.title}\n')
            f16.write(f' masse atomique        :    {d.atomic_mass[0]:11.7f}\n')
            f16.write(f' A                     :          {d.mass_number[0]:5.1f}\n')
            f16.write(f' q                     :          {d.charge_state[0]:5.1f}\n')
            r0mm = d.r0[0] * 1000.0
            teta0d = d.theta0[0] * 180.0 / np.pi
            z0mm = d.z0[0] * 1000.0
            pr0m = d.pr0[0] * 1000.0
            pz0m = d.pz0[0] * 1000.0
            f16.write(f' r0     : {r0mm:6.1f}\n')
            f16.write(f' teta0  : {teta0d:6.1f}\n')
            f16.write(f' z0     : {z0mm:6.1f}\n')
            f16.write(f' pr0    : {pr0m:8.2f}\n')
            f16.write(f' pz0    : {pz0m:8.2f}\n')
            f16.write(f' w0     : {d.w0_mev_per_u:14.10f}\n')
            f16.write(f' phi0   : {d.phi0_deg[0]:11.2f}\n')
            f16.write(f' fhf    : {d.rf_frequency:9.0f}\n')
            f16.write(f' harm   : {d.harmonic:9.0f}\n')
            f16.write(f' pot_d  : {-d.v_dee_right:9.1f}\n')
            f16.write(f' pot_g  : {-d.v_dee_left:9.1f}\n')
            f16.write(f' db/b   : {d.db_over_b:9.6f}\n')
            if d.b_factor is not None:
                f16.write(f' B fact : {d.b_factor:9.6f} T\n')
            dtnano = d.dt * 1e9
            f16.write(f' nb part: {d.n_particles:14d}\n')
            f16.write(f' nb iter: {d.n_steps:14d}\n')
            f16.write(f' nb i/t : {d.steps_per_turn:14d}\n')
            f16.write(f' dt(ns) : {dtnano:9.6f}\n')
            f16.write('    ite   tr   teta      x       y      xc'
                      '      yc        bz '
                      '    es      er           r      tpr'
                      '       w           dphi\n')

    d.theta_prev_deg = d.theta_deg
    rr = np.sqrt(d.x[0]**2 + d.y[0]**2)
    d.theta_deg = np.arccos(np.clip(d.x[0] / rr, -1.0, 1.0))
    if d.y[0] < 0.0:
        d.theta_deg = 2.0 * np.pi - d.theta_deg
    d.theta_deg = 180.0 / np.pi * d.theta_deg

    if d.steps_per_turn > 0:
        if d.theta_deg < d.theta_prev_deg:
            d.turn += 1
    else:
        if d.theta_deg > d.theta_prev_deg:
            d.turn -= 1

    d.theta_total_deg = d.theta_deg + 360.0 * (d.turn - 1)

    if (d.step % d.output_every != 0 and d.step != 1 and
            d.step != d.n_steps and d.step != (d.n_steps - 1)):
        return

    phihf = 180.0 / np.pi * d.omega_rf * d.t
    dph = phihf - d.harmonic * d.theta_total_deg - 90.0
    wmevpa = d.rest_energy[0] * (d.gamma[0] - 1.0) / MEV / d.mass_number[0]
    xmm = d.x[0] * 1000.0
    ymm = d.y[0] * 1000.0
    zmm = d.z[0] * 1000.0
    rmm = rr * 1000.0
    pp = np.sqrt(d.px[0]**2 + d.py[0]**2)
    ppp = np.sqrt(d.px[0]**2 + d.py[0]**2 + d.pz[0]**2)
    pr_val = (xmm * d.px[0] + ymm * d.py[0]) / rmm
    tpr = np.arcsin(pr_val / ppp) * 1000.0
    xc1mm = (d.x[0] + d.py[0] * d.dt / d.boris_e_coef[0] / d.bz[0]) * 1000.0
    yc1mm = (d.y[0] - d.px[0] * d.dt / d.boris_e_coef[0] / d.bz[0]) * 1000.0
    eskvcm = (d.ex[0] * d.px[0] + d.ey[0] * d.py[0]) / pp / 100000.0
    erkvcm = (d.ex[0] * d.py[0] - d.ey[0] * d.px[0]) / pp / 100000.0
    ezkvcm = d.ez[0] / 100000.0
    tpz = np.arcsin(d.pz[0] / ppp) * 1000.0

    if d.step % d.output_every != 0 and d.step != 1:
        pass  # skip normal output, go to end-of-run output
    else:
        if d.z_output_file != 'no' and 'z' in d.outputs:
            if abs(d.z[0]) > d.nuz_z_cut:
                anuz_c = (-1.0 / d.mass[0] / d.gamma[0] * (d.harmonic / d.omega_rf)**2 / d.z[0] *
                          d.charge[0] * (d.ez[0] + C_LIGHT / d.gamma[0] *
                          (d.px[0] * d.by[0] - d.py[0] * d.bx[0])))
                d.nuz = np.sqrt(abs(anuz_c))
            else:
                d.nuz = 0.0
            d.outputs['z'].write(
                f' {d.step:6d}  {d.turn:3d}  {d.theta_deg:5.1f}'
                f' {xmm:8.2f} {ymm:8.2f} {zmm:8.2f} {rmm:8.2f}'
                f'   {d.bx[0]:7.4f} {d.by[0]:7.4f} {d.bz[0]:7.4f}'
                f' {ezkvcm:7.3f}     {tpz:6.1f} {d.nuz:12.5f}'
                f' {d.sorter_voltage_seen:12.5f}\n')

        if d.xy_output_file != 'no' and 'xy' in d.outputs:
            d.outputs['xy'].write(
                f' {d.step:7d}  {d.turn:3d}  {d.theta_deg:5.1f}'
                f' {xmm:7.1f} {ymm:7.1f} {xc1mm:7.1f} {yc1mm:7.1f}'
                f'   {d.bz[0]:7.4f} {eskvcm:7.3f} {erkvcm:7.3f}'
                f'     {rmm:7.1f}  {tpr:6.1f} {wmevpa:9.5f} {dph:15.2f}\n')

    # Corrélations (sautees avec les sorties z/xy, comme le goto 800 du Fortran)
    if ((d.step % d.output_every == 0 or d.step == 1) and d.correlations != 0
            and not (d.correlation_file == 'no' and d.matrix_file == 'no' and d.dc_matrix_file == 'no')):
        for k in range(d.n_families):
            nbr = d.family_size[k]
            i1 = d.family_start[k]
            i2 = i1 + nbr
            d.p[i1] = np.sqrt(d.px[i1]**2 + d.py[i1]**2 + d.pz[i1]**2)
            rc = -d.mass_c[k] * d.p[i1] / d.charge[i1] / d.bz[i1]
            # Calcul phase space (simplifié, même logique que caldxs)
            d.ps_cartesian[0, :nbr] = d.x[i1:i2] - d.x[i1]
            d.ps_cartesian[4, :nbr] = d.y[i1:i2] - d.y[i1]
            ri1 = np.sqrt(d.x[i1]**2 + d.y[i1]**2)
            tetai1 = np.arccos(np.clip(d.x[i1] / ri1, -1.0, 1.0))
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

            d.sigma[:] = 0.0
            for i in range(6):
                for j in range(6):
                    d.sigma[i, j] = np.sum(d.ps[i, :nbr] * d.ps[j, :nbr] *
                                         d.alive[:nbr].astype(float)) / d.n_alive

            tetai1_deg = tetai1 * 180.0 / np.pi
            ri1mm = ri1 * 1000.0
            d.dc_transform[:] = 0.0
            d.dc_transform[0, 0] = -1.0
            d.dc_transform[1, 1] = -1.0
            for i in range(2, 6):
                d.dc_transform[i, i] = 1.0
            d.dc_transform[0, 5] = ri1 / d.p[i1]
            d.dc_transform[4, 1] = -d.dc_transform[0, 5]
            d.sigma_dc = d.dc_transform @ d.sigma @ d.dc_transform.T

            if d.correlation_file != 'no' and 'correlation' in d.outputs:
                gai1 = np.sqrt(1.0 + d.p[i1]**2)
                dr = np.sqrt(d.sigma[0, 0]) * SQRT8
                drmm = dr * 1000.0
                ds = np.sqrt(d.sigma[4, 4]) * SQRT8
                dsmm = ds * 1000.0
                dz_val = np.sqrt(d.sigma[2, 2]) * SQRT8
                dzmm = dz_val * 1000.0
                dphi = ds * d.omega_rf / C_LIGHT / d.p[i1] * 180.0 / np.pi * gai1
                drdpr = 8.0 * np.sqrt(abs(d.sigma[0, 0] * d.sigma[1, 1] - d.sigma[0, 1]**2))
                epsr = drdpr * 1e6 / d.p[i1]
                dzdpz = 8.0 * np.sqrt(abs(d.sigma[2, 2] * d.sigma[3, 3] - d.sigma[2, 3]**2))
                epsz = dzdpz * 1e6 / d.p[i1]
                dp = np.sqrt(d.sigma[5, 5]) * SQRT8
                dsdp = 8.0 * np.sqrt(abs(d.sigma[4, 4] * d.sigma[5, 5] - d.sigma[4, 5]**2))
                dwsw = dp * d.p[i1] / gai1 / (gai1 - 1.0)
                wmpa = d.rest_energy[k] * (gai1 - 1.0) / MEV / d.mass_number[k]
                dwmpa = dwsw * wmpa
                dwdphi = (d.omega_rf / C_LIGHT / d.p[i1] * 180.0 / np.pi * gai1 *
                          d.p[i1] / gai1 / (gai1 - 1.0) * wmpa * dsdp * 1000.0)
                d.outputs['correlation'].write(
                    f' {d.turn:3d}  {ri1mm:7.2f}  {tetai1_deg:5.1f}'
                    f' {drmm:6.2f} {dsmm:6.2f} {dphi:7.2f}'
                    f' {epsr:6.2f} {dwsw:12.7e} {wmpa:6.3f}'
                    f' {dwmpa:12.7e} {dwdphi:12.7e}'
                    f'  {dzmm:6.2f} {epsz:6.2f}\n')

    # Sauvegarde fin de simulation
    if d.step == (d.n_steps - 1):
        d.end_x = d.x[0]
        d.end_y = d.y[0]
        d.end_z = d.z[0]
        d.end_r = rr
        d.end_theta_deg = d.theta_deg
        d.end_phase_deg = phihf
        if d.bunch_save_file != 'no':
            print('attention : sauvegarde paquet d ions...')
            with open(d.bunch_save_file, 'wb') as f26:
                for i in range(d.n_particles):
                    f26.write(struct.pack('ddd', d.x[i], d.y[i], d.z[i]))

    if d.step == d.n_steps:
        pxf = d.px[0]
        pyf = d.py[0]
        pzf = d.pz[0]
        pf = np.sqrt(pxf**2 + pyf**2 + pzf**2)
        copzf = np.cos(pzf)
        alphaf = np.arccos(pxf / pf / copzf)
        if pyf < 0.0:
            alphaf = 2.0 * np.pi - alphaf
        tprf = ((np.pi / 2.0 + d.end_theta_deg * np.pi / 180.0) - alphaf) * 1000.0
        wf = wmevpa
        tpzf = np.arcsin(pzf / pf) * 1000.0
        if d.bunch_save_file != 'no':
            with open(d.bunch_save_file, 'ab') as f26:
                for i in range(d.n_particles):
                    f26.write(struct.pack('ddd', d.px[i], d.py[i], d.pz[i]))

        if d.xy_output_file != 'no' and 'xy' in d.outputs:
            f16 = d.outputs['xy']
            f16.write('---------------------------------\n')
            f16.write('info pour le mouvement inverse...\n')
            f16.write('r0 teta0 z0 pr0 pz0 w0mev phi0 :\n')
            f16.write(f'{d.end_r} {d.end_theta_deg} {d.end_z} {tprf}\n')
            f16.write(f'{tpzf} {wf} {d.end_phase_deg}\n')
            f16.write('info pour continuer le mouvement ...\n')
            f16.write(f'{rr} {d.theta_deg} {d.z[0]} {tprf}\n')
            f16.write(f'{tpzf} {wf} {phihf}\n')
            f16.write('info utile pour DEDUCE_SIGMA ...\n')
            f16.write('amat a qn w0mev : \n')
            f16.write(f'{d.atomic_mass[0]} {d.mass_number[0]} {d.charge_state[0]} {wf}\n')
            f16.write(f'champ magnetique (val. absolue) :\n{-d.bz[0]}\n')
            f16.write(f'fhf harmonique :\n{d.rf_frequency} {d.harmonic}\n')
