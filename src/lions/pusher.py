"""Relativistic Boris pusher."""


import numpy as np



def boris_push(bunch, run, traps):
    """Advance momenta and positions by one Boris step. With trapped_follow == 0 only live particles move."""
    if traps.follow == 0:
        m = bunch.alive
        vec_m = bunch.boris_b_coef[m] / bunch.gamma[m]
        a1_m = vec_m * bunch.bx[m]
        a2_m = vec_m * bunch.by[m]
        a3_m = vec_m * bunch.bz[m]
        v1 = a1_m * a1_m
        v2 = a2_m * a2_m
        v3 = a3_m * a3_m
        vec_m = 1.0 / (1.0 + v1 + v2 + v3)
        s12 = a1_m * a2_m
        s21 = s12 - a3_m
        s12 = s12 + a3_m
        s13 = a1_m * a3_m
        s31 = s13 + a2_m
        s13 = s13 - a2_m
        s23 = a2_m * a3_m
        s32 = s23 - a1_m
        s23 = s23 + a1_m

        fx_m = bunch.boris_e_coef[m] * bunch.ex[m]
        fy_m = bunch.boris_e_coef[m] * bunch.ey[m]
        fz_m = bunch.boris_e_coef[m] * bunch.ez[m]

        v4 = vec_m * (fx_m + 2.0 * bunch.px[m])
        v5 = vec_m * (fy_m + 2.0 * bunch.py[m])
        v6 = vec_m * (fz_m + 2.0 * bunch.pz[m])
        bunch.px[m] = bunch.px[m] - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx_m
        bunch.py[m] = bunch.py[m] + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy_m
        bunch.pz[m] = bunch.pz[m] + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz_m

        bunch.gamma[m] = np.sqrt(1.0 + bunch.px[m]**2 + bunch.py[m]**2 + bunch.pz[m]**2)
        vec_m = run.c_dt / bunch.gamma[m]
        bunch.x[m] = bunch.x[m] + vec_m * bunch.px[m]
        bunch.y[m] = bunch.y[m] + vec_m * bunch.py[m]
        bunch.z[m] = bunch.z[m] + vec_m * bunch.pz[m]
    else:
        vec = bunch.boris_b_coef / bunch.gamma
        a1 = vec * bunch.bx
        a2 = vec * bunch.by
        a3 = vec * bunch.bz
        v1 = a1 * a1
        v2 = a2 * a2
        v3 = a3 * a3
        vec = 1.0 / (1.0 + v1 + v2 + v3)
        s12 = a1 * a2
        s21 = s12 - a3
        s12 = s12 + a3
        s13 = a1 * a3
        s31 = s13 + a2
        s13 = s13 - a2
        s23 = a2 * a3
        s32 = s23 - a1
        s23 = s23 + a1

        fx = bunch.boris_e_coef * bunch.ex
        fy = bunch.boris_e_coef * bunch.ey
        fz = bunch.boris_e_coef * bunch.ez

        v4 = vec * (fx + 2.0 * bunch.px)
        v5 = vec * (fy + 2.0 * bunch.py)
        v6 = vec * (fz + 2.0 * bunch.pz)
        bunch.px = bunch.px - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx
        bunch.py = bunch.py + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy
        bunch.pz = bunch.pz + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz

        bunch.gamma = np.sqrt(1.0 + bunch.px**2 + bunch.py**2 + bunch.pz**2)
        vec = run.c_dt / bunch.gamma
        bunch.x = bunch.x + vec * bunch.px
        bunch.y = bunch.y + vec * bunch.py
        bunch.z = bunch.z + vec * bunch.pz
