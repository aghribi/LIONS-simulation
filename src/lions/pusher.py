"""Relativistic Boris pusher."""


import numpy as np



def boris_push(d):
    """Avance les particules avec le Boris pusher."""
    if d.trapped_follow == 0:
        m = d.alive  # masque
        vec_m = d.boris_b_coef[m] / d.gamma[m]
        a1_m = vec_m * d.bx[m]
        a2_m = vec_m * d.by[m]
        a3_m = vec_m * d.bz[m]
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

        fx_m = d.boris_e_coef[m] * d.ex[m]
        fy_m = d.boris_e_coef[m] * d.ey[m]
        fz_m = d.boris_e_coef[m] * d.ez[m]

        v4 = vec_m * (fx_m + 2.0 * d.px[m])
        v5 = vec_m * (fy_m + 2.0 * d.py[m])
        v6 = vec_m * (fz_m + 2.0 * d.pz[m])
        d.px[m] = d.px[m] - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx_m
        d.py[m] = d.py[m] + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy_m
        d.pz[m] = d.pz[m] + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz_m

        d.gamma[m] = np.sqrt(1.0 + d.px[m]**2 + d.py[m]**2 + d.pz[m]**2)
        vec_m = d.c_dt / d.gamma[m]
        d.x[m] = d.x[m] + vec_m * d.px[m]
        d.y[m] = d.y[m] + vec_m * d.py[m]
        d.z[m] = d.z[m] + vec_m * d.pz[m]
    else:
        # Pas de masque, toutes les particules
        vec = d.boris_b_coef / d.gamma
        a1 = vec * d.bx
        a2 = vec * d.by
        a3 = vec * d.bz
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

        fx = d.boris_e_coef * d.ex
        fy = d.boris_e_coef * d.ey
        fz = d.boris_e_coef * d.ez

        v4 = vec * (fx + 2.0 * d.px)
        v5 = vec * (fy + 2.0 * d.py)
        v6 = vec * (fz + 2.0 * d.pz)
        d.px = d.px - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx
        d.py = d.py + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy
        d.pz = d.pz + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz

        d.gamma = np.sqrt(1.0 + d.px**2 + d.py**2 + d.pz**2)
        vec = d.c_dt / d.gamma
        d.x = d.x + vec * d.px
        d.y = d.y + vec * d.py
        d.z = d.z + vec * d.pz
