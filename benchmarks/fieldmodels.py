"""Field models on the LIONS map nodes for case H (DESIGN.md section
23). F0 bilinear and F2, F3 splines are in case_g1 (bilinear_field,
smooth_field); this module adds the node arrays and F1, the bicubic
Hermite field.
"""

import numpy as np


def node_arrays(grids, attrs, every):
    """LIONS map nodes (full turn, endpoint included, LIONS sign): r,
    theta, Bz, dBz/dr, dBz/dtheta.
    """
    n_sec, nth = int(attrs['sector_count']), grids['bz'].shape[1]
    if nth % every:
        raise ValueError(
            f'{nth} master columns per sector are not a multiple of {every}'
        )
    ntet = n_sec * nth // every + 1
    cols = (np.arange(ntet) * every) % nth
    th = attrs['theta_step_rad'] * every * np.arange(ntet)
    r = attrs['ri_m'] + attrs['rstep_m'] * np.arange(attrs['nrad'])
    return (
        r,
        th,
        *(-grids[k][:, cols] for k in ('bz', 'dbz_dr', 'dbz_dtheta')),
    )


def _basis(u):
    """Cubic Hermite basis on [0, 1] and its derivative: values (h00,
    h01), slopes (h10, h11).
    """
    u2, u3 = u * u, u * u * u
    value = (2 * u3 - 3 * u2 + 1, -2 * u3 + 3 * u2)
    slope = (u3 - 2 * u2 + u, u3 - u2)
    d_value = (6 * u2 - 6 * u, -6 * u2 + 6 * u)
    d_slope = (3 * u2 - 4 * u + 1, 3 * u2 - 2 * u)
    return value, slope, d_value, d_slope


def hermite_surface(r, th, bz, dbr, dbt):
    """S(rr, t) -> (Bz, dBz/dr, dBz/dtheta): bicubic Hermite, C1, from
    node values and first derivatives.

    The cross derivative d2Bz/drdtheta comes from second-order
    differences of dBz/dtheta in r.
    """
    dr, dth = r[1] - r[0], th[1] - th[0]
    drt = np.gradient(dbt, r, axis=0, edge_order=2)
    nr, nt = bz.shape

    def surface(rr, t):
        i = np.clip(((rr - r[0]) / dr).astype(int), 0, nr - 2)
        j = np.clip((t / dth).astype(int), 0, nt - 2)
        u, v = (rr - r[i]) / dr, (t - th[j]) / dth
        hu, su, dhu, dsu = _basis(u)
        hv, sv, dhv, dsv = _basis(v)
        s = ds_r = ds_t = 0.0
        for a in (0, 1):
            for b in (0, 1):
                f, fr, ft, frt = (
                    bz[i + a, j + b],
                    dbr[i + a, j + b],
                    dbt[i + a, j + b],
                    drt[i + a, j + b],
                )
                s = (
                    s
                    + f * hu[a] * hv[b]
                    + fr * dr * su[a] * hv[b]
                    + ft * dth * hu[a] * sv[b]
                    + frt * dr * dth * su[a] * sv[b]
                )
                ds_r = ds_r + (
                    f * dhu[a] * hv[b] / dr
                    + fr * dsu[a] * hv[b]
                    + ft * dth * dhu[a] * sv[b] / dr
                    + frt * dth * dsu[a] * sv[b]
                )
                ds_t = ds_t + (
                    f * hu[a] * dhv[b] / dth
                    + fr * dr * su[a] * dhv[b] / dth
                    + ft * hu[a] * dsv[b]
                    + frt * dr * su[a] * dsv[b]
                )
        return s, ds_r, ds_t

    return surface


def hermite_field(r, th, bz, dbr, dbt):
    """b(x, y, z) in T: F1, Bz from the Hermite surface, Br and Btheta
    from its derivatives (first order in z).
    """
    surface = hermite_surface(r, th, bz, dbr, dbt)

    def field(x, y, z):
        rr = np.hypot(x, y)
        t = np.arctan2(y, x) % (2 * np.pi)
        bz_, ds_r, ds_t = surface(rr, t)
        br, bt = z * ds_r, z * ds_t / rr
        c, s = np.cos(t), np.sin(t)
        return br * c - bt * s, br * s + bt * c, bz_

    return field
