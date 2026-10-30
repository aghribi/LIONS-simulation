"""Region and element masks: box, central and outer regions, deflectors, dipoles, grids, traps."""


import numpy as np



def check_box(d):
    """Vérifie quelles particules sont dans la boîte de simulation."""
    d.r = np.sqrt(d.x * d.x + d.y * d.y)
    d.in_box = ((d.r > d.map_rmin) & (d.r < (d.map_rmax - 0.02)) &
             (d.z > d.zmin) & (d.z < d.zmax))
    d.in_box = d.in_box & d.alive
    d.n_in_box = np.count_nonzero(d.in_box)


def check_central(d):
    """Vérifie les particules dans la zone centrale."""
    d.in_central = d.in_box & (d.r < d.central_rmax) & (np.abs(d.y) < (d.central_ymax - d.central_dy))
    d.n_central = np.count_nonzero(d.in_central)


def check_outer(d, vecx, vecy):
    """Vérifie les particules dans la zone périphérique (vecx, vecy : repere du gap)."""
    d.in_outer = ((vecx > d.outer_xmin) & (vecx < d.outer_xmax) &
                (vecy > d.outer_ymin) & (vecy < d.outer_ymax) &
                d.in_box & ~d.in_central)
    d.n_outer = np.count_nonzero(d.in_outer)


def check_deflector(d, i):
    """Vérifie les particules dans un déflecteur électrique."""
    # i est 0-based en Python
    cond_in = (((d.x - d.deflector_in_x[i]) * d.deflector_in_uy[i] - (d.y - d.deflector_in_y[i]) * d.deflector_in_ux[i]) <= 0.0)
    cond_out = (((d.x - d.deflector_out_x[i]) * d.deflector_out_uy[i] - (d.y - d.deflector_out_y[i]) * d.deflector_out_ux[i]) >= 0.0)

    if i > 0:
        d.in_element = d.extracted & cond_in & cond_out & d.alive
        d.n_in_element = np.count_nonzero(d.in_element)
        return

    if d.extraction_start_step == 0:
        d.in_element = (d.r > d.deflector_r[i]) & cond_in & cond_out
        d.in_element = d.in_element | (d.extracted & cond_in & cond_out)
        d.extracted = d.extracted | d.in_element
    else:
        d.in_element = (d.step >= d.extraction_start_step) & cond_in & cond_out
        d.in_element = d.in_element | (d.extracted & cond_in & cond_out)
        d.extracted = d.extracted | d.in_element

    d.extracted = d.extracted & d.alive
    d.in_element = d.in_element & d.alive
    d.n_in_element = np.count_nonzero(d.in_element)


def check_dipole(d, i):
    """Vérifie les particules dans un dipôle magnétique."""
    cond_in = (((d.x - d.dipole_x[i]) * d.dipole_in_uy[i] - (d.y - d.dipole_y[i]) * d.dipole_in_ux[i]) <= 0.0)
    cond_out = (((d.x - d.dipole_x[i]) * d.dipole_out_uy[i] - (d.y - d.dipole_y[i]) * d.dipole_out_ux[i]) >= 0.0)
    d.in_element = d.extracted & cond_in & cond_out & d.alive
    d.n_in_element = np.count_nonzero(d.in_element)


def check_grid(d, i):
    """Vérifie les particules aux grilles d'extraction."""
    cond_cross = (((d.x - d.grid_x[i]) * d.grid_uy[i] -
                   (d.y - d.grid_y[i]) * d.grid_ux[i]) <= 0.0)

    if d.n_deflectors == 0 and d.n_dipoles == 0 and d.n_grids == 1:
        if d.step < d.extraction_start_step:
            d.crossed_grid[:] = False
            return
        d.crossed_grid = (d.grid_reached == 0) & d.alive & cond_cross
        reset_mask = (d.grid_reached == 1) & (
            ((d.x - d.grid_x[i]) * d.grid_uy[i] -
             (d.y - d.grid_y[i]) * d.grid_ux[i]) > 0.0)
        d.grid_reached[reset_mask] = 0
        return

    d.crossed_grid = d.extracted & d.alive & (d.grid_reached == i) & cond_cross


def apply_traps(d):
    """Gère le piégeage des ions dans les obstacles."""
    if d.n_traps == 0:
        return
    mask_z = (d.z < d.zmin) | (d.z > d.zmax)
    d.alive[mask_z] = False

    if d.step <= d.trap_last_step:
        for k in range(d.n_traps):
            rpie2 = d.trap_r[k]**2
            mask_pie = ((d.x - d.trap_x[k])**2 + (d.y - d.trap_y[k])**2) <= rpie2
            d.alive[mask_pie] = False

    d.n_alive = np.count_nonzero(d.alive)
    if d.trapped_follow == 1:
        not_ok = ~d.alive
        d.x[not_ok] = d.x[0]
        d.y[not_ok] = d.y[0]
        d.z[not_ok] = d.z[0]
        d.px[not_ok] = d.px[0]
        d.py[not_ok] = d.py[0]
        d.pz[not_ok] = d.pz[0]
        d.gamma[not_ok] = d.gamma[0]
