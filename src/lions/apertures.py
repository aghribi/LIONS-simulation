"""Region and element masks: box, central and outer regions, deflectors,
dipoles, grids, traps.
"""

import numpy as np


def check_box(bunch, field, box):
    """Mask of live particles inside the field map and the z limits."""
    bunch.r = np.sqrt(bunch.x * bunch.x + bunch.y * bunch.y)
    bunch.in_box = (
        (bunch.r > field.rmin)
        & (bunch.r < (field.rmax - 0.02))
        & (bunch.z > box.zmin)
        & (bunch.z < box.zmax)
    )
    bunch.in_box = bunch.in_box & bunch.alive
    bunch.n_in_box = np.count_nonzero(bunch.in_box)


def check_central(bunch, central):
    """Mask of particles in the central region."""
    bunch.in_central = (
        bunch.in_box
        & (bunch.r < central.rmax)
        & (np.abs(bunch.y) < (central.ymax - central.dy))
    )
    bunch.n_central = np.count_nonzero(bunch.in_central)


def check_outer(bunch, outer, vecx, vecy):
    """Mask of particles in the outer region (vecx, vecy: gap frame)."""
    bunch.in_outer = (
        (vecx > outer.xmin)
        & (vecx < outer.xmax)
        & (vecy > outer.ymin)
        & (vecy < outer.ymax)
        & bunch.in_box
        & ~bunch.in_central
    )
    bunch.n_outer = np.count_nonzero(bunch.in_outer)


def check_deflector(bunch, extraction, run, i):
    """Mask of particles inside electrostatic deflector i (0-based)."""
    cond_in = (
        (bunch.x - extraction.deflectors[i].in_x)
        * extraction.deflectors[i].in_uy
        - (bunch.y - extraction.deflectors[i].in_y)
        * extraction.deflectors[i].in_ux
    ) <= 0.0
    cond_out = (
        (bunch.x - extraction.deflectors[i].out_x)
        * extraction.deflectors[i].out_uy
        - (bunch.y - extraction.deflectors[i].out_y)
        * extraction.deflectors[i].out_ux
    ) >= 0.0

    if i > 0:
        bunch.in_element = bunch.extracted & cond_in & cond_out & bunch.alive
        bunch.n_in_element = np.count_nonzero(bunch.in_element)
        return

    if extraction.start_step == 0:
        bunch.in_element = (
            (bunch.r > extraction.deflectors[i].r) & cond_in & cond_out
        )
        bunch.in_element = bunch.in_element | (
            bunch.extracted & cond_in & cond_out
        )
        bunch.extracted = bunch.extracted | bunch.in_element
    else:
        bunch.in_element = (
            (run.step >= extraction.start_step) & cond_in & cond_out
        )
        bunch.in_element = bunch.in_element | (
            bunch.extracted & cond_in & cond_out
        )
        bunch.extracted = bunch.extracted | bunch.in_element

    bunch.extracted = bunch.extracted & bunch.alive
    bunch.in_element = bunch.in_element & bunch.alive
    bunch.n_in_element = np.count_nonzero(bunch.in_element)


def check_dipole(bunch, extraction, i):
    """Mask of particles inside magnetic dipole i (0-based)."""
    cond_in = (
        (bunch.x - extraction.dipoles[i].x) * extraction.dipoles[i].in_uy
        - (bunch.y - extraction.dipoles[i].y) * extraction.dipoles[i].in_ux
    ) <= 0.0
    cond_out = (
        (bunch.x - extraction.dipoles[i].x) * extraction.dipoles[i].out_uy
        - (bunch.y - extraction.dipoles[i].y) * extraction.dipoles[i].out_ux
    ) >= 0.0
    bunch.in_element = bunch.extracted & cond_in & cond_out & bunch.alive
    bunch.n_in_element = np.count_nonzero(bunch.in_element)


def check_grid(bunch, extraction, run, i):
    """Mask of particles crossing extraction grid i (0-based)."""
    cond_cross = (
        (bunch.x - extraction.grids[i].x) * extraction.grids[i].uy
        - (bunch.y - extraction.grids[i].y) * extraction.grids[i].ux
    ) <= 0.0

    if (
        len(extraction.deflectors) == 0
        and len(extraction.dipoles) == 0
        and len(extraction.grids) == 1
    ):
        if run.step < extraction.start_step:
            bunch.crossed_grid[:] = False
            return
        bunch.crossed_grid = (
            (bunch.grid_reached == 0) & bunch.alive & cond_cross
        )
        reset_mask = (bunch.grid_reached == 1) & (
            (
                (bunch.x - extraction.grids[i].x) * extraction.grids[i].uy
                - (bunch.y - extraction.grids[i].y) * extraction.grids[i].ux
            )
            > 0.0
        )
        bunch.grid_reached[reset_mask] = 0
        return

    bunch.crossed_grid = (
        bunch.extracted & bunch.alive & (bunch.grid_reached == i) & cond_cross
    )


def apply_traps(bunch, box, traps, run):
    """Remove particles outside the z limits or inside a trap circle."""
    if len(traps.circles) == 0:
        return
    mask_z = (bunch.z < box.zmin) | (bunch.z > box.zmax)
    bunch.alive[mask_z] = False

    if run.step <= traps.last_step:
        for k in range(len(traps.circles)):
            r2 = traps.circles[k].r ** 2
            inside = (
                (bunch.x - traps.circles[k].x) ** 2
                + (bunch.y - traps.circles[k].y) ** 2
            ) <= r2
            bunch.alive[inside] = False

    bunch.n_alive = np.count_nonzero(bunch.alive)
    if traps.follow == 1:
        not_ok = ~bunch.alive
        bunch.x[not_ok] = bunch.x[0]
        bunch.y[not_ok] = bunch.y[0]
        bunch.z[not_ok] = bunch.z[0]
        bunch.px[not_ok] = bunch.px[0]
        bunch.py[not_ok] = bunch.py[0]
        bunch.pz[not_ok] = bunch.pz[0]
        bunch.gamma[not_ok] = bunch.gamma[0]
