"""Input deck, field maps and potentials; Fortran list-directed reading."""


import numpy as np

from .fortran_text import FortranFile, fortran_float, parse_floats
from .deck import Deck, read_deck
from .state import CENTRAL_NX, CENTRAL_NY, MAP_NR, MAP_NTHETA, OUTER_NX, OUTER_NY
from .constants import AMU, C_LIGHT, ELECTRON_MASS_U, E_CHARGE, MEV
from .beam import distribute_bunch
from .fields import add_harmonics


def prompt(prompt=''):
    """Print the prompt if any, then read one stripped line from stdin."""
    if prompt:
        print(prompt)
    return input().strip()


def _read_block(f, nrows, ncols):
    """nrows x ncols values in list-directed format; ValueError if the file ends early."""
    rows = []
    for _ in range(nrows):
        tokens = f.read_tokens(ncols)
        if len(tokens) < ncols:
            raise ValueError('file ends before the end of the block')
        rows.append([fortran_float(v) for v in tokens])
    return np.array(rows)


def read_bfield_map(path):
    """Median-plane field map: Bz, d2Bz, dBz/dx, dBz/dy on (r, theta), extended to 360 deg.

    Full map: ntheta columns over 360 deg. 1/8 map (ntheta = 91): 0 to 45 deg, read
    with the theta index reversed and extended by the 8-fold symmetry. Arrays are
    returned at the state size (MAP_NR, MAP_NTHETA).
    """
    with open(path) as raw:
        f = FortranFile(raw)
        title = f.readline().strip()
        vals = f.readline().split()
        rmin, rmax = fortran_float(vals[0]), fortran_float(vals[1])
        nr, ntheta = int(vals[2]), int(vals[3])
        if not (1 < nr <= MAP_NR and 1 < ntheta <= MAP_NTHETA):
            raise ValueError(f'{path}: map size nr={nr}, ntheta={ntheta} outside '
                             f'{MAP_NR} x {MAP_NTHETA}')
        bz, d2bz, dbzdx, dbzdy = (np.zeros((MAP_NR, MAP_NTHETA)) for _ in range(4))
        if ntheta != 91:
            print('full field map...')
            for arr in (bz, d2bz, dbzdx, dbzdy):
                arr[:nr, :ntheta] = _read_block(f, nr, ntheta)
            dtheta = 2.0 * np.pi / (ntheta - 1.0)
            ntheta_all = ntheta
        else:
            print('1/8 field map...')
            bz[:nr, 90::-1] = _read_block(f, nr, 91)
            d2bz[:nr, 90::-1] = _read_block(f, nr, 91)
            try:
                dbzx_read = _read_block(f, nr, 91)
                dbzy_read = _read_block(f, nr, 91)
            except (ValueError, IndexError):
                print('field map without dBz/dx, dBz/dy: Bx = By = 0, motion stays in the median plane')
                dbzx_read = dbzy_read = np.zeros((nr, 91))
            # 8-fold symmetry, as in the Fortran (1-based ranges in comments)
            for arr in (bz, d2bz):
                arr[:nr, 180:90:-1] = arr[:nr, :90]        # (181:92:-1) = (1:90)
                arr[:nr, 181:361] = arr[:nr, 1:181]        # (182:361) = (2:181)
                arr[:nr, 361:721] = arr[:nr, 1:361]        # (362:721) = (2:361)
            r2 = np.sqrt(2.0) / 2.0
            dbzdx[:nr, 90::-1] = r2 * (dbzx_read + dbzy_read)
            dbzdy[:nr, 90::-1] = r2 * (dbzx_read - dbzy_read)
            dbzdx[:nr, 180:90:-1] = dbzdy[:nr, :90]
            dbzdy[:nr, 180:90:-1] = dbzdx[:nr, :90]
            dbzdx[:nr, 181:361] = -dbzdx[:nr, 179::-1]
            dbzdy[:nr, 181:361] = dbzdy[:nr, 179::-1]
            dbzdx[:nr, 361:721] = dbzdx[:nr, 359::-1]
            dbzdy[:nr, 361:721] = -dbzdy[:nr, 359::-1]
            dtheta = np.pi / 4.0 / (ntheta - 1.0)
            ntheta_all = 721
    return {'title': title, 'rmin': rmin, 'rmax': rmax, 'nr': nr, 'ntheta': ntheta,
            'ntheta_all': ntheta_all, 'dr': (rmax - rmin) / (nr - 1.0), 'dtheta': dtheta,
            'bz': bz, 'd2bz': d2bz, 'dbzdx': dbzdx, 'dbzdy': dbzdy}


def read_central_potential(path, swap=False):
    """Central-region potential: Fortran unformatted file, grid header in mm.

    Records: (xmin, xmax, ymin, ymax as float64, nx, ny as int32), then nx rows of
    pot, then nx rows of d2pot (optional, zero when missing). swap: big-endian file.
    """
    fdtype = np.dtype('>f8' if swap else '<f8')
    idtype = np.dtype('>i4' if swap else '<i4')
    with open(path, 'rb') as f:
        f.read(4)
        header = np.frombuffer(f.read(32), dtype=fdtype)
        nx, ny = (int(v) for v in np.frombuffer(f.read(8), dtype=idtype))
        f.read(4)
        if not (1 < nx <= CENTRAL_NX and 1 < ny <= CENTRAL_NY):
            raise ValueError(f'{path}: grid nx={nx}, ny={ny} outside {CENTRAL_NX} x {CENTRAL_NY}'
                             ' (wrong byte order?)')
        rows = []
        for _ in range(2 * nx):
            f.read(4)
            raw = f.read(ny * 8)
            f.read(4)
            if len(raw) < ny * 8:
                break
            rows.append(np.frombuffer(raw, dtype=fdtype))
    if len(rows) < nx:
        raise ValueError(f'{path}: {len(rows)} potential rows, {nx} expected')
    pot = np.array(rows[:nx])
    if len(rows) == 2 * nx:
        d2pot = np.array(rows[nx:])
    else:
        print('central potential without d2pot: Ez = 0 in the central region')
        d2pot = np.zeros((nx, ny))
    xmin, xmax, ymin, ymax = (float(v) / 1000.0 for v in header)
    return {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax,
            'nx': nx, 'ny': ny, 'pot': pot, 'd2pot': d2pot}


def read_outer_potential(path):
    """Outer-region potential (text): xmin xmax ymin ymax nx ny in m, then pot and d2pot rows."""
    with open(path) as raw:
        f = FortranFile(raw)
        vals = f.read_tokens(6)
        xmin, xmax, ymin, ymax = (fortran_float(v) for v in vals[:4])
        nx, ny = int(fortran_float(vals[4])), int(fortran_float(vals[5]))
        if not (1 < nx <= OUTER_NX and 1 < ny <= OUTER_NY):
            raise ValueError(f'{path}: grid nx={nx}, ny={ny} outside {OUTER_NX} x {OUTER_NY}')
        pot = _read_block(f, nx, ny)
        d2pot = _read_block(f, nx, ny)
    return {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax,
            'nx': nx, 'ny': ny, 'pot': pot, 'd2pot': d2pot}


def load(d, deck):
    """Set up d from a deck (Deck or path): field maps, potentials, bunch and output files."""
    d.outputs.clear()
    if isinstance(deck, Deck):
        d.deck_file = ''
    else:
        d.deck_file = str(deck)
        deck = read_deck(deck)

    d.title = deck.title
    d.steps_per_turn = deck.steps_per_turn
    d.design_r_final = deck.design_r_final
    d.design_turns = deck.design_turns
    d.rf_frequency = deck.rf_frequency
    d.dee_half_angle = deck.dee_half_angle_deg
    d.harmonic = deck.harmonic
    d.db_over_b = deck.db_over_b
    d.phase_right = deck.phase_right_deg * np.pi / 180.0
    d.phase_left = deck.phase_left_deg * np.pi / 180.0
    d.v_dee_right = -deck.v_dee_right
    d.v_dee_left = -deck.v_dee_left
    for order, h in enumerate((deck.harmonic1, deck.harmonic2)):
        d.harmonic_n_coef[order] = 0 if h is None else len(h.coefs)
        if h is not None:
            d.harmonic_coef[:len(h.coefs), order] = h.coefs
            d.harmonic_phase[order] = h.phase_deg
    d.gap_offset = 0.0
    if deck.gap_offset is not None:
        print(f'skewed gap, offset yomega: {deck.gap_offset}')
        d.gap_offset = deck.gap_offset

    d.zmin, d.zmax, d.central_rmax = deck.zmin, deck.zmax, deck.central_rmax
    d.extraction = deck.extraction
    d.cut_bunch = deck.cut_bunch
    d.grids_on = deck.grids_on
    d.grid_start = deck.grid_start
    d.n_deflectors = 0
    d.n_dipoles = 0

    if d.extraction < 0.0:
        if d.cut_bunch <= 0.0:
            d.extraction_start_step = 0
            print('extraction: the bunch will be cut')
        else:
            d.extraction_start_step = int(d.cut_bunch + 0.5)
            print('extraction: the bunch will be kept whole')
        d.n_deflectors = len(deck.deflectors)
        for i, e in enumerate(deck.deflectors):
            d.deflector_r[i] = e.r
            d.deflector_e[i] = e.e
            d.deflector_in_x[i], d.deflector_in_y[i] = e.in_x, e.in_y
            d.deflector_out_x[i], d.deflector_out_y[i] = e.out_x, e.out_y
            d.deflector_in_angle[i] = e.in_angle_deg * np.pi / 180.0
            d.deflector_out_angle[i] = e.out_angle_deg * np.pi / 180.0
            d.deflector_in_ux[i] = np.cos(d.deflector_in_angle[i])
            d.deflector_in_uy[i] = np.sin(d.deflector_in_angle[i])
            d.deflector_out_ux[i] = np.cos(d.deflector_out_angle[i])
            d.deflector_out_uy[i] = np.sin(d.deflector_out_angle[i])
        d.n_dipoles = len(deck.dipoles)
        for i, m in enumerate(deck.dipoles):
            d.dipole_r[i] = m.r
            d.dipole_gradient[i] = m.gradient
            d.dipole_x[i], d.dipole_y[i] = m.x, m.y
            d.dipole_in_angle[i] = m.in_angle_deg * np.pi / 180.0
            d.dipole_out_angle[i] = m.out_angle_deg * np.pi / 180.0
            d.dipole_in_ux[i] = np.cos(d.dipole_in_angle[i])
            d.dipole_in_uy[i] = np.sin(d.dipole_in_angle[i])
            d.dipole_out_ux[i] = np.cos(d.dipole_out_angle[i])
            d.dipole_out_uy[i] = np.sin(d.dipole_out_angle[i])
        if d.grids_on > 0.0:
            print('extraction grids')
            d.grid_start_step = int(d.grid_start + 0.5)
            d.n_grids = len(deck.grids)
            for kg, g in enumerate(deck.grids):
                d.grid_angle_deg[kg] = g.angle_deg
                d.grid_x[kg] = g.x_mm / 1000.0
                d.grid_y[kg] = g.y_mm / 1000.0
                d.grid_xc[kg] = g.xc_mm / 1000.0
                d.grid_yc[kg] = g.yc_mm / 1000.0
                d.grid_length[kg] = np.sqrt((d.grid_x[kg] - d.grid_xc[kg])**2 +
                                            (d.grid_y[kg] - d.grid_yc[kg])**2)
                d.grid_ux[kg] = (d.grid_x[kg] - d.grid_xc[kg]) / d.grid_length[kg]
                d.grid_uy[kg] = (d.grid_y[kg] - d.grid_yc[kg]) / d.grid_length[kg]
            d.grid_file_prefix = deck.grid_file_prefix
            for kg in range(d.n_grids):
                d.grid_files[kg] = d.grid_file_prefix + '.' + str(kg + 1)
                fgrid = open(d.grid_files[kg], 'w')
                d.outputs[f'grid{kg}'] = fgrid
                fgrid.write('theta x y xc yc : \n')
                fgrid.write(f'{d.grid_angle_deg[kg]} {d.grid_x[kg]} {d.grid_y[kg]}'
                            f' {d.grid_xc[kg]} {d.grid_yc[kg]}\n')

    d.dt = d.harmonic / d.rf_frequency / d.steps_per_turn
    d.omega_rf = 2.0 * np.pi * d.rf_frequency
    d.omega_rev = d.omega_rf / d.harmonic
    d.c_dt = C_LIGHT * d.dt
    d.dee_half_angle = d.dee_half_angle * np.pi / 180.0
    d.cut_bunch = d.cut_bunch * np.pi / 180.0
    d.grids_on = d.grids_on * np.pi / 180.0
    harm = int(d.harmonic + 0.5)
    d.harmonic_sign = 1.0
    if (harm // 2) * 2 == harm:
        d.harmonic_sign = -1.0

    d.n_families = len(deck.families)
    d.n_particles = deck.total_particles
    if sum(fam.size for fam in deck.families) != d.n_particles:
        raise ValueError(f'{d.deck_file}: family sizes add up to '
                         f'{sum(fam.size for fam in deck.families)}, n = {d.n_particles}')
    d.allocate()
    d.extracted[:] = False
    d.grid_reached[:] = 0
    d.alive[:] = True
    d.n_alive = d.n_particles
    d.n_traps = 0
    d.trapped_follow = 1
    d.family_start[0] = 0
    d.w0_mev_per_u = deck.families[0].w0_mev_per_u
    for k, fam in enumerate(deck.families):
        d.family_size[k] = fam.size
        d.atomic_mass[k] = fam.atomic_mass
        d.mass_number[k] = fam.mass_number
        d.charge_state[k] = fam.charge_state
        d.r0[k] = fam.r0
        d.z0[k] = fam.z0
        d.phi0_deg[k] = fam.phi0_deg
        d.theta0[k] = fam.theta0_deg * (np.pi / 180.0)
        d.pr0[k] = fam.pr0_mrad * 0.001
        d.pz0[k] = fam.pz0_mrad * 0.001
        d.family_start[k + 1] = d.family_start[k] + d.family_size[k]
        i1 = d.family_start[k]
        i2 = d.family_start[k] + d.family_size[k]
        d.x[i1:i2] = d.r0[k] * np.cos(d.theta0[k])
        d.y[i1:i2] = d.r0[k] * np.sin(d.theta0[k])
        d.z[i1:i2] = d.z0[k]
        d.w0[k] = fam.w0_mev_per_u * MEV * d.mass_number[k]
        ion_mass = AMU * (d.atomic_mass[k] - ELECTRON_MASS_U * d.charge_state[k])
        d.mass[i1:i2] = ion_mass
        d.mass_c[i1:i2] = d.mass[i1:i2] * C_LIGHT
        d.rest_energy[i1:i2] = d.mass[i1:i2] * C_LIGHT**2
        d.gamma[i1:i2] = 1.0 + d.w0[k] / d.rest_energy[i1:i2]
        alpha = (np.pi / 2.0 + d.theta0[k]) - d.pr0[k]
        d.p[i1:i2] = np.sqrt((d.gamma[i1:i2] - 1.0) * (d.gamma[i1:i2] + 1.0))
        d.px[i1:i2] = d.p[i1:i2] * np.cos(alpha) * np.cos(d.pz0[k])
        d.py[i1:i2] = d.p[i1:i2] * np.sin(alpha) * np.cos(d.pz0[k])
        d.pz[i1:i2] = d.p[i1:i2] * np.sin(d.pz0[k])
        d.charge[i1:i2] = d.charge_state[k] * E_CHARGE
        d.boris_b_coef[i1:i2] = d.charge[i1:i2] * d.dt / 2.0 / d.mass[i1:i2]
        d.boris_e_coef[i1:i2] = d.charge[i1:i2] * d.dt / d.mass_c[i1:i2]

    if d.grids_on > 0.0 and d.n_grids > 0:
        for kg in range(d.n_grids):
            f = d.outputs[f'grid{kg}']
            f.write('particles (n A amat q):\n')
            f.write(f'{d.n_particles} {d.mass_number[0]} {d.atomic_mass[0]} {d.charge_state[0]}\n')
            f.write('RF (steps per turn, frequency, harmonic):\n')
            f.write(f'{d.steps_per_turn} {d.rf_frequency} {d.harmonic}\n')

    d.correlations = deck.correlations
    d.sigma_file = "no correlation file"
    sigma_matrices = []
    if d.correlations != 0:
        d.sigma_file = deck.sigma_file
        d.phase_space = deck.phase_space
        d.correlation_file = deck.correlation_file
        d.matrix_file = deck.matrix_file
        d.dc_matrix_file = deck.dc_matrix_file
        if d.correlation_file != 'no':
            d.outputs['correlation'] = open(d.correlation_file, 'w')
            d.outputs['correlation'].write(
                'turn     r     theta    dr    ds      dphi'
                '   epsr    dwsw        wmpa     dw          dwdphi'
                '       dz     epsz\n')
        if d.matrix_file != 'no':
            d.outputs['matrix'] = open(d.matrix_file, 'w')
        if d.dc_matrix_file != 'no':
            d.outputs['dc_matrix'] = open(d.dc_matrix_file, 'w')
            d.outputs['dc_matrix'].write(
                'turn     r     theta   drdc  dsdc  dphidc'
                ' epsrdc  dwswdc        wmpa   dwdc        dwdphidc\n')
        # Sigma matrices are read now, the bunch is built after the field map is loaded
        with open(d.sigma_file) as fsigma:
            for k in range(d.n_families):
                sigma_matrices.append(np.array([[fortran_float(v) for v in fsigma.readline().split()[:6]]
                                                for _ in range(6)]))

    d.t0 = d.phi0_deg[0] * np.pi / 180.0 / d.omega_rf
    d.t = d.t0
    d.n_steps = deck.n_steps
    d.output_every = deck.output_every
    d.z_output_file = deck.z_output_file
    d.xy_output_file = deck.xy_output_file
    d.bfield_file = deck.bfield_file
    d.central_potential_file = deck.central_potential_file
    d.outer_potential_file = deck.outer_potential_file
    d.b_factor = deck.b_factor
    d.bunch_load_file = deck.bunch_load_file
    d.bunch_save_file = deck.bunch_save_file

    bmap = read_bfield_map(d.bfield_file)
    d.map_rmin, d.map_rmax = bmap['rmin'], bmap['rmax']
    d.map_nr, d.map_ntheta = bmap['nr'], bmap['ntheta']
    d.map_dtheta = bmap['dtheta']
    ntheta_all = bmap['ntheta_all']
    d.bz_map[:] = bmap['bz']
    d.d2bz_map[:] = bmap['d2bz']
    d.dbzdx_map[:] = bmap['dbzdx']
    d.dbzdy_map[:] = bmap['dbzdy']

    if d.harmonic_n_coef[0] != 0 or d.harmonic_n_coef[1] != 0:
        print('field map plus first and/or second harmonic')
        add_harmonics(d.map_rmin, d.map_rmax, d.map_nr, ntheta_all, d.bz_map, d.d2bz_map,
             d.dbzdx_map, d.dbzdy_map, d.harmonic_n_coef, d.harmonic_coef, d.harmonic_phase)

    d.map_dr = (d.map_rmax - d.map_rmin) / (d.map_nr - 1.0)

    # Without b_factor the map is used as read (as in the Fortran).
    # With b_factor it is scaled to max|Bz| = b_factor (T).
    if d.b_factor is not None:
        bz_max = np.max(np.abs(d.bz_map[:d.map_nr, :ntheta_all]))
        if bz_max == 0.0:
            raise ValueError(f'magnetic field map is zero ({d.bfield_file}): b_factor cannot be applied')
        scale = d.b_factor / bz_max
        d.bz_map *= scale
        d.d2bz_map *= scale
        d.dbzdx_map *= scale
        d.dbzdy_map *= scale
        print(f'magnetic field map: Bmax read = {bz_max:.6e} T, B factor = {d.b_factor:.6f} T')

    d.bz_map *= (1.0 + d.db_over_b)
    d.d2bz_map *= (1.0 + d.db_over_b)

    if d.central_potential_file != 'no':
        pot = read_central_potential(d.central_potential_file, d.swap_bytes == 'y')
        d.central_xmin, d.central_xmax = pot['xmin'], pot['xmax']
        d.central_ymin, d.central_ymax = pot['ymin'], pot['ymax']
        nx1, ny1 = pot['nx'], pot['ny']
        d.central_pot[:nx1, :ny1] = pot['pot']
        d.central_d2pot[:nx1, :ny1] = pot['d2pot']
    else:
        print('no central region potential file')
        vals = parse_floats(prompt('central grid: xmin1 xmax1 ymin1 ymax1 nx1 ny1:'))
        d.central_xmin, d.central_xmax, d.central_ymin, d.central_ymax = vals[:4]
        nx1 = int(vals[4])
        ny1 = int(vals[5])
        d.central_pot[:] = 0.0
        d.central_d2pot[:] = 0.0

    d.central_dx = (d.central_xmax - d.central_xmin) / (nx1 - 1.0)
    d.central_dy = (d.central_ymax - d.central_ymin) / (ny1 - 1.0)

    pot = read_outer_potential(d.outer_potential_file)
    d.outer_xmin, d.outer_xmax = pot['xmin'], pot['xmax']
    d.outer_ymin, d.outer_ymax = pot['ymin'], pot['ymax']
    nx2, ny2 = pot['nx'], pot['ny']
    d.outer_pot[:nx2, :ny2] = pot['pot']
    d.outer_d2pot[:nx2, :ny2] = pot['d2pot']

    d.outer_dx = (d.outer_xmax - d.outer_xmin) / (nx2 - 1.0)
    d.outer_dy = (d.outer_ymax - d.outer_ymin) / (ny2 - 1.0)

    distribute_bunch(d, sigma_matrices)

    r = np.sqrt(d.x[0]**2 + d.y[0]**2)
    d.theta_deg = np.arccos(d.x[0] / r)
    if d.y[0] < 0.0:
        d.theta_deg = 2.0 * np.pi - d.theta_deg
    d.theta_deg = 180.0 / np.pi * d.theta_deg
    d.turn = 0
    if d.steps_per_turn < 0:
        d.turn = -d.n_steps // d.steps_per_turn

    if d.z_output_file != 'no':
        d.outputs['z'] = open(d.z_output_file, 'w')
    if d.xy_output_file != 'no':
        d.outputs['xy'] = open(d.xy_output_file, 'w')
