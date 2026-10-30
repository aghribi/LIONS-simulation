"""Input deck, field maps and potentials; Fortran list-directed reading."""


import numpy as np

from .fortran_text import FortranFile, fortran_float, parse_floats
from .deck import Deck, read_deck
from .state import (Bunch, DeflectorElement, DipoleElement, GridElement, IonFamily,
                    CENTRAL_NX, CENTRAL_NY, MAP_NR, MAP_NTHETA, OUTER_NX, OUTER_NY)
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
    d.outputs.open.clear()
    if isinstance(deck, Deck):
        d.run.deck_file = ''
    else:
        d.run.deck_file = str(deck)
        deck = read_deck(deck)

    run, rf, field, ext = d.run, d.rf, d.field, d.extraction
    run.title = deck.title
    run.steps_per_turn = deck.steps_per_turn
    run.design_r_final = deck.design_r_final
    run.design_turns = deck.design_turns
    rf.frequency = deck.rf_frequency
    rf.harmonic = deck.harmonic
    rf.phase_right = deck.phase_right_deg * np.pi / 180.0
    rf.phase_left = deck.phase_left_deg * np.pi / 180.0
    rf.v_right = -deck.v_dee_right
    rf.v_left = -deck.v_dee_left
    field.db_over_b = deck.db_over_b
    for order, h in enumerate((deck.harmonic1, deck.harmonic2)):
        field.harmonic_n_coef[order] = 0 if h is None else len(h.coefs)
        if h is not None:
            field.harmonic_coef[:len(h.coefs), order] = h.coefs
            field.harmonic_phase[order] = h.phase_deg
    rf.gap_offset = 0.0
    if deck.gap_offset is not None:
        print(f'skewed gap, offset yomega: {deck.gap_offset}')
        rf.gap_offset = deck.gap_offset

    d.box.zmin, d.box.zmax = deck.zmin, deck.zmax
    d.central.rmax = deck.central_rmax
    ext.mode = deck.extraction
    ext.cut_bunch = deck.cut_bunch
    ext.grids_on = deck.grids_on
    ext.grid_start = deck.grid_start
    ext.deflectors, ext.dipoles, ext.grids = [], [], []

    if ext.mode < 0.0:
        if ext.cut_bunch <= 0.0:
            ext.start_step = 0
            print('extraction: the bunch will be cut')
        else:
            ext.start_step = int(ext.cut_bunch + 0.5)
            print('extraction: the bunch will be kept whole')
        ext.deflectors = [DeflectorElement(e.r, e.e, e.in_x, e.in_y, e.in_angle_deg * np.pi / 180.0,
                                           e.out_x, e.out_y, e.out_angle_deg * np.pi / 180.0)
                          for e in deck.deflectors]
        ext.dipoles = [DipoleElement(m.r, m.gradient, m.x, m.y, m.in_angle_deg * np.pi / 180.0,
                                     m.out_angle_deg * np.pi / 180.0) for m in deck.dipoles]
        if ext.grids_on > 0.0:
            print('extraction grids')
            ext.grid_start_step = int(ext.grid_start + 0.5)
            ext.grid_file_prefix = deck.grid_file_prefix
            ext.grids = [GridElement(g.angle_deg, g.x_mm / 1000.0, g.y_mm / 1000.0,
                                     g.xc_mm / 1000.0, g.yc_mm / 1000.0,
                                     f'{deck.grid_file_prefix}.{kg + 1}')
                         for kg, g in enumerate(deck.grids)]
            for kg, grid in enumerate(ext.grids):
                fgrid = open(grid.file, 'w')
                d.outputs.open[f'grid{kg}'] = fgrid
                fgrid.write('theta x y xc yc : \n')
                fgrid.write(f'{grid.angle_deg} {grid.x} {grid.y} {grid.xc} {grid.yc}\n')

    run.dt = rf.harmonic / rf.frequency / run.steps_per_turn
    rf.omega = 2.0 * np.pi * rf.frequency
    rf.omega_rev = rf.omega / rf.harmonic
    run.c_dt = C_LIGHT * run.dt
    rf.dee_half_angle = deck.dee_half_angle_deg * np.pi / 180.0
    ext.cut_bunch = ext.cut_bunch * np.pi / 180.0
    ext.grids_on = ext.grids_on * np.pi / 180.0
    harm = int(rf.harmonic + 0.5)
    rf.harmonic_sign = -1.0 if (harm // 2) * 2 == harm else 1.0

    n_particles = deck.total_particles
    if sum(fam.size for fam in deck.families) != n_particles:
        raise ValueError(f'{run.deck_file}: family sizes add up to '
                         f'{sum(fam.size for fam in deck.families)}, n = {n_particles}')
    b = d.bunch = Bunch(n_particles)
    b.alive[:] = True
    b.n_alive = n_particles
    d.traps.circles = []
    d.traps.follow = 1
    b.w0_mev_per_u = deck.families[0].w0_mev_per_u
    start = 0
    for fam in deck.families:
        f = IonFamily(start=start, size=fam.size, atomic_mass=fam.atomic_mass,
                      mass_number=fam.mass_number, charge_state=fam.charge_state,
                      r0=fam.r0, theta0=fam.theta0_deg * (np.pi / 180.0), z0=fam.z0,
                      pr0=fam.pr0_mrad * 0.001, pz0=fam.pz0_mrad * 0.001, phi0_deg=fam.phi0_deg)
        f.w0 = fam.w0_mev_per_u * MEV * f.mass_number
        b.families.append(f)
        i1, i2 = start, start + fam.size
        start = i2
        b.x[i1:i2] = f.r0 * np.cos(f.theta0)
        b.y[i1:i2] = f.r0 * np.sin(f.theta0)
        b.z[i1:i2] = f.z0
        b.mass[i1:i2] = AMU * (f.atomic_mass - ELECTRON_MASS_U * f.charge_state)
        b.mass_c[i1:i2] = b.mass[i1:i2] * C_LIGHT
        b.rest_energy[i1:i2] = b.mass[i1:i2] * C_LIGHT**2
        b.gamma[i1:i2] = 1.0 + f.w0 / b.rest_energy[i1:i2]
        alpha = (np.pi / 2.0 + f.theta0) - f.pr0
        b.p[i1:i2] = np.sqrt((b.gamma[i1:i2] - 1.0) * (b.gamma[i1:i2] + 1.0))
        b.px[i1:i2] = b.p[i1:i2] * np.cos(alpha) * np.cos(f.pz0)
        b.py[i1:i2] = b.p[i1:i2] * np.sin(alpha) * np.cos(f.pz0)
        b.pz[i1:i2] = b.p[i1:i2] * np.sin(f.pz0)
        b.charge[i1:i2] = f.charge_state * E_CHARGE
        b.boris_b_coef[i1:i2] = b.charge[i1:i2] * run.dt / 2.0 / b.mass[i1:i2]
        b.boris_e_coef[i1:i2] = b.charge[i1:i2] * run.dt / b.mass_c[i1:i2]

    f0 = b.families[0]
    for kg in range(len(ext.grids)):
        fgrid = d.outputs.open[f'grid{kg}']
        fgrid.write('particles (n A amat q):\n')
        fgrid.write(f'{n_particles} {f0.mass_number} {f0.atomic_mass} {f0.charge_state}\n')
        fgrid.write('RF (steps per turn, frequency, harmonic):\n')
        fgrid.write(f'{run.steps_per_turn} {rf.frequency} {rf.harmonic}\n')

    files = d.files
    b.correlations = deck.correlations
    files.sigma = "no correlation file"
    sigma_matrices = []
    if b.correlations != 0:
        files.sigma = deck.sigma_file
        b.phase_space = deck.phase_space
        files.correlation = deck.correlation_file
        files.matrix = deck.matrix_file
        files.dc_matrix = deck.dc_matrix_file
        if files.correlation != 'no':
            d.outputs.open['correlation'] = open(files.correlation, 'w')
            d.outputs.open['correlation'].write(
                'turn     r     theta    dr    ds      dphi'
                '   epsr    dwsw        wmpa     dw          dwdphi'
                '       dz     epsz\n')
        if files.matrix != 'no':
            d.outputs.open['matrix'] = open(files.matrix, 'w')
        if files.dc_matrix != 'no':
            d.outputs.open['dc_matrix'] = open(files.dc_matrix, 'w')
            d.outputs.open['dc_matrix'].write(
                'turn     r     theta   drdc  dsdc  dphidc'
                ' epsrdc  dwswdc        wmpa   dwdc        dwdphidc\n')
        # Sigma matrices are read now, the bunch is built after the field map is loaded
        with open(files.sigma) as fsigma:
            for _ in b.families:
                sigma_matrices.append(np.array([[fortran_float(v) for v in fsigma.readline().split()[:6]]
                                                for _ in range(6)]))

    run.t0 = f0.phi0_deg * np.pi / 180.0 / rf.omega
    run.t = run.t0
    run.n_steps = deck.n_steps
    run.output_every = deck.output_every
    files.z_output = deck.z_output_file
    files.xy_output = deck.xy_output_file
    files.bfield = deck.bfield_file
    files.central_potential = deck.central_potential_file
    files.outer_potential = deck.outer_potential_file
    field.b_factor = deck.b_factor
    files.bunch_load = deck.bunch_load_file
    files.bunch_save = deck.bunch_save_file

    bmap = read_bfield_map(d.files.bfield)
    d.field.rmin, d.field.rmax = bmap['rmin'], bmap['rmax']
    d.field.nr, d.field.ntheta = bmap['nr'], bmap['ntheta']
    d.field.dtheta = bmap['dtheta']
    ntheta_all = bmap['ntheta_all']
    d.field.bz[:] = bmap['bz']
    d.field.d2bz[:] = bmap['d2bz']
    d.field.dbzdx[:] = bmap['dbzdx']
    d.field.dbzdy[:] = bmap['dbzdy']

    if d.field.harmonic_n_coef[0] != 0 or d.field.harmonic_n_coef[1] != 0:
        print('field map plus first and/or second harmonic')
        add_harmonics(d.field.rmin, d.field.rmax, d.field.nr, ntheta_all, d.field.bz, d.field.d2bz,
             d.field.dbzdx, d.field.dbzdy, d.field.harmonic_n_coef, d.field.harmonic_coef, d.field.harmonic_phase)

    d.field.dr = (d.field.rmax - d.field.rmin) / (d.field.nr - 1.0)

    # Without b_factor the map is used as read (as in the Fortran).
    # With b_factor it is scaled to max|Bz| = b_factor (T).
    if d.field.b_factor is not None:
        bz_max = np.max(np.abs(d.field.bz[:d.field.nr, :ntheta_all]))
        if bz_max == 0.0:
            raise ValueError(f'magnetic field map is zero ({d.files.bfield}): b_factor cannot be applied')
        scale = d.field.b_factor / bz_max
        d.field.bz *= scale
        d.field.d2bz *= scale
        d.field.dbzdx *= scale
        d.field.dbzdy *= scale
        print(f'magnetic field map: Bmax read = {bz_max:.6e} T, B factor = {d.field.b_factor:.6f} T')

    d.field.bz *= (1.0 + d.field.db_over_b)
    d.field.d2bz *= (1.0 + d.field.db_over_b)

    if d.files.central_potential != 'no':
        pot = read_central_potential(d.files.central_potential, d.run.swap_bytes == 'y')
        d.central.xmin, d.central.xmax = pot['xmin'], pot['xmax']
        d.central.ymin, d.central.ymax = pot['ymin'], pot['ymax']
        nx1, ny1 = pot['nx'], pot['ny']
        d.central.pot[:nx1, :ny1] = pot['pot']
        d.central.d2pot[:nx1, :ny1] = pot['d2pot']
    else:
        print('no central region potential file')
        vals = parse_floats(prompt('central grid: xmin1 xmax1 ymin1 ymax1 nx1 ny1:'))
        d.central.xmin, d.central.xmax, d.central.ymin, d.central.ymax = vals[:4]
        nx1 = int(vals[4])
        ny1 = int(vals[5])
        d.central.pot[:] = 0.0
        d.central.d2pot[:] = 0.0

    d.central.dx = (d.central.xmax - d.central.xmin) / (nx1 - 1.0)
    d.central.dy = (d.central.ymax - d.central.ymin) / (ny1 - 1.0)

    pot = read_outer_potential(d.files.outer_potential)
    d.outer.xmin, d.outer.xmax = pot['xmin'], pot['xmax']
    d.outer.ymin, d.outer.ymax = pot['ymin'], pot['ymax']
    nx2, ny2 = pot['nx'], pot['ny']
    d.outer.pot[:nx2, :ny2] = pot['pot']
    d.outer.d2pot[:nx2, :ny2] = pot['d2pot']

    d.outer.dx = (d.outer.xmax - d.outer.xmin) / (nx2 - 1.0)
    d.outer.dy = (d.outer.ymax - d.outer.ymin) / (ny2 - 1.0)

    distribute_bunch(d, sigma_matrices)

    r = np.sqrt(d.bunch.x[0]**2 + d.bunch.y[0]**2)
    d.run.theta_deg = np.arccos(d.bunch.x[0] / r)
    if d.bunch.y[0] < 0.0:
        d.run.theta_deg = 2.0 * np.pi - d.run.theta_deg
    d.run.theta_deg = 180.0 / np.pi * d.run.theta_deg
    d.run.turn = 0
    if d.run.steps_per_turn < 0:
        d.run.turn = -d.run.n_steps // d.run.steps_per_turn

    if d.files.z_output != 'no':
        d.outputs.open['z'] = open(d.files.z_output, 'w')
    if d.files.xy_output != 'no':
        d.outputs.open['xy'] = open(d.files.xy_output, 'w')
