"""Input deck, field maps and potentials; Fortran list-directed reading."""

import sys

import numpy as np

from .state import CENTRAL_NX, CENTRAL_NY, MAP_NR, MAP_NTHETA, OUTER_NX, OUTER_NY
from .constants import AMU, C_LIGHT, ELECTRON_MASS_U, E_CHARGE, MEV
from .beam import distribute_bunch
from .fields import add_harmonics


def prompt(prompt=''):
    """Lit une ligne d'entrée utilisateur."""
    if prompt:
        print(prompt)
    return input().strip()


def fortran_float(s):
    """Convertit un nombre Fortran (notation d/D) en float Python."""
    return float(s.replace('d', 'e').replace('D', 'e'))


def parse_floats(s):
    """Parse une chaîne en liste de floats (supporte la notation Fortran d0)."""
    return [fortran_float(x) for x in s.split()]


def parse_ints(s):
    """Parse une chaîne en liste d'ints."""
    return [int(fortran_float(x)) for x in s.split()]


class FortranFile:
    """Émule la lecture Fortran list-directed (read(unit,*)).

    - Saute les lignes vides
    - Peut lire N tokens sur plusieurs lignes
    - Gère les chaînes entre quotes (simple ou double)
    """

    def __init__(self, f):
        self.f = f
        self._buffer = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    @staticmethod
    def _tokenize(line):
        """Découpe une ligne en tokens, en gérant les chaînes quotées."""
        tokens = []
        i = 0
        while i < len(line):
            while i < len(line) and line[i] in ' \t,':
                i += 1
            if i >= len(line):
                break
            if line[i] in "'\"":
                quote = line[i]
                j = i + 1
                while j < len(line) and line[j] != quote:
                    j += 1
                tokens.append(line[i + 1:j])
                i = j + 1
            else:
                j = i
                while j < len(line) and line[j] not in ' \t,':
                    j += 1
                tokens.append(line[i:j])
                i = j
        return tokens

    def readline(self):
        """Lit la prochaine ligne non vide (comme Fortran)."""
        while True:
            line = self.f.readline()
            if not line:
                return ''
            if line.strip():
                return line

    def read_tokens(self, n):
        """Lit exactement n tokens, en continuant sur les lignes suivantes si nécessaire."""
        result = list(self._buffer)
        self._buffer = []
        while len(result) < n:
            line = self.readline().strip()
            if not line:
                break
            result.extend(self._tokenize(line))
        if len(result) > n:
            self._buffer = result[n:]
            result = result[:n]
        return result


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


def load(d, ficin):
    """Initialisation : lecture du fichier de données ficin, cartes de champ, etc."""
    d.outputs.clear()


    d.deck_file = ficin
    with open(d.deck_file, 'r') as f20_raw:
      f20 = FortranFile(f20_raw)
      with f20:  # dummy context for indentation
        d.title = f20.readline().strip().strip('"').strip("'")
        vals = f20.readline().split()
        d.steps_per_turn = int(vals[0])
        vals = f20.readline().split()
        d.design_r_final = fortran_float(vals[0])
        d.design_turns = int(vals[1])
        vals = f20.readline().split()
        d.rf_frequency = fortran_float(vals[0])
        d.dee_half_angle = fortran_float(vals[2])
        d.harmonic = fortran_float(vals[3])
        d.db_over_b = fortran_float(vals[4])
        vals = f20.readline().split()
        d.phase_right = fortran_float(vals[0]) * np.pi / 180.0
        d.phase_left = fortran_float(vals[1]) * np.pi / 180.0
        vals = f20.readline().split()
        d.v_dee_right = -fortran_float(vals[0])
        d.v_dee_left = -fortran_float(vals[1])
        vals = f20.readline().split()
        d.harmonic_n_coef[0] = int(vals[0])
        if d.harmonic_n_coef[0] != 0:
            vals = f20.readline().split()
            for i in range(d.harmonic_n_coef[0]):
                d.harmonic_coef[i, 0] = fortran_float(vals[i])
            d.harmonic_phase[0] = fortran_float(f20.readline().split()[0])
        vals = f20.readline().split()
        d.harmonic_n_coef[1] = int(vals[0])
        if d.harmonic_n_coef[1] != 0:
            vals = f20.readline().split()
            for i in range(d.harmonic_n_coef[1]):
                d.harmonic_coef[i, 1] = fortran_float(vals[i])
            d.harmonic_phase[1] = fortran_float(f20.readline().split()[0])

        d.gap_offset = 0.0
        if d.dee_half_angle < 0.0:
            print('attention : gap en biais...')
            print('on remet teta1 positif')
            print('et on va lire le decalage...')
            d.dee_half_angle = -d.dee_half_angle
            d.gap_offset = fortran_float(f20.readline().split()[0])
            print(f'decalage yomega : {d.gap_offset}')

        vals = f20.readline().split()
        d.zmin = fortran_float(vals[0])
        d.zmax = fortran_float(vals[1])
        d.central_rmax = fortran_float(vals[2])
        vals = f20.readline().split()
        d.extraction = fortran_float(vals[0])
        d.cut_bunch = fortran_float(vals[1])
        d.grids_on = fortran_float(vals[2])
        d.grid_start = fortran_float(vals[3])
        d.n_deflectors = 0
        d.n_dipoles = 0

        if d.extraction < 0.0:
            if d.cut_bunch <= 0.0:
                d.extraction_start_step = 0
                print('extraction : le paquet sera coupe')
            else:
                d.extraction_start_step = int(d.cut_bunch + 0.5)
                print('extraction : le paquet sera entier')
            d.n_deflectors = int(f20.readline().split()[0])
            if d.n_deflectors > 0:
                for i in range(d.n_deflectors):
                    vals = f20.readline().split()
                    d.deflector_r[i] = fortran_float(vals[0])
                    d.deflector_e[i] = fortran_float(vals[1])
                    d.deflector_in_x[i] = fortran_float(vals[2])
                    d.deflector_in_y[i] = fortran_float(vals[3])
                    d.deflector_out_x[i] = fortran_float(vals[4])
                    d.deflector_out_y[i] = fortran_float(vals[5])
                    tetdfid = fortran_float(vals[6])
                    tetdfod = fortran_float(vals[7])
                    d.deflector_in_angle[i] = tetdfid * np.pi / 180.0
                    d.deflector_out_angle[i] = tetdfod * np.pi / 180.0
                    d.deflector_in_ux[i] = np.cos(d.deflector_in_angle[i])
                    d.deflector_in_uy[i] = np.sin(d.deflector_in_angle[i])
                    d.deflector_out_ux[i] = np.cos(d.deflector_out_angle[i])
                    d.deflector_out_uy[i] = np.sin(d.deflector_out_angle[i])
            d.n_dipoles = int(f20.readline().split()[0])
            if d.n_dipoles > 0:
                for i in range(d.n_dipoles):
                    vals = f20.readline().split()
                    d.dipole_r[i] = fortran_float(vals[0])
                    d.dipole_gradient[i] = fortran_float(vals[1])
                    d.dipole_x[i] = fortran_float(vals[2])
                    d.dipole_y[i] = fortran_float(vals[3])
                    tetdmid = fortran_float(vals[4])
                    tetdmod = fortran_float(vals[5])
                    d.dipole_in_angle[i] = tetdmid * np.pi / 180.0
                    d.dipole_out_angle[i] = tetdmod * np.pi / 180.0
                    d.dipole_in_ux[i] = np.cos(d.dipole_in_angle[i])
                    d.dipole_in_uy[i] = np.sin(d.dipole_in_angle[i])
                    d.dipole_out_ux[i] = np.cos(d.dipole_out_angle[i])
                    d.dipole_out_uy[i] = np.sin(d.dipole_out_angle[i])
            if d.grids_on > 0.0:
                print('grilles extraction...')
                d.grid_start_step = int(d.grid_start + 0.5)
                d.n_grids = int(f20.readline().split()[0])
                for kg in range(d.n_grids):
                    vals = f20.readline().split()
                    d.grid_angle_deg[kg] = fortran_float(vals[0])
                    d.grid_x[kg] = fortran_float(vals[1]) / 1000.0
                    d.grid_y[kg] = fortran_float(vals[2]) / 1000.0
                    d.grid_xc[kg] = fortran_float(vals[3]) / 1000.0
                    d.grid_yc[kg] = fortran_float(vals[4]) / 1000.0
                    d.grid_length[kg] = np.sqrt((d.grid_x[kg] - d.grid_xc[kg])**2 +
                                          (d.grid_y[kg] - d.grid_yc[kg])**2)
                    d.grid_ux[kg] = (d.grid_x[kg] - d.grid_xc[kg]) / d.grid_length[kg]
                    d.grid_uy[kg] = (d.grid_y[kg] - d.grid_yc[kg]) / d.grid_length[kg]
                d.grid_file_prefix = f20.readline().strip()
                for kg in range(d.n_grids):
                    ficfingex = str(kg + 1)
                    d.grid_files[kg] = d.grid_file_prefix + '.' + ficfingex
                    fgex = open(d.grid_files[kg], 'w')
                    d.outputs[f'grid{kg}'] = fgex
                    fgex.write('teta x y xc yc : \n')
                    fgex.write(f'{d.grid_angle_deg[kg]} {d.grid_x[kg]} {d.grid_y[kg]}'
                               f' {d.grid_xc[kg]} {d.grid_yc[kg]}\n')

        d.dt = d.harmonic / d.rf_frequency / d.steps_per_turn
        d.omega_rf = 2.0 * np.pi * d.rf_frequency
        d.omega_rev = d.omega_rf / d.harmonic
        d.c_dt = C_LIGHT * d.dt
        d.dee_half_angle = d.dee_half_angle * np.pi / 180.0
        d.cut_bunch = d.cut_bunch * np.pi / 180.0
        d.grids_on = d.grids_on * np.pi / 180.0
        iharm = int(d.harmonic + 0.5)
        d.harmonic_sign = 1.0
        if (iharm // 2) * 2 == iharm:
            d.harmonic_sign = -1.0

        d.n_families = int(f20.readline().split()[0])
        if d.n_families > 1:
            d.n_particles = int(f20.readline().split()[0])
        d.family_start[0] = 0  # 0-based
        n_loaded = 0

        for k in range(d.n_families):
            vals = f20.readline().split()
            d.family_size[k] = int(vals[0])
            d.atomic_mass[k] = fortran_float(vals[1])
            d.mass_number[k] = fortran_float(vals[2])
            d.charge_state[k] = fortran_float(vals[3])
            if k == 0:
                if d.n_families == 1:
                    d.n_particles = d.family_size[0]
                d.allocate()
                d.extracted[:] = False
                d.grid_reached[:] = 0
                d.alive[:] = True
                d.n_alive = d.n_particles
                d.n_traps = 0
                d.trapped_follow = 1

            vals = f20.readline().split()
            d.r0[k] = fortran_float(vals[0])
            d.theta0[k] = fortran_float(vals[1])
            d.z0[k] = fortran_float(vals[2])
            d.pr0[k] = fortran_float(vals[3])
            d.pz0[k] = fortran_float(vals[4])
            w0mev = fortran_float(vals[5])
            d.phi0_deg[k] = fortran_float(vals[6])
            if k == 0:
                d.w0_mev_per_u = w0mev

            d.theta0[k] *= np.pi / 180.0
            d.pr0[k] *= 0.001
            d.pz0[k] *= 0.001
            n_loaded += d.family_size[k]
            if n_loaded > d.n_particles:
                print(f'n np : {d.n_particles} {n_loaded}')
                print('changer n ou np...')
                sys.exit(1)
            d.family_start[k + 1] = d.family_start[k] + d.family_size[k]
            i1 = d.family_start[k]
            i2 = d.family_start[k] + d.family_size[k]
            d.x[i1:i2] = d.r0[k] * np.cos(d.theta0[k])
            d.y[i1:i2] = d.r0[k] * np.sin(d.theta0[k])
            d.z[i1:i2] = d.z0[k]
            d.w0[k] = w0mev * MEV * d.mass_number[k]
            am0k = AMU * (d.atomic_mass[k] - ELECTRON_MASS_U * d.charge_state[k])
            d.mass[i1:i2] = am0k
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
                f.write('caracteristiques particules :\n')
                f.write(f'{d.n_particles} {d.mass_number[0]} {d.atomic_mass[0]} {d.charge_state[0]}\n')
                f.write('caracteristiques HF :\n')
                f.write(f'{d.steps_per_turn} {d.rf_frequency} {d.harmonic}\n')

        d.correlations = int(f20.readline().split()[0])
        d.sigma_file = "pas de fichier correl..."
        sigma_matrices = []
        if d.correlations != 0:
            d.sigma_file = f20.readline().strip()
            d.phase_space = f20.readline().strip()
            f25 = open(d.sigma_file, 'r')
            vals = f20.read_tokens(3)
            d.correlation_file = vals[0]
            d.matrix_file = vals[1]
            d.dc_matrix_file = vals[2]
            if d.correlation_file != 'no':
                d.outputs['correlation'] = open(d.correlation_file, 'w')
                d.outputs['correlation'].write(
                    '  tr     r      teta    dr    ds      dphi'
                    '   epsr    dwsw        wmpa     dw          dwdphi'
                    '       dz     epsz\n')
            if d.matrix_file != 'no':
                d.outputs['matrix'] = open(d.matrix_file, 'w')
            if d.dc_matrix_file != 'no':
                d.outputs['dc_matrix'] = open(d.dc_matrix_file, 'w')
                d.outputs['dc_matrix'].write(
                    '  tr     r      teta   drdc  dsdc  dphidc'
                    ' epsrdc  dwswdc        wmpa   dwdc        dwdphidc\n')

            # Pre-lire les matrices sigma (calcul differe apres chargement champ magnetique)
            for k in range(d.n_families):
                cor_mat = np.zeros((6, 6), dtype=np.float64)
                for i in range(6):
                    vals_sig = f25.readline().split()
                    for j in range(6):
                        cor_mat[i, j] = fortran_float(vals_sig[j])
                sigma_matrices.append(cor_mat)
            f25.close()

        if n_loaded != d.n_particles:
            print(f'n np : {d.n_particles} {n_loaded}')
            print('changer n ou np...')
            sys.exit(1)

        d.t0 = d.phi0_deg[0] * np.pi / 180.0 / d.omega_rf
        d.t = d.t0

        vals = f20.readline().split()
        d.n_steps = int(vals[0])
        d.output_every = int(vals[1])
        vals = f20.read_tokens(2)
        d.z_output_file = vals[0]
        d.xy_output_file = vals[1]
        vals = f20.read_tokens(3)
        d.bfield_file = vals[0]
        d.central_potential_file = vals[1]
        d.outer_potential_file = vals[2]
        # Fin optionnelle : [b_factor] [ficloa ficsto]. Le Fortran ne lit que
        # ficloa ficsto ; b_factor n'est reconnu que s'il est numerique.
        tail = []
        while True:
            tok = f20.read_tokens(1)
            if not tok:
                break
            tail += tok
        d.b_factor = None
        if tail:
            try:
                d.b_factor = fortran_float(tail[0])
                tail = tail[1:]
            except ValueError:
                pass
        d.bunch_load_file, d.bunch_save_file = (tail + ['no', 'no'])[:2]

    bmap = read_bfield_map(d.bfield_file)
    d.map_rmin, d.map_rmax = bmap['rmin'], bmap['rmax']
    d.map_nr, d.map_ntheta = bmap['nr'], bmap['ntheta']
    d.map_dtheta = bmap['dtheta']
    ntetall = bmap['ntheta_all']
    d.bz_map[:] = bmap['bz']
    d.d2bz_map[:] = bmap['d2bz']
    d.dbzdx_map[:] = bmap['dbzdx']
    d.dbzdy_map[:] = bmap['dbzdy']

    if d.harmonic_n_coef[0] != 0 or d.harmonic_n_coef[1] != 0:
        print('carte mag. + rajout harmonique 1 et/ou 2...')
        add_harmonics(d.map_rmin, d.map_rmax, d.map_nr, ntetall, d.bz_map, d.d2bz_map,
             d.dbzdx_map, d.dbzdy_map, d.harmonic_n_coef, d.harmonic_coef, d.harmonic_phase)

    d.map_dr = (d.map_rmax - d.map_rmin) / (d.map_nr - 1.0)

    # Sans b_factor, la carte est utilisee telle quelle (comme le Fortran).
    # Avec b_factor, elle est ramenee a max|Bz| = b_factor (T).
    if d.b_factor is not None:
        bz_max = np.max(np.abs(d.bz_map[:d.map_nr, :ntetall]))
        if bz_max == 0.0:
            raise ValueError(f'carte magnetique nulle ({d.bfield_file}) : b_factor inapplicable')
        scale = d.b_factor / bz_max
        d.bz_map *= scale
        d.d2bz_map *= scale
        d.dbzdx_map *= scale
        d.dbzdy_map *= scale
        print(f'carte magnetique : Bmax_lu = {bz_max:.6e} T, facteur B = {d.b_factor:.6f} T')

    # Correction relative dB/B
    d.bz_map *= (1.0 + d.db_over_b)
    d.d2bz_map *= (1.0 + d.db_over_b)

    # Lecture potentiel centre
    if d.central_potential_file != 'no':
        pot = read_central_potential(d.central_potential_file, d.swap_bytes == 'y')
        d.central_xmin, d.central_xmax = pot['xmin'], pot['xmax']
        d.central_ymin, d.central_ymax = pot['ymin'], pot['ymax']
        nx1, ny1 = pot['nx'], pot['ny']
        d.central_pot[:nx1, :ny1] = pot['pot']
        d.central_d2pot[:nx1, :ny1] = pot['d2pot']
    else:
        print('pas de fichier centre cyclo...')
        vals = parse_floats(prompt('xmin1 xmax1 ymin1 ymax1 nx1 ny1 ?'))
        d.central_xmin, d.central_xmax, d.central_ymin, d.central_ymax = vals[:4]
        nx1 = int(vals[4])
        ny1 = int(vals[5])
        d.central_pot[:] = 0.0
        d.central_d2pot[:] = 0.0

    d.central_dx = (d.central_xmax - d.central_xmin) / (nx1 - 1.0)
    d.central_dy = (d.central_ymax - d.central_ymin) / (ny1 - 1.0)

    # Lecture potentiel périphérique
    pot = read_outer_potential(d.outer_potential_file)
    d.outer_xmin, d.outer_xmax = pot['xmin'], pot['xmax']
    d.outer_ymin, d.outer_ymax = pot['ymin'], pot['ymax']
    nx2, ny2 = pot['nx'], pot['ny']
    d.outer_pot[:nx2, :ny2] = pot['pot']
    d.outer_d2pot[:nx2, :ny2] = pot['d2pot']

    d.outer_dx = (d.outer_xmax - d.outer_xmin) / (nx2 - 1.0)
    d.outer_dy = (d.outer_ymax - d.outer_ymin) / (ny2 - 1.0)

    distribute_bunch(d, sigma_matrices)

    rr = np.sqrt(d.x[0]**2 + d.y[0]**2)
    d.theta_deg = np.arccos(d.x[0] / rr)
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
