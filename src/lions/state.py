"""Simulation state of LIONS (port of the Fortran module lion_m).

Names and their Fortran counterparts: docs/glossary.md. Memory holds SI units,
angles in rad, momenta as beta gamma, unless the name says otherwise.
"""

import numpy as np

# Array sizes of the Fortran module
MAX_FAMILIES = 100
CENTRAL_NX, CENTRAL_NY = 500, 300
OUTER_NX, OUTER_NY = 50, 1000
MAP_NR, MAP_NTHETA = 250, 721
MAX_DEFLECTORS = 10
MAX_DIPOLES = 10
MAX_GRIDS = 10
MAX_TRAPS = 100
MAX_BIDIMS = 10
MAX_BARY = 10800


class State:
    """All data of one simulation. Per-particle arrays are created by allocate()."""

    PER_PARTICLE = ('x', 'y', 'z', 'px', 'py', 'pz', 'gamma', 'p', 'r', 'theta',
                    'mass', 'charge', 'mass_c', 'rest_energy',
                    'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'boris_b_coef', 'boris_e_coef',
                    'proj_x', 'proj_y', 'proj_offset', 'proj_angle', 'proj_phase')
    MASKS = ('alive', 'in_box', 'in_central', 'in_outer', 'in_element', 'in_window',
             'extracted', 'crossed_grid')

    def __init__(self):
        # Run and files
        self.title = ''
        self.deck_file = ''
        self.bfield_file = ''
        self.central_potential_file = ''
        self.outer_potential_file = ''
        self.z_output_file = ''
        self.xy_output_file = ''
        self.sigma_file = ''
        self.phase_space = ''
        self.correlation_file = ''
        self.matrix_file = ''
        self.dc_matrix_file = ''
        self.bunch_load_file = 'no'
        self.bunch_save_file = 'no'
        self.outputs = {}  # open output files, closed by cli.track
        self.headers_written = False
        self.grid_file_prefix = ''
        self.grid_files = [''] * MAX_GRIDS
        self.swap_bytes = ''
        self.steps_per_turn = 0
        self.n_steps = 0
        self.output_every = 0
        self.step = 0
        self.t0 = 0.0
        self.t = 0.0
        self.dt = 0.0
        self.c_dt = 0.0
        self.turn = 0
        self.theta_deg = 0.0
        self.theta_prev_deg = 0.0
        self.theta_total_deg = 0.0
        self.design_r_final = 0.0
        self.design_turns = 0
        self.correlations = 0

        # RF and magnet
        self.rf_frequency = 0.0
        self.omega_rf = 0.0
        self.omega_rev = 0.0
        self.harmonic = 0.0
        self.harmonic_sign = 0.0
        self.v_dee_right = 0.0
        self.v_dee_left = 0.0
        self.phase_right = 0.0
        self.phase_left = 0.0
        self.dee_half_angle = 0.0
        self.gap_offset = 0.0
        self.db_over_b = 0.0
        self.b_factor = None  # T: max|Bz| imposed on the map; None = map used as read
        self.bz_map = np.zeros((MAP_NR, MAP_NTHETA))
        self.d2bz_map = np.zeros((MAP_NR, MAP_NTHETA))
        self.dbzdx_map = np.zeros((MAP_NR, MAP_NTHETA))
        self.dbzdy_map = np.zeros((MAP_NR, MAP_NTHETA))
        self.map_nr = 0
        self.map_ntheta = 0
        self.map_rmin = 0.0
        self.map_rmax = 0.0
        self.map_dr = 0.0
        self.map_dtheta = 0.0
        self.harmonic_n_coef = np.zeros(10, dtype=np.int32)
        self.harmonic_coef = np.zeros((10, 10))
        self.harmonic_phase = np.zeros(10)
        self.central_pot = np.zeros((CENTRAL_NX, CENTRAL_NY))
        self.central_d2pot = np.zeros((CENTRAL_NX, CENTRAL_NY))
        self.central_xmin = 0.0
        self.central_xmax = 0.0
        self.central_ymin = 0.0
        self.central_ymax = 0.0
        self.central_dx = 0.0
        self.central_dy = 0.0
        self.central_rmax = 0.0
        self.outer_pot = np.zeros((OUTER_NX, OUTER_NY))
        self.outer_d2pot = np.zeros((OUTER_NX, OUTER_NY))
        self.outer_xmin = 0.0
        self.outer_xmax = 0.0
        self.outer_ymin = 0.0
        self.outer_ymax = 0.0
        self.outer_dx = 0.0
        self.outer_dy = 0.0
        self.zmin = 0.0
        self.zmax = 0.0

        # Particle families
        self.n_particles = 0
        self.n_families = 0
        self.family_size = np.zeros(MAX_FAMILIES, dtype=np.int32)
        self.family_start = np.zeros(MAX_FAMILIES, dtype=np.int32)
        self.atomic_mass = np.zeros(MAX_FAMILIES)
        self.mass_number = np.zeros(MAX_FAMILIES)
        self.charge_state = np.zeros(MAX_FAMILIES)
        self.r0 = np.zeros(MAX_FAMILIES)
        self.theta0 = np.zeros(MAX_FAMILIES)
        self.z0 = np.zeros(MAX_FAMILIES)
        self.pr0 = np.zeros(MAX_FAMILIES)
        self.pz0 = np.zeros(MAX_FAMILIES)
        self.w0 = np.zeros(MAX_FAMILIES)
        self.phi0_deg = np.zeros(MAX_FAMILIES)
        self.w0_mev_per_u = 0.0
        self.n_alive = 0
        self.n_in_box = 0
        self.n_central = 0
        self.n_outer = 0
        self.n_in_element = 0

        # Correlations
        self.sigma = np.zeros((6, 6))
        self.sigma_dc = np.zeros((6, 6))
        self.dc_transform = np.zeros((6, 6))

        # Extraction
        self.extraction = 0.0
        self.cut_bunch = 0.0
        self.extraction_start_step = 0
        self.grids_on = 0.0
        self.grid_start = 0.0
        self.grid_start_step = 0
        self.n_deflectors = 0
        self.deflector_r = np.zeros(MAX_DEFLECTORS)
        self.deflector_e = np.zeros(MAX_DEFLECTORS)
        self.deflector_in_x = np.zeros(MAX_DEFLECTORS)
        self.deflector_in_y = np.zeros(MAX_DEFLECTORS)
        self.deflector_in_angle = np.zeros(MAX_DEFLECTORS)
        self.deflector_in_ux = np.zeros(MAX_DEFLECTORS)
        self.deflector_in_uy = np.zeros(MAX_DEFLECTORS)
        self.deflector_out_x = np.zeros(MAX_DEFLECTORS)
        self.deflector_out_y = np.zeros(MAX_DEFLECTORS)
        self.deflector_out_angle = np.zeros(MAX_DEFLECTORS)
        self.deflector_out_ux = np.zeros(MAX_DEFLECTORS)
        self.deflector_out_uy = np.zeros(MAX_DEFLECTORS)
        self.n_dipoles = 0
        self.dipole_r = np.zeros(MAX_DIPOLES)
        self.dipole_gradient = np.zeros(MAX_DIPOLES)
        self.dipole_x = np.zeros(MAX_DIPOLES)
        self.dipole_y = np.zeros(MAX_DIPOLES)
        self.dipole_in_angle = np.zeros(MAX_DIPOLES)
        self.dipole_in_ux = np.zeros(MAX_DIPOLES)
        self.dipole_in_uy = np.zeros(MAX_DIPOLES)
        self.dipole_out_angle = np.zeros(MAX_DIPOLES)
        self.dipole_out_ux = np.zeros(MAX_DIPOLES)
        self.dipole_out_uy = np.zeros(MAX_DIPOLES)
        self.n_grids = 0
        self.grid_angle_deg = np.zeros(MAX_GRIDS)
        self.grid_x = np.zeros(MAX_GRIDS)
        self.grid_y = np.zeros(MAX_GRIDS)
        self.grid_xc = np.zeros(MAX_GRIDS)
        self.grid_yc = np.zeros(MAX_GRIDS)
        self.grid_length = np.zeros(MAX_GRIDS)
        self.grid_ux = np.zeros(MAX_GRIDS)
        self.grid_uy = np.zeros(MAX_GRIDS)

        # Traps
        self.n_traps = 0
        self.trap_x = np.zeros(MAX_TRAPS)
        self.trap_y = np.zeros(MAX_TRAPS)
        self.trap_r = np.zeros(MAX_TRAPS)
        self.trap_last_step = 0
        self.trapped_follow = 0

        # Vertical sorter
        self.sorter_rmin = 0.0
        self.sorter_rmax = 0.0
        self.sorter_theta1 = 0.0
        self.sorter_theta2 = 0.0
        self.sorter_voltage = 0.0
        self.sorter_gap = 0.0
        self.sorter_phase = 0.0
        self.sorter_dv = 0.0
        self.sorter_nuz = 0.0
        self.sorter_voltage_seen = 0.0
        self.nuz = 0.0
        self.nuz_z_cut = 0.0

        # Bidims (2D histograms)
        self.bidim_first = 0
        self.bidim_every = 0
        self.bidim_start_step = 0
        self.centre_mode = ''
        self.bidim_ax0 = np.zeros(MAX_BIDIMS)
        self.bidim_ax1 = np.zeros(MAX_BIDIMS)
        self.bidim_ay0 = np.zeros(MAX_BIDIMS)
        self.bidim_ay1 = np.zeros(MAX_BIDIMS)
        self.bary_x = np.zeros(MAX_BARY)
        self.bary_y = np.zeros(MAX_BARY)
        self.bary_size = 0
        self.bary_count = 0

        # Particle 0 at step n_steps - 1, for the end-of-run output
        self.end_x = 0.0
        self.end_y = 0.0
        self.end_z = 0.0
        self.end_r = 0.0
        self.end_theta_deg = 0.0
        self.end_phase_deg = 0.0

        # Per-particle arrays, created by allocate()
        for name in self.PER_PARTICLE + self.MASKS:
            setattr(self, name, None)
        self.grid_reached = None
        self.ps = self.ps_cartesian = self.ps_projected = None

    def allocate(self):
        """Create the per-particle arrays for n_particles particles."""
        n = self.n_particles
        for name in self.PER_PARTICLE:
            setattr(self, name, np.zeros(n))
        for name in self.MASKS:
            setattr(self, name, np.zeros(n, dtype=bool))
        self.grid_reached = np.zeros(n, dtype=np.int32)
        self.ps = np.zeros((6, n))
        self.ps_cartesian = np.zeros((6, n))
        self.ps_projected = np.zeros((6, n))
