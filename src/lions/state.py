"""Simulation state of LIONS, grouped in components.

Names and their Fortran counterparts: docs/glossary.md. Memory holds SI units,
angles in rad, momenta as beta gamma, unless the name says otherwise.
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Field and potential arrays keep the Fortran sizes: the interpolation reads
# neighbours up to i + 2 and j + 1, and the zero padding is part of the results.
CENTRAL_NX, CENTRAL_NY = 500, 300
OUTER_NX, OUTER_NY = 50, 1000
# MAP_NTHETA is above the Fortran 721 so that 1/3 deg full maps (1081 columns, PSI Ring) fit; smaller maps
# read the same values, the extra columns are padding.
MAP_NR, MAP_NTHETA = 250, 1081
MAX_BIDIMS = 10
MAX_BARY = 10800


@dataclass
class Run:
    title: str = ''
    deck_file: str = ''
    swap_bytes: str = ''
    steps_per_turn: int = 0
    n_steps: int = 0
    output_every: int = 0
    step: int = 0
    t0: float = 0.0  # s
    t: float = 0.0  # s
    dt: float = 0.0  # s
    c_dt: float = 0.0  # m
    turn: int = 0  # turn counter of particle 0
    theta_deg: float = 0.0  # azimuth of particle 0
    theta_prev_deg: float = 0.0
    theta_total_deg: float = 0.0
    design_r_final: float = 0.0  # m, read and not used for tracking
    design_turns: int = 0


@dataclass
class Files:
    bfield: str = ''
    central_potential: str = ''
    outer_potential: str = ''
    z_output: str = ''
    xy_output: str = ''
    sigma: str = ''
    correlation: str = ''
    matrix: str = ''
    dc_matrix: str = ''
    bunch_load: str = 'no'
    bunch_save: str = 'no'


@dataclass
class RF:
    frequency: float = 0.0  # Hz
    omega: float = 0.0  # rad/s
    omega_rev: float = 0.0  # rad/s
    harmonic: float = 0.0
    harmonic_sign: float = 0.0  # -1 for an even harmonic, +1 for an odd one
    v_right: float = 0.0  # V, sign as used in the field
    v_left: float = 0.0
    phase_right: float = 0.0  # rad
    phase_left: float = 0.0
    dee_half_angle: float = 0.0  # rad
    gap_offset: float = 0.0  # m


@dataclass
class FieldMap:
    bz: np.ndarray = field(
        default_factory=lambda: np.zeros((MAP_NR, MAP_NTHETA))
    )
    d2bz: np.ndarray = field(
        default_factory=lambda: np.zeros((MAP_NR, MAP_NTHETA))
    )
    dbzdx: np.ndarray = field(
        default_factory=lambda: np.zeros((MAP_NR, MAP_NTHETA))
    )
    dbzdy: np.ndarray = field(
        default_factory=lambda: np.zeros((MAP_NR, MAP_NTHETA))
    )
    nr: int = 0
    ntheta: int = 0
    rmin: float = 0.0  # m
    rmax: float = 0.0
    dr: float = 0.0
    dtheta: float = 0.0  # rad
    harmonic_n_coef: np.ndarray = field(
        default_factory=lambda: np.zeros(10, dtype=np.int32)
    )
    harmonic_coef: np.ndarray = field(
        default_factory=lambda: np.zeros((10, 10))
    )
    harmonic_phase: np.ndarray = field(default_factory=lambda: np.zeros(10))
    db_over_b: float = 0.0
    b_factor: Optional[float] = (
        None  # T: max|Bz| imposed on the map; None = map used as read
    )
    interpolation: str = 'bilinear'  # 'hermite': bicubic Hermite Bz from the arrays below (fields.prepare_hermite)
    azimuth: str = 'arccos'  # 'atan2': particle azimuth for the field lookup from atan2(y, x)
    dbzdr: Optional[np.ndarray] = None  # T/m
    dbzdtheta: Optional[np.ndarray] = None  # T/rad
    d2bzdrdtheta: Optional[np.ndarray] = None  # T/(m rad)


@dataclass
class Region:
    """Potential on a Cartesian grid (m): central region, or outer region in the gap frame."""

    nx: int
    ny: int
    pot: np.ndarray = None
    d2pot: np.ndarray = None
    xmin: float = 0.0
    xmax: float = 0.0
    ymin: float = 0.0
    ymax: float = 0.0
    dx: float = 0.0
    dy: float = 0.0
    rmax: float = 0.0  # central region only: central / outer boundary

    def __post_init__(self):
        self.pot = np.zeros((self.nx, self.ny))
        self.d2pot = np.zeros((self.nx, self.ny))


@dataclass
class Box:
    zmin: float = 0.0  # m
    zmax: float = 0.0


@dataclass
class IonFamily:
    start: int = 0  # index of the first particle
    size: int = 0
    atomic_mass: float = 0.0  # u
    mass_number: float = 0.0
    charge_state: float = 0.0
    r0: float = 0.0  # m
    theta0: float = 0.0  # rad
    z0: float = 0.0  # m
    pr0: float = 0.0  # rad
    pz0: float = 0.0  # rad
    w0: float = 0.0  # J, kinetic energy of the ion
    phi0_deg: float = 0.0


PER_PARTICLE = (
    'x',
    'y',
    'z',
    'px',
    'py',
    'pz',
    'gamma',
    'p',
    'r',
    'theta',
    'mass',
    'charge',
    'mass_c',
    'rest_energy',
    'bx',
    'by',
    'bz',
    'ex',
    'ey',
    'ez',
    'boris_b_coef',
    'boris_e_coef',
    'proj_x',
    'proj_y',
    'proj_offset',
    'proj_angle',
    'proj_phase',
)
MASKS = (
    'alive',
    'in_box',
    'in_central',
    'in_outer',
    'in_element',
    'in_window',
    'extracted',
    'crossed_grid',
)


class Bunch:
    """Per-particle arrays (SI, p = beta gamma), masks, families and sigma matrices."""

    def __init__(self, n_particles=0):
        self.n_particles = n_particles
        n = n_particles
        for name in PER_PARTICLE:
            setattr(self, name, np.zeros(n))
        for name in MASKS:
            setattr(self, name, np.zeros(n, dtype=bool))
        self.grid_reached = np.zeros(n, dtype=np.int32)
        self.ps = np.zeros((6, n))
        self.ps_cartesian = np.zeros((6, n))
        self.ps_projected = np.zeros((6, n))
        self.families = []
        self.w0_mev_per_u = 0.0
        self.n_alive = 0
        self.n_in_box = 0
        self.n_central = 0
        self.n_outer = 0
        self.n_in_element = 0
        self.correlations = 0
        self.phase_space = ''
        self.sigma = np.zeros((6, 6))
        self.sigma_dc = np.zeros((6, 6))
        self.dc_transform = np.zeros((6, 6))


@dataclass
class DeflectorElement:
    r: float  # m
    e: float  # V/m
    in_x: float
    in_y: float
    in_angle: float  # rad
    out_x: float
    out_y: float
    out_angle: float

    @property
    def in_ux(self):
        return np.cos(self.in_angle)

    @property
    def in_uy(self):
        return np.sin(self.in_angle)

    @property
    def out_ux(self):
        return np.cos(self.out_angle)

    @property
    def out_uy(self):
        return np.sin(self.out_angle)


@dataclass
class DipoleElement:
    r: float  # m
    gradient: float  # T/m
    x: float
    y: float
    in_angle: float  # rad
    out_angle: float

    @property
    def in_ux(self):
        return np.cos(self.in_angle)

    @property
    def in_uy(self):
        return np.sin(self.in_angle)

    @property
    def out_ux(self):
        return np.cos(self.out_angle)

    @property
    def out_uy(self):
        return np.sin(self.out_angle)


@dataclass
class GridElement:
    angle_deg: float
    x: float  # m
    y: float
    xc: float
    yc: float
    file: str = ''

    @property
    def length(self):
        return np.sqrt((self.x - self.xc) ** 2 + (self.y - self.yc) ** 2)

    @property
    def ux(self):
        return (self.x - self.xc) / self.length

    @property
    def uy(self):
        return (self.y - self.yc) / self.length


@dataclass
class Extraction:
    mode: float = 0.0  # deck rex: < 0 enables the elements
    cut_bunch: float = 0.0
    start_step: int = 0
    grids_on: float = 0.0
    grid_start: float = 0.0
    grid_start_step: int = 0
    grid_file_prefix: str = ''
    deflectors: list = field(default_factory=list)
    dipoles: list = field(default_factory=list)
    grids: list = field(default_factory=list)


@dataclass
class TrapCircle:
    x: float  # m
    y: float
    r: float


@dataclass
class Traps:
    circles: list = field(default_factory=list)
    last_step: int = 0  # traps active up to this step
    follow: int = 0  # 1: trapped particles follow particle 0


@dataclass
class VerticalSorter:
    rmin: float = 0.0  # m
    rmax: float = 0.0
    theta1: float = 0.0  # rad
    theta2: float = 0.0
    voltage: float = 0.0  # V
    gap: float = 0.0  # m, 0 = no vertical sorting
    phase: float = 0.0  # rad
    dv: float = 0.0
    nuz: float = 0.0
    nuz_z_cut: float = 0.0  # m
    voltage_seen: float = 0.0  # V, by particle 0


@dataclass
class Bidims:
    first: int = 0  # number of the first histogram, 0 = none
    every: int = 0
    start_step: int = 0
    centre_mode: str = ''
    ax0: np.ndarray = field(default_factory=lambda: np.zeros(MAX_BIDIMS))
    ax1: np.ndarray = field(default_factory=lambda: np.zeros(MAX_BIDIMS))
    ay0: np.ndarray = field(default_factory=lambda: np.zeros(MAX_BIDIMS))
    ay1: np.ndarray = field(default_factory=lambda: np.zeros(MAX_BIDIMS))
    bary_x: np.ndarray = field(default_factory=lambda: np.zeros(MAX_BARY))
    bary_y: np.ndarray = field(default_factory=lambda: np.zeros(MAX_BARY))
    bary_size: int = 0
    bary_count: int = 0


@dataclass
class Outputs:
    open: dict = field(default_factory=dict)  # open files, closed by cli.track
    headers_written: bool = False
    nuz: float = 0.0  # vertical tune of particle 0, written in the z output
    # Particle 0 at step n_steps - 1, for the end-of-run summary
    end_x: float = 0.0
    end_y: float = 0.0
    end_z: float = 0.0
    end_r: float = 0.0
    end_theta_deg: float = 0.0
    end_phase_deg: float = 0.0


class State:
    """All data of one simulation, by component."""

    def __init__(self):
        self.run = Run()
        self.files = Files()
        self.rf = RF()
        self.field = FieldMap()
        self.central = Region(CENTRAL_NX, CENTRAL_NY)
        self.outer = Region(OUTER_NX, OUTER_NY)
        self.box = Box()
        self.bunch = Bunch()
        self.extraction = Extraction()
        self.traps = Traps()
        self.sorter = VerticalSorter()
        self.bidims = Bidims()
        self.outputs = Outputs()
