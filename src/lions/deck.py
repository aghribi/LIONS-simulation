"""Input deck: the LIONS/Fortran text format as a dataclass, one reader and one writer.

Values are kept in the units of the file (m, deg, mrad, MeV/u, mm for grids, V, Hz).
Fortran names of each field: docs/guide/02_inputs.md and docs/glossary.md.
"""

from dataclasses import dataclass, field
from typing import Optional

from .fortran_text import FortranFile, fortran_float


@dataclass
class FieldHarmonic:
    """Polynomial term added to the field map (first or second harmonic)."""
    coefs: list
    phase_deg: float


@dataclass
class Deflector:
    r: float              # m, radius above which the first deflector catches particles
    e: float              # V/m
    in_x: float           # m, entrance line point
    in_y: float
    out_x: float          # m, exit line point
    out_y: float
    in_angle_deg: float
    out_angle_deg: float


@dataclass
class Dipole:
    r: float              # m, reference radius
    gradient: float       # T/m
    x: float              # m, centre
    y: float
    in_angle_deg: float
    out_angle_deg: float


@dataclass
class Grid:
    angle_deg: float
    x_mm: float           # point on the grid line
    y_mm: float
    xc_mm: float          # centre of the grid line direction
    yc_mm: float


@dataclass
class Family:
    size: int
    atomic_mass: float    # u
    mass_number: float
    charge_state: float
    r0: float             # m
    theta0_deg: float
    z0: float             # m
    pr0_mrad: float
    pz0_mrad: float
    w0_mev_per_u: float
    phi0_deg: float


@dataclass
class Deck:
    title: str
    steps_per_turn: int
    design_r_final: float          # m, read and not used for tracking
    design_turns: int              # read and not used for tracking
    rf_frequency: float            # Hz
    rf_voltage: float              # V, read and not used (dee voltages below)
    dee_half_angle_deg: float      # > 0
    harmonic: float
    db_over_b: float
    phase_right_deg: float
    phase_left_deg: float
    v_dee_right: float             # V
    v_dee_left: float              # V
    zmin: float                    # m
    zmax: float                    # m
    central_rmax: float            # m
    families: list
    n_steps: int
    output_every: int
    z_output_file: str
    xy_output_file: str
    bfield_file: str
    central_potential_file: str
    outer_potential_file: str
    harmonic1: Optional[FieldHarmonic] = None
    harmonic2: Optional[FieldHarmonic] = None
    gap_offset: Optional[float] = None   # m; set = skewed gap (written as a negative dee angle)
    extraction: float = 1.0        # < 0 enables deflectors, dipoles and grids
    cut_bunch: float = 0.0         # <= 0: bunch cut at the deflector radius; else start step
    grids_on: float = 0.0          # > 0 with extraction < 0: grids are read
    grid_start: float = 0.0        # grids active from this step
    deflectors: list = field(default_factory=list)
    dipoles: list = field(default_factory=list)
    grids: list = field(default_factory=list)
    grid_file_prefix: str = ''
    n_particles: Optional[int] = None    # written only with several families
    correlations: int = 0
    sigma_file: str = 'no'
    phase_space: str = '6d'
    correlation_file: str = 'no'
    matrix_file: str = 'no'
    dc_matrix_file: str = 'no'
    b_factor: Optional[float] = None
    bunch_load_file: str = 'no'
    bunch_save_file: str = 'no'

    @property
    def total_particles(self):
        return self.n_particles if len(self.families) > 1 else self.families[0].size


def _floats(line):
    return [fortran_float(v) for v in line.split()]


def read_deck(path):
    """Read an input deck in the Fortran list-directed format."""
    with open(path) as raw:
        f = FortranFile(raw)
        title = f.readline().strip().strip('"').strip("'")
        steps_per_turn = int(f.readline().split()[0])
        design_r_final, design_turns = _floats(f.readline())[:2]
        rf_frequency, rf_voltage, dee_angle, harmonic, db_over_b = _floats(f.readline())[:5]
        phase_right, phase_left = _floats(f.readline())[:2]
        v_right, v_left = _floats(f.readline())[:2]

        harmonics = []
        for _ in range(2):
            n = int(f.readline().split()[0])
            if n:
                coefs = _floats(f.readline())[:n]
                harmonics.append(FieldHarmonic(coefs, _floats(f.readline())[0]))
            else:
                harmonics.append(None)

        gap_offset = None
        if dee_angle < 0.0:
            dee_angle = -dee_angle
            gap_offset = _floats(f.readline())[0]

        zmin, zmax, central_rmax = _floats(f.readline())[:3]
        extraction, cut_bunch, grids_on, grid_start = _floats(f.readline())[:4]
        deflectors, dipoles, grids, grid_file_prefix = [], [], [], ''
        if extraction < 0.0:
            n = int(f.readline().split()[0])
            deflectors = [Deflector(*_floats(f.readline())[:8]) for _ in range(n)]
            n = int(f.readline().split()[0])
            dipoles = [Dipole(*_floats(f.readline())[:6]) for _ in range(n)]
            if grids_on > 0.0:
                n = int(f.readline().split()[0])
                grids = [Grid(*_floats(f.readline())[:5]) for _ in range(n)]
                grid_file_prefix = f.read_tokens(1)[0]

        n_families = int(f.readline().split()[0])
        n_particles = int(f.readline().split()[0]) if n_families > 1 else None
        families = []
        for _ in range(n_families):
            size, atomic_mass, mass_number, charge_state = _floats(f.readline())[:4]
            families.append(Family(int(size), atomic_mass, mass_number, charge_state,
                                   *_floats(f.readline())[:7]))

        correlations = int(f.readline().split()[0])
        sigma_file, phase_space = 'no', '6d'
        correlation_file = matrix_file = dc_matrix_file = 'no'
        if correlations != 0:
            sigma_file = f.read_tokens(1)[0]
            phase_space = f.read_tokens(1)[0]
            correlation_file, matrix_file, dc_matrix_file = f.read_tokens(3)

        n_steps, output_every = (int(v) for v in f.readline().split()[:2])
        z_output_file, xy_output_file = f.read_tokens(2)
        bfield_file, central_potential_file, outer_potential_file = f.read_tokens(3)

        # Optional tail: [b_factor] [ficloa ficsto]. The Fortran reads only
        # ficloa ficsto; b_factor is taken only when the first token is numeric.
        tail = []
        while True:
            tok = f.read_tokens(1)
            if not tok:
                break
            tail += tok
        b_factor = None
        if tail:
            try:
                b_factor = fortran_float(tail[0])
                tail = tail[1:]
            except ValueError:
                pass
        bunch_load_file, bunch_save_file = (tail + ['no', 'no'])[:2]

    return Deck(
        title=title, steps_per_turn=steps_per_turn,
        design_r_final=design_r_final, design_turns=int(design_turns),
        rf_frequency=rf_frequency, rf_voltage=rf_voltage, dee_half_angle_deg=dee_angle,
        harmonic=harmonic, db_over_b=db_over_b,
        phase_right_deg=phase_right, phase_left_deg=phase_left,
        v_dee_right=v_right, v_dee_left=v_left,
        zmin=zmin, zmax=zmax, central_rmax=central_rmax, families=families,
        n_steps=n_steps, output_every=output_every,
        z_output_file=z_output_file, xy_output_file=xy_output_file,
        bfield_file=bfield_file, central_potential_file=central_potential_file,
        outer_potential_file=outer_potential_file,
        harmonic1=harmonics[0], harmonic2=harmonics[1], gap_offset=gap_offset,
        extraction=extraction, cut_bunch=cut_bunch, grids_on=grids_on, grid_start=grid_start,
        deflectors=deflectors, dipoles=dipoles, grids=grids, grid_file_prefix=grid_file_prefix,
        n_particles=n_particles, correlations=correlations, sigma_file=sigma_file,
        phase_space=phase_space, correlation_file=correlation_file,
        matrix_file=matrix_file, dc_matrix_file=dc_matrix_file,
        b_factor=b_factor, bunch_load_file=bunch_load_file, bunch_save_file=bunch_save_file)


def _num(v):
    """Shortest text that reads back to the same float or int."""
    return repr(v) if isinstance(v, int) else repr(float(v))


def _line(*values):
    return ' '.join(_num(v) for v in values)


def _quoted(name):
    # Quoted so that Fortran list-directed input keeps '/' and spaces in paths
    return f"'{name}'"


def write_deck(deck, path):
    """Write a deck that both read_deck and the Fortran reference read back unchanged."""
    angle = -deck.dee_half_angle_deg if deck.gap_offset is not None else deck.dee_half_angle_deg
    lines = [
        _quoted(deck.title),
        _line(deck.steps_per_turn),
        _line(deck.design_r_final, deck.design_turns),
        _line(deck.rf_frequency, deck.rf_voltage, angle, deck.harmonic, deck.db_over_b),
        _line(deck.phase_right_deg, deck.phase_left_deg),
        _line(deck.v_dee_right, deck.v_dee_left),
    ]
    for h in (deck.harmonic1, deck.harmonic2):
        if h is None:
            lines.append('0')
        else:
            lines += [_line(len(h.coefs)), _line(*h.coefs), _line(h.phase_deg)]
    if deck.gap_offset is not None:
        lines.append(_line(deck.gap_offset))
    lines.append(_line(deck.zmin, deck.zmax, deck.central_rmax))
    lines.append(_line(deck.extraction, deck.cut_bunch, deck.grids_on, deck.grid_start))
    if deck.extraction < 0.0:
        lines.append(_line(len(deck.deflectors)))
        lines += [_line(e.r, e.e, e.in_x, e.in_y, e.out_x, e.out_y, e.in_angle_deg,
                        e.out_angle_deg) for e in deck.deflectors]
        lines.append(_line(len(deck.dipoles)))
        lines += [_line(m.r, m.gradient, m.x, m.y, m.in_angle_deg, m.out_angle_deg)
                  for m in deck.dipoles]
        if deck.grids_on > 0.0:
            lines.append(_line(len(deck.grids)))
            lines += [_line(g.angle_deg, g.x_mm, g.y_mm, g.xc_mm, g.yc_mm) for g in deck.grids]
            lines.append(_quoted(deck.grid_file_prefix))
    lines.append(_line(len(deck.families)))
    if len(deck.families) > 1:
        lines.append(_line(deck.n_particles))
    for fam in deck.families:
        lines.append(_line(fam.size, fam.atomic_mass, fam.mass_number, fam.charge_state))
        lines.append(_line(fam.r0, fam.theta0_deg, fam.z0, fam.pr0_mrad, fam.pz0_mrad,
                           fam.w0_mev_per_u, fam.phi0_deg))
    lines.append(_line(deck.correlations))
    if deck.correlations != 0:
        lines += [_quoted(deck.sigma_file), deck.phase_space,
                  ' '.join(_quoted(n) for n in (deck.correlation_file, deck.matrix_file,
                                                deck.dc_matrix_file))]
    lines.append(_line(deck.n_steps, deck.output_every))
    lines.append(' '.join(_quoted(n) for n in (deck.z_output_file, deck.xy_output_file)))
    lines.append(' '.join(_quoted(n) for n in (deck.bfield_file, deck.central_potential_file,
                                               deck.outer_potential_file)))
    if deck.b_factor is not None:
        lines.append(_line(deck.b_factor))
    if (deck.bunch_load_file, deck.bunch_save_file) != ('no', 'no'):
        lines.append(' '.join(_quoted(n) for n in (deck.bunch_load_file, deck.bunch_save_file)))
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
