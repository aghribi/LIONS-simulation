"""Control-room model of a cyclotron for the missions: macroparticles from the source to the target, turn by turn.

Tier T0: closed forms along the machine, no tracked field map. It keeps what a tuner meets: frequency and field are coupled, the
harmonic sets the energy gain, phase slips, losses at the septum. It is not LIONS, and only the quantities tested in
tests/test_game.py are checked. Fixed random numbers per seed keep the readings smooth for humans and optimisers.
"""
from dataclasses import dataclass

import numpy as np

from lions.constants import C_LIGHT as C, E_CHARGE
from lions.game import physics

U_MEV = 931.494    # MeV per atomic mass unit
MACROPARTICLES = 120
PROBE_BINS = 40


@dataclass
class Ion:
    name: str
    mass_number: float
    charge: float

    @property
    def rest_mev(self):
        return physics.PROTON_MEV if self.mass_number == 1 and self.charge == 1 else self.mass_number * U_MEV

    @property
    def mass_kg(self):
        return self.rest_mev * 1e6 * E_CHARGE / C**2

    def gamma(self, w_mev):
        return 1 + w_mev / self.rest_mev

    def rigidity(self, w_mev):
        """B rho [T m]."""
        g = self.gamma(w_mev)
        return np.sqrt(g * g - 1) * self.mass_kg * C / (self.charge * E_CHARGE)

    def f_rev(self, w_mev, b):
        return self.charge * E_CHARGE * b / (2 * np.pi * self.gamma(w_mev) * self.mass_kg)


@dataclass
class Machine:
    """Fixed hardware. b0 and f_rf are knobs."""
    radius: float           # m, extraction radius
    sectors: int
    flutter: float
    spiral_deg: float
    dees: int
    dee_angle_deg: float
    w_injection: float      # MeV per nucleon
    z_injection: float      # mm, vertical amplitude at injection
    z_aperture: float       # mm, half height of the dee gap
    deflector_scale: float  # kV per MV of p v / q
    r_size: float           # mm, rms radial size at extraction
    septum: float           # mm, thickness
    deflector_gain: float   # mm of offset at the target per relative deflector error
    spot_base: float        # mm, rms spot at the best quadrupole setting
    spot_quad_width: float  # relative quadrupole error that doubles the spot area
    drift: float            # m, degrader to target
    theta_accept: float     # rad, beamline acceptance after the degrader
    space_charge: float     # Qz^2 per uA of source current at beta = 0.1: illustrative
    max_lost_kw: float      # kW of lost beam the machine tolerates


@dataclass
class Errors:
    """Hidden state of this run of the machine."""
    field: float = 0.0
    steer_x: float = 0.0    # mm at the target
    steer_y: float = 0.0
    quad: float = 0.0       # relative
    phase_deg: float = 0.0  # injection timing


@dataclass
class Beam:
    extracted_ua: float
    transmission: float
    w_mean: float           # MeV, total
    w_spread: float         # relative rms
    energies: np.ndarray    # MeV, total, extracted macroparticles
    turns: int
    lost_kw: float
    harmonic: int
    probe_r: np.ndarray
    probe: np.ndarray       # relative current against radius
    phase: np.ndarray       # deg, mean RF phase against radius
    x: float                # mm at the target
    y: float
    spot: float             # mm rms


def nominal(machine, ion, w_extraction):
    """Central field [T] and revolution frequency [Hz] of an isochronous machine extracting at w_extraction [MeV, total]."""
    b_ext = ion.rigidity(w_extraction) / machine.radius
    b0 = b_ext / ion.gamma(w_extraction) * ion.gamma(machine.w_injection * ion.mass_number)
    return b0, ion.f_rev(machine.w_injection * ion.mass_number, b0)


def machine_energy(machine, ion, b0):
    """Kinetic energy [MeV] at the extraction radius of an isochronous field: B rho = gamma B0 R gives beta c = omega0 R."""
    beta = min(ion.charge * E_CHARGE * b0 * machine.radius / (ion.mass_kg * C), 0.99)
    return (1 / np.sqrt(1 - beta * beta) - 1) * ion.rest_mev


def energy_gain(machine, ion, dee_kv, harmonic, phase):
    """MeV per turn: 2 n_dees q V sin(h alpha / 2) cos(phi)."""
    transit = np.sin(harmonic * np.radians(machine.dee_angle_deg) / 2)
    return 2 * machine.dees * ion.charge * dee_kv * 1e-3 * transit * np.cos(phase)


def deflector_nominal(machine, ion, w_mev):
    """Deflector voltage [kV], proportional to p v / q."""
    g = ion.gamma(w_mev)
    return machine.deflector_scale * (g - 1 / g) * ion.rest_mev / ion.charge


def run(machine, ion, knobs, errors, draws):
    """One pass. The harmonic is not a knob: it is the integer nearest f_rf / f_rev at injection. knobs: rf_mhz, dee_kv, phase_deg, buncher, source_ua, b0, trim, deflector_kv, steer_x, steer_y, quad, slit.

    draws: fixed standard normals (MACROPARTICLES, 3) for this seed.
    """
    phi_u = np.linspace(-np.pi, np.pi, MACROPARTICLES, endpoint=False)
    phi = phi_u - knobs['buncher'] * np.sin(phi_u) + np.radians(knobs['phase_deg'] + errors.phase_deg)     # buncher: phi - X sin(phi)
    phi = np.remainder(phi + np.pi, 2 * np.pi) - np.pi
    w0 = machine.w_injection * ion.mass_number
    w = np.full(MACROPARTICLES, w0, dtype=float)
    z = np.full(MACROPARTICLES, machine.z_injection) * np.abs(1 + 0.3 * draws[:, 0])
    alive = np.ones(MACROPARTICLES, dtype=bool)
    arrived = np.zeros(MACROPARTICLES, dtype=bool)
    b0 = knobs['b0'] * (1 + errors.field)
    g0 = ion.gamma(w0)
    f_rf = knobs['rf_mhz'] * 1e6
    h = max(round(f_rf / ion.f_rev(w0, b0)), 1)
    spiral = np.radians(machine.spiral_deg)
    focusing = machine.sectors**2 / (machine.sectors**2 - 1) * machine.flutter * (1 + 2 * np.tan(spiral)**2)

    edges = np.linspace(0, machine.radius, PROBE_BINS + 1)
    probe, phase_sum, counts = np.zeros(PROBE_BINS), np.zeros(PROBE_BINS), np.zeros(PROBE_BINS)
    lost_power = 0.0     # MeV per macroparticle, summed
    turn = 0
    w_prev = w.copy()
    full_gain = abs(energy_gain(machine, ion, knobs['dee_kv'], h, 0.0))
    max_turns = 3 * machine_energy(machine, ion, b0) / max(full_gain, 1e-6) + 50
    while alive.any() and turn < max_turns:
        turn += 1
        a = alive
        w_prev[a] = w[a]
        gained = energy_gain(machine, ion, knobs['dee_kv'], h, phi[a])
        w[a] = np.maximum(w[a] + gained, 0.5 * w0)
        stalled = w[a] <= 0.5 * w0        # decelerated back to the centre
        g = ion.gamma(w[a])
        b = b0 * (1 + knobs['trim'] * (g / g0 - 1))
        bg = np.sqrt(g * g - 1)
        rho = bg * ion.mass_kg * C / (ion.charge * E_CHARGE * b)
        phi[a] += 2 * np.pi * (f_rf * 2 * np.pi * g * ion.mass_kg / (ion.charge * E_CHARGE * b) - h)
        beta = bg / g
        k = knobs['trim'] * (g * g - 1)
        qz2 = focusing - k - machine.space_charge * knobs['source_ua'] * 0.1 / beta
        z[a] *= np.where(qz2 > 0, 1.0, np.exp(np.minimum(2 * np.pi * np.sqrt(np.abs(qz2)), 2.0)))     # defocused: grows per turn
        over_z = z[a] > machine.z_aperture
        wrong = np.cos(phi[a]) < 0.2        # behind the RF: it stops gaining and is lost on the dee structure
        i = np.minimum((rho / machine.radius * PROBE_BINS).astype(int), PROBE_BINS - 1)
        np.add.at(probe, i, 1.0)
        np.add.at(phase_sum, i, np.degrees(np.remainder(phi[a] + np.pi, 2 * np.pi) - np.pi))
        np.add.at(counts, i, 1.0)
        done = rho >= machine.radius
        lost = stalled | over_z | (wrong & ~done)
        index = np.flatnonzero(a)
        lost_power += w[index[lost]].sum()
        alive[index[lost | done]] = False
        arrived[index[done & ~lost]] = True
    lost_power += w[alive].sum()     # still circulating at the turn limit: counted as lost
    ok = arrived
    energies = w[ok]
    n_in = MACROPARTICLES
    if ok.any():
        # Septum: a particle ends the last turn anywhere within one turn separation; the blade removes its thickness.
        g_ext = ion.gamma(energies.mean())
        gain = np.abs(energies - w_prev[ok]).mean()
        separation = 1e3 * machine.radius * g_ext / (g_ext + 1) * gain / energies.mean()     # mm per turn
        survive = np.clip(1 - (machine.septum + machine.r_size) / max(separation + machine.r_size, 1e-6), 0, 1)
        lost_power += energies.sum() * (1 - survive)
    else:
        survive = 0.0
    keep = np.ones(len(energies), dtype=bool)
    if len(energies) and knobs['slit'] > 0:
        keep = np.abs(energies / energies.mean() - 1) <= knobs['slit']
    transmission = (keep.sum() / n_in) * survive
    current = knobs['source_ua'] * transmission * (1.0 if ok.any() else 0.0)
    kept = energies[keep] if len(energies) else energies
    w_mean = float(kept.mean()) if len(kept) else 0.0
    spread = float(kept.std() / w_mean) if len(kept) > 1 else 0.0
    rel = 0.0 if not len(kept) else knobs['deflector_kv'] / deflector_nominal(machine, ion, w_mean) - 1
    x = machine.deflector_gain * rel + knobs['steer_x'] + errors.steer_x
    y = knobs['steer_y'] + errors.steer_y
    q_err = (knobs['quad'] + errors.quad) / machine.spot_quad_width
    spot = machine.spot_base * float(np.sqrt(1 + q_err * q_err))
    lost_kw = lost_power / n_in * knobs['source_ua'] / ion.charge * 1e-3       # MeV x uA / q = W; to kW
    return Beam(current, transmission, w_mean, spread, kept, turn, float(lost_kw), h,
                0.5 * (edges[:-1] + edges[1:]), probe / n_in, np.divide(phase_sum, counts, out=np.zeros(PROBE_BINS), where=counts > 0),
                float(x), float(y), spot)
