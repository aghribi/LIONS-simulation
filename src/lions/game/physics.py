"""Tier T0 of the trainer: closed forms for a proton in a cyclotron (Lecture chapter 1). Tested against LIONS."""

import math

from lions.constants import AMU, C_LIGHT as C, E_CHARGE, ELECTRON_MASS_U

HYDROGEN_U = 1.00782503  # atomic mass of hydrogen (u), as the LIONS test decks
PROTON_KG = AMU * (HYDROGEN_U - ELECTRON_MASS_U)
PROTON_MEV = PROTON_KG * C**2 / E_CHARGE / 1e6


def gamma(w_mev):
    return 1 + w_mev / PROTON_MEV


def momentum(w_mev):
    """Normalised momentum beta gamma."""
    g = gamma(w_mev)
    return math.sqrt(g * g - 1)


def radius(w_mev, b_tesla):
    """Orbit radius rho = p / (q B) [m] in a uniform field."""
    return momentum(w_mev) * PROTON_KG * C / (E_CHARGE * abs(b_tesla))


def revolution_frequency(w_mev, b_tesla):
    """f = q B / (2 pi gamma m) [Hz]."""
    return E_CHARGE * abs(b_tesla) / (2 * math.pi * gamma(w_mev) * PROTON_KG)


def rotation(b_tesla):
    """+1 counter-clockwise seen from +z, -1 clockwise: a positive ion turns clockwise when Bz > 0."""
    return -1 if b_tesla > 0 else 1


def isochronous_index(w_mev):
    """Field index k = (r / B) dB/dr of the isochronous field B = gamma B0."""
    return gamma(w_mev) ** 2 - 1


def qz_squared(k, n_sectors, flutter, spiral_rad):
    """Smooth approximation: Qz^2 = -k + N^2 / (N^2 - 1) F (1 + 2 tan^2 xi)."""
    return -k + n_sectors**2 / (n_sectors**2 - 1) * flutter * (
        1 + 2 * math.tan(spiral_rad) ** 2
    )


def qr_squared(k, n_sectors, flutter, spiral_rad):
    """Smooth approximation: Qr^2 = 1 + k + 3 N^2 / ((N^2 - 1)(N^2 - 4)) F (1 + tan^2 xi)."""
    n2 = n_sectors**2
    return (
        1
        + k
        + 3
        * n2
        / ((n2 - 1) * (n2 - 4))
        * flutter
        * (1 + math.tan(spiral_rad) ** 2)
    )


# Beam on matter. Water, protons: Bortfeld, Med. Phys. 24, 2024 (1997). Highland form of multiple scattering.
BRAGG_ALPHA, BRAGG_P = 0.0022, 1.77  # R [cm] = alpha E [MeV]^p
STRAGGLING_A, STRAGGLING_B = 0.012, 0.935  # sigma [cm] = a R [cm]^b
WATER_X0_CM = (
    36.08  # radiation length [UNVERIFIED: the usual value, from memory]
)
GY_PER_FLUENCE = 1.602e-10  # Gy per (MeV cm2/g) per (proton/cm2), from 1 MeV = 1.602e-13 J and 1 g = 1e-3 kg


def proton_range(w_mev):
    """Range in water [cm], Bragg-Kleeman rule."""
    return BRAGG_ALPHA * w_mev**BRAGG_P


def energy_at_range(range_cm):
    return (range_cm / BRAGG_ALPHA) ** (1 / BRAGG_P)


def straggling(range_cm):
    return STRAGGLING_A * range_cm**STRAGGLING_B


def scattering_angle(w_mev, thickness_cm):
    """Highland: rms projected angle [rad] of a proton through water, with the energy at the entrance."""
    if thickness_cm <= 0:
        return 0.0
    g = gamma(w_mev)
    beta = math.sqrt(1 - 1 / g**2)
    p_mev = beta * g * PROTON_MEV
    x = thickness_cm / WATER_X0_CM
    return 13.6 / (beta * p_mev) * math.sqrt(x) * (1 + 0.038 * math.log(x))
