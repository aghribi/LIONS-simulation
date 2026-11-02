"""Reduced model of a pass, stage B of docs/training/PHYSICAL.md: T0's turn loop reading tables from LIONS on a field map, in place of closed forms.

Same arguments and Beam as machine.run. Rules and how each is verified: PHYSICAL.md. No phase cut, stall, clamp or septum in the turn loop. The labelled definitions
(|z| aperture, r at the map edge, extraction radius) are those of the LIONS verdict. Below the lowest valid gain, gain and orbit are held or scaled: unverified.
"""
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy.interpolate import CubicSpline

from lions.constants import C_LIGHT as C
from lions.game import machine as mc

GRID = 4000
EXTRACTED, Z_APERTURE, OUT_OF_MAP, CIRCULATING = 'extracted', 'z_aperture', 'out_of_map', 'turn_limit'


class Tables:
    """Closed orbits against momentum (beta gamma), on the unscaled map. The orbit of p in a field s B is the orbit of p / s in B."""

    def __init__(self, path):
        table = pq.read_table(path)
        meta = {k.decode(): v.decode() for k, v in table.schema.metadata.items()}
        t = table.to_pydict()
        self.rest_mev, self.w_start, self.crest_deg = float(meta['rest_mev']), float(meta['w_start_mev']), float(meta['crest_deg'])
        self.r_ext, self.r_min, self.b0_nominal, self.harmonic = float(meta['r_ext_m']), float(meta['r_min_m']), float(meta['b0_nominal_t']), int(meta['harmonic'])
        lnp, lnl, lnr = np.log(t['p']), np.log(t['path_m']), np.log(t['r_section_m'])
        self.lnp_min = lnp[0] + np.log(self.r_min) - lnr[0]        # r proportional to p below the table: r = v / omega
        self.lnp = np.linspace(self.lnp_min, lnp[-1], GRID)
        low = self.lnp < lnp[0]
        self.lnl = np.where(low, lnl[0] + self.lnp - lnp[0], CubicSpline(lnp, lnl)(np.clip(self.lnp, lnp[0], None)))
        self.lnp_ext = float(np.interp(np.log(self.r_ext), lnr, lnp))
        self.qz_lnp, self.qz = lnp, np.array(t['qz'])
        self.gain_lnw, self.gain = np.log(t['w_mev']), np.array(t['gain_mev_per_kv'])

    def path(self, lnp):
        return np.exp(np.interp(lnp, self.lnp, self.lnl))


class Reduced:
    def __init__(self, path):
        self.tables = Tables(Path(path))

    def run(self, machine, ion, knobs, errors, draws, limit_factor=1):
        """One pass, as machine.run. limit_factor: the numerical turn limit, in units of T0's; a particle still circulating there counts as lost power."""
        out = self.track(machine, ion, knobs, errors, draws, limit_factor)
        return mc.extract(machine, ion, knobs, errors, out['w_last'], out['w_prev'], out['fate'] == EXTRACTED, out['w'][out['fate'] != EXTRACTED].sum(), out['turn'], out['h'],
                          out['probe'], out['phase_sum'], out['probe'])

    def fates(self, machine, ion, knobs, errors, draws):
        out = self.track(machine, ion, knobs, errors, draws)
        return out['fate'], out['w'], out['turns']

    def track(self, machine, ion, knobs, errors, draws, limit_factor=1):
        t = self.tables
        n = mc.MACROPARTICLES
        b0 = knobs['b0'] * (1 + errors.field)
        s = b0 / t.b0_nominal
        w_inj = machine.w_injection * ion.mass_number
        g0 = ion.gamma(w_inj)
        f_rf = knobs['rf_mhz'] * 1e6
        h = max(round(f_rf / ion.f_rev(w_inj, b0)), 1)
        if h != t.harmonic:
            raise ValueError(f'harmonic {h}: the gain table was measured at {t.harmonic}')
        dee_kv = knobs['dee_kv']
        phi = mc.start_phases(knobs, errors)
        z = mc.start_heights(machine, draws)
        lns = np.log(s)
        w_min = t.rest_mev * (np.sqrt(1 + (np.exp(t.lnp_min) * s)**2) - 1)     # r = r_min
        p = np.sqrt((1 + t.w_start / t.rest_mev)**2 - 1)
        action = np.interp(np.log(p) - lns, t.qz_lnp, t.qz) * p / t.path(np.log(p) - lns)

        w, w_last, w_prev = np.full(n, t.w_start), np.full(n, t.w_start), np.full(n, t.w_start)
        alive = np.ones(n, dtype=bool)
        fate = np.full(n, CIRCULATING, dtype=object)
        turns = np.zeros(n, dtype=int)
        probe, phase_sum = np.zeros(mc.PROBE_BINS), np.zeros(mc.PROBE_BINS)
        limit = int(np.ceil(mc.turn_limit(machine, ion, knobs, b0, h))) * limit_factor
        turn = 0
        while alive.any() and turn < limit:
            turn += 1
            a = np.flatnonzero(alive)
            g = 1 + w[a] / t.rest_mev
            bg = np.sqrt(g * g - 1)
            lnq = np.log(bg) - lns
            # A crossing of theta = 0: the probes read, the fate is decided, then the turn is made.
            w_prev[a], w_last[a] = w_last[a], w[a]
            turns[a] += 1
            ga = ion.gamma(w[a])
            rho = ion.rigidity(w[a]) / (b0 * (1 + knobs['trim'] * (ga / g0 - 1)))
            i = np.minimum((rho / machine.radius * mc.PROBE_BINS).astype(int), mc.PROBE_BINS - 1)
            np.add.at(probe, i, 1.0)
            np.add.at(phase_sum, i, np.degrees(np.remainder(phi[a] + np.pi, 2 * np.pi) - np.pi))
            length = t.path(lnq)
            height = z[a] * np.sqrt(action / (np.interp(lnq, t.qz_lnp, t.qz) * bg / length))
            over_z = height > machine.z_aperture
            done = lnq >= t.lnp_ext
            fate[a[done]], fate[a[over_z]] = EXTRACTED, Z_APERTURE
            alive[a[done | over_z]] = False
            on = ~(done | over_z)
            a, phi_a, bg, g, length = a[on], phi[a[on]], bg[on], g[on], length[on]
            phi[a] = phi_a + 2 * np.pi * (f_rf * length / (bg / g * C) - h)
            after = w[a] + np.interp(np.log(w[a]), t.gain_lnw, t.gain) * dee_kv * np.cos(phi_a)
            below = after <= w_min
            w[a[~below]] = after[~below]
            alive[a[below]] = False
            fate[a[below]] = OUT_OF_MAP
        return {'fate': fate, 'w': w, 'w_last': w_last, 'w_prev': w_prev, 'turn': turn, 'turns': turns, 'h': h, 'probe': probe, 'phase_sum': phase_sum}
