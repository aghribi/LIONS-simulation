"""The tutorial helper: one next step, from what the control room reads. It never looks at the hidden errors.

These rules are also the scripted baseline of the Lab track. They are one reasonable procedure written for this model, not the CIME
operator procedure (SPEC decision D4).
"""
import math
from dataclasses import dataclass

import numpy as np


@dataclass
class Advice:
    text: str
    why: str
    knob: str = ''          # key of the knob to turn; empty: nothing to turn
    direction: int = 0
    fine: bool = False
    beam: bool = False      # hold the beam on


class Coach:
    def __init__(self, mission):
        self.m = mission
        self.keys = {s['key']: s for s in mission.specs}
        self.quad_sign, self.last_spot, self.holding = 1, None, False

    def turn(self, key, direction, text, why, fine=False):
        return Advice(text, why, key, direction, fine)

    def phase_profile(self):
        b = self.m.beam
        seen = np.flatnonzero(b.probe > 0)
        return b.phase[seen] if len(seen) > 1 else np.zeros(2)

    def rf_steps(self):
        """Get the beam out: frequency first, then the injection phase, from the phase probe."""
        phase = self.phase_profile()
        drift = phase[-1] - phase[0]
        third = max(len(phase) // 3, 1)
        extracted = self.m.has_beam()
        slope_in = phase[third] - phase[0]
        slope_out = phase[-1] - phase[-third - 1]
        if 'trim' in self.keys and extracted and abs(slope_out - slope_in) > 10:
            return self.turn('trim', 1 if slope_out > slope_in else -1, 'TURN TRIM',
                             'The phase bends with radius: the field does not follow gamma.', fine=abs(slope_out - slope_in) < 30)
        if abs(drift) > (10 if extracted else 15):
            return self.turn('rf_mhz', -1 if drift > 0 else 1, 'LOWER THE RF' if drift > 0 else 'RAISE THE RF',
                             f'The phase {"climbs" if drift > 0 else "falls"} with radius: the RF is too {"high" if drift > 0 else "low"} for this field.',
                             fine=abs(drift) < 40)
        mean = float(phase.mean())
        if 'phase_deg' in self.keys and abs(mean) > 8:
            return self.turn('phase_deg', -1 if mean > 0 else 1, 'TURN THE INJECTION PHASE',
                             'Flat phase, but off the crest: the beam gains less and is more easily lost.', fine=abs(mean) < 20)
        return None

    def energy_step(self, wanted, tolerance):
        w = self.m.beam.w_mean
        if w > 0 and abs(w / wanted - 1) > tolerance:
            up = w < wanted
            return self.turn('b0', 1 if up else -1, 'RAISE THE FIELD' if up else 'LOWER THE FIELD',
                             'The radius is fixed, so the field sets the energy. The RF must follow it.', fine=abs(w / wanted - 1) < 3 * tolerance)
        return None

    def axis_steps(self, limit):
        r = self.m.readings
        if not self.m.has_beam():
            return None
        for key, reading, label in (('steer_x', r['x_mm'], 'STEER X'), ('steer_y', r['y_mm'], 'STEER Y')):
            if key in self.keys and abs(reading) > 0.5 * limit:
                if key == 'steer_x' and abs(reading) > 12 and 'deflector_kv' in self.keys:
                    return self.turn('deflector_kv', -1 if reading > 0 else 1, 'TURN THE DEFLECTOR', 'The beam leaves far off axis: the deflector does not match its rigidity.')
                return self.turn(key, -1 if reading > 0 else 1, label, f'The beam is {abs(reading):.1f} mm off axis.', fine=abs(reading) < 2 * limit)
        return None

    def quad_hunt(self, want_bigger, text, why):
        """The spot has no sign: turn, watch, turn back if it went the wrong way."""
        spot = self.m.spot()
        if self.last_spot is not None and abs(spot - self.last_spot) > 1e-6 and (spot > self.last_spot) != want_bigger:
            self.quad_sign = -self.quad_sign
        self.last_spot = spot
        return self.turn('quad', self.quad_sign, text, why)

    def advise(self):
        ready = self.m.ready()
        if self.holding and ready and not getattr(self.m, 'release_due', lambda: False)():     # in spec: stop retuning
            return Advice('KEEP HOLDING SPACE', 'In spec. Watch the readings.', beam=True)
        self.holding = False
        step = None if ready else self.rf_steps()
        if not step:
            step = self.mission_steps() if self.m.has_beam() else self.lost_beam()
        self.holding = step.beam
        return step

    def lost_beam(self):
        reach = self.m.reach()
        if 'buncher' in self.keys and abs(self.m.knobs['buncher'] - 1) > 0.2:
            return self.turn('buncher', 1 if self.m.knobs['buncher'] < 1 else -1, 'SET THE BUNCHER NEAR 1', 'Bunched beam fits the RF acceptance.')
        if reach < 0.9 and 'b0' in self.keys:
            return self.turn('b0', 1, 'RAISE THE FIELD', 'Flat phase but the beam stops inside: check the field.')
        return Advice('WATCH THE RADIAL PROBE', 'The beam stops inside the machine.')

    def loss_steps(self):
        m = self.m
        if m.beam.lost_kw > 0.8 * m.machine.max_lost_kw:
            if 'dee_kv' in self.keys and m.knobs['dee_kv'] < self.keys['dee_kv']['max']:
                return self.turn('dee_kv', 1, 'RAISE THE DEE VOLTAGE', 'More gain per turn separates the last turns: less beam on the septum.')
            if 'source_ua' in self.keys:
                return self.turn('source_ua', -1, 'LOWER THE SOURCE', 'Too much beam lost in the machine.')
        return None

    def room_for_more(self):
        """Raising the source must not undo the loss step: a dead band between the two."""
        return self.m.beam.lost_kw < 0.65 * self.m.machine.max_lost_kw

    def mission_steps(self):
        if 'buncher' in self.keys and abs(self.m.knobs['buncher'] - 1) > 0.2:
            return self.turn('buncher', 1 if self.m.knobs['buncher'] < 1 else -1, 'SET THE BUNCHER NEAR 1', 'A bunched beam fits the RF acceptance: more current for the same source.')
        return getattr(self, f'steps_{self.m.key}')()

    def steps_therapy(self):
        m, p = self.m, self.m.p
        error = m.peak - p['depth_cm']
        if abs(error) > 0.5 * p['half_width_cm']:
            return self.turn('degrader_cm', 1 if error > 0 else -1, 'MORE DEGRADER' if error > 0 else 'LESS DEGRADER',
                             f'The Bragg peak is at {m.peak:.1f} cm, the tumour at {p["depth_cm"]:.1f} cm.', fine=abs(error) < 1.5)
        step = self.axis_steps(p['lateral_mm'])
        if step:
            return step
        if m.release_due():
            return Advice('RELEASE SPACE NOW', 'The prescription is reached.')
        return Advice('HOLD SPACE, RELEASE AT THE YELLOW MARK', 'Peak on the tumour, on axis.', beam=True)

    def steps_isotope(self):
        m, p = self.m, self.m.p
        step = self.energy_step(p['w_ref'], 0.5 * p['w_tolerance']) or self.loss_steps()
        if step:
            return step
        if m.target_kw > 0.95 * p['target_kw']:
            return self.turn('source_ua', -1, 'LOWER THE SOURCE', 'The water is near its power limit.')
        if m.beam.extracted_ua < 1.1 * p['current_min'] and m.target_kw < 0.85 * p['target_kw'] and self.room_for_more():
            return self.turn('source_ua', 1, 'RAISE THE SOURCE', 'More current, more activity.')
        step = self.axis_steps(p['window_mm'])
        if step:
            return step
        return Advice('HOLD SPACE: IRRADIATE', 'On target, in limits. Activity builds up and saturates.', beam=True)

    def steps_experiment(self):
        m, p, b = self.m, self.m.p, self.m.beam
        step = self.energy_step(p['w_target'], 0.5 * p['w_tolerance'])
        if step:
            return step
        if b.w_spread > 0.9 * p['spread_max']:
            return self.turn('slit', -1, 'CLOSE THE ENERGY SLIT', 'The spread is above what the experiment accepts.', fine=m.knobs['slit'] < 0.004)
        if b.extracted_ua > 0.9 * p['current_max']:
            return self.turn('source_ua', -1, 'LOWER THE SOURCE', 'Above the current the detectors accept.')
        if b.extracted_ua < 1.5 * p['current_min']:
            return self.turn('source_ua', 1, 'RAISE THE SOURCE', 'Below the current the experiment needs.')
        step = self.axis_steps(p['position_mm'])
        if step:
            return step
        if b.spot > 0.9 * p['spot_max']:
            return self.quad_hunt(False, 'TURN THE QUADRUPOLE', 'The spot is too big: find the focus.')
        return Advice('HOLD SPACE: TAKE DATA', 'All lights green. Keep them green.', beam=True)

    def steps_isol(self):
        m, p, b = self.m, self.m.p, self.m.beam
        step = self.loss_steps()
        if step:
            return step
        if m.density() > 0.8 * p['density_limit']:
            return self.quad_hunt(True, 'TURN THE QUADRUPOLE', 'The power density is too high: open the spot.')
        if m.power() > 0.95 * p['power_limit_kw']:
            return self.turn('source_ua', -1, 'LOWER THE SOURCE', 'Near the power limit of the target.')
        if m.power() < 1.1 * p['power_min_kw'] and self.room_for_more():
            return self.turn('source_ua', 1, 'RAISE THE SOURCE', 'Power is what makes the isotopes.')
        step = self.axis_steps(0.5 * p['target_radius_mm'])
        if step:
            return step
        return Advice('HOLD SPACE: IRRADIATE', 'Power and spot within limits.', beam=True)

    def steps_qualification(self):
        m, p = self.m, self.m.p
        w = m.w_sample()
        if abs(w / p['w_sample'] - 1) > 0.5 * p['w_sample_tolerance']:
            more = w > p['w_sample']
            return self.turn('degrader_cm', 1 if more else -1, 'MORE DEGRADER' if more else 'LESS DEGRADER', f'{w:.0f} MeV on the sample, {p["w_sample"]:.0f} wanted.')
        step = self.axis_steps(3.0)
        if step:
            return step
        if m.edge_ratio() < p['uniformity'] + 0.02:
            return self.quad_hunt(True, 'TURN THE QUADRUPOLE', f'Edge over centre is {m.edge_ratio():.2f}: open the spot.')
        flux = m.flux()
        if flux > 0.8 * p['flux_max']:
            return self.turn('source_ua', -1, 'LOWER THE SOURCE', 'The flux is above the window.', fine=flux < 2 * p['flux_max'])
        if flux < 1.5 * p['flux_min']:
            return self.turn('source_ua', 1, 'RAISE THE SOURCE', 'The flux is below the window.', fine=flux > 0.5 * p['flux_min'])
        return Advice('HOLD SPACE: IRRADIATE', 'Energy, flux and uniformity in spec.', beam=True)
