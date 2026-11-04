"""The tutorial helper: one next step, from what the control room reads.
It never looks at the hidden errors.

These rules are also the scripted baseline of the Lab track. They are
one reasonable procedure written for this model, not the CIME operator
procedure (SPEC decision D4).
"""

import math
from dataclasses import dataclass

import numpy as np

from lions.game import physics

KNOB_HELP = {
    'rf_mhz': (
        'RF: the kick rhythm. A whole multiple of the turn rate, or the '
        'phase slips.'
    ),
    'b0': (
        'FIELD: sets the energy at extraction and the turn rate. The RF '
        'must follow.'
    ),
    'trim': (
        'TRIM: how fast the field grows with gamma. Right: same turn time '
        'at all radii.'
    ),
    'dee_kv': (
        'DEE VOLTAGE: gain per turn. More gain: last turns apart, less on '
        'the septum.'
    ),
    'phase_deg': 'PHASE: when the bunch arrives. On the crest it gains most.',
    'buncher': (
        'BUNCHER: packs the beam into the RF acceptance. Too strong '
        'overbunches.'
    ),
    'source_ua': (
        'SOURCE: current in. Current out, losses and power all scale with it.'
    ),
    'deflector_kv': (
        'DEFLECTOR: pulls the last turn out. It must match the beam rigidity.'
    ),
    'steer_x': 'STEER X: moves the beam sideways on the target.',
    'steer_y': 'STEER Y: moves the beam up and down on the target.',
    'quad': (
        'QUADRUPOLE: focuses. Smallest spot at the best value, wider on '
        'either side.'
    ),
    'slit': 'ENERGY SLIT: cuts the energy tails. Less spread, less current.',
    'degrader_cm': (
        'DEGRADER: more material, lower energy, more scattering, fewer '
        'protons.'
    ),
}


@dataclass
class Advice:
    """One step of advice: text, reason, and the knob or the keys it
    points to.
    """

    text: str
    why: str
    knob: str = ''  # key of the knob to turn; empty: nothing to turn
    direction: int = 0
    fine: bool = False
    beam: bool = False  # hold the beam on
    x: float = 0.0  # the keys to press, for the stages of level 0
    y: float = 0.0
    fire: bool = False


class Coach:
    """The helper of a mission: the next step from the control-room
    readings.
    """

    def __init__(self, mission):
        self.m = mission
        self.keys = {s['key']: s for s in mission.specs}
        self.quad_sign, self.last_spot, self.holding = 1, None, False

    def turn(self, key, direction, text, why, fine=False):
        """Advice to turn a knob in a direction."""
        return Advice(text, why, key, direction, fine)

    def phase_profile(self):
        """The phase probe at the radii with beam, or two zeros."""
        b = self.m.beam
        seen = np.flatnonzero(b.probe > 0)
        return b.phase[seen] if len(seen) > 1 else np.zeros(2)

    def rf_steps(self):
        """Get the beam out: frequency first, then the injection phase,
        from the phase probe.
        """
        phase = self.phase_profile()
        drift = phase[-1] - phase[0]
        third = max(len(phase) // 3, 1)
        extracted = self.m.has_beam()
        slope_in = phase[third] - phase[0]
        slope_out = phase[-1] - phase[-third - 1]
        if (
            'trim' in self.keys
            and extracted
            and abs(slope_out - slope_in) > 10
        ):
            return self.turn(
                'trim',
                1 if slope_out > slope_in else -1,
                'TURN TRIM',
                (
                    'The phase bends with radius: the field does not '
                    'follow gamma.'
                ),
                fine=abs(slope_out - slope_in) < 30,
            )
        if abs(drift) > (10 if extracted else 15):
            return self.turn(
                'rf_mhz',
                -1 if drift > 0 else 1,
                'LOWER THE RF' if drift > 0 else 'RAISE THE RF',
                (
                    f'The phase {"climbs" if drift > 0 else "falls"} with '
                    'radius: the RF is too '
                    f'{"high" if drift > 0 else "low"} for this field.'
                ),
                fine=abs(drift) < 40,
            )
        mean = float(phase.mean())
        if 'phase_deg' in self.keys and abs(mean) > 8:
            return self.turn(
                'phase_deg',
                -1 if mean > 0 else 1,
                'TURN THE INJECTION PHASE',
                (
                    'Flat phase, but off the crest: the beam gains less '
                    'and is more easily lost.'
                ),
                fine=abs(mean) < 20,
            )
        return None

    def energy_step(self, wanted, tolerance):
        """Advice on the field when the mean energy is off wanted by
        more than tolerance, relative.
        """
        w = self.m.beam.w_mean
        if w > 0 and abs(w / wanted - 1) > tolerance:
            up = w < wanted
            return self.turn(
                'b0',
                1 if up else -1,
                'RAISE THE FIELD' if up else 'LOWER THE FIELD',
                (
                    'The radius is fixed, so the field sets the energy. '
                    'The RF must follow it.'
                ),
                fine=abs(w / wanted - 1) < 3 * tolerance,
            )
        return None

    def axis_steps(self, limit):
        """Advice to steer, or to fix the deflector, when the beam is
        off axis, else None.
        """
        r = self.m.readings
        if not self.m.has_beam():
            return None
        for key, reading, label in (
            ('steer_x', r['x_mm'], 'STEER X'),
            ('steer_y', r['y_mm'], 'STEER Y'),
        ):
            if key in self.keys and abs(reading) > 0.5 * limit:
                if (
                    key == 'steer_x'
                    and abs(reading) > 12
                    and 'deflector_kv' in self.keys
                ):
                    return self.turn(
                        'deflector_kv',
                        -1 if reading > 0 else 1,
                        'TURN THE DEFLECTOR',
                        (
                            'The beam leaves far off axis: the deflector '
                            'does not match its rigidity.'
                        ),
                    )
                return self.turn(
                    key,
                    -1 if reading > 0 else 1,
                    label,
                    f'The beam is {abs(reading):.1f} mm off axis.',
                    fine=abs(reading) < 2 * limit,
                )
        return None

    def quad_hunt(self, want_bigger, text, why):
        """The spot has no sign: turn, watch, turn back if it went the
        wrong way.
        """
        spot = self.m.spot()
        if (
            self.last_spot is not None
            and abs(spot - self.last_spot) > 1e-6
            and (spot > self.last_spot) != want_bigger
        ):
            self.quad_sign = -self.quad_sign
        self.last_spot = spot
        return self.turn('quad', self.quad_sign, text, why)

    def advise(self):
        """The next advice: hold when in spec, else the RF steps, then
        those of the mission.
        """
        ready = self.m.ready()
        if (
            self.holding
            and ready
            and not getattr(self.m, 'release_due', lambda: False)()
        ):  # in spec: stop retuning
            return Advice(
                'KEEP HOLDING SPACE', 'In spec. Watch the readings.', beam=True
            )
        self.holding = False
        step = None if ready else self.rf_steps()
        if not step:
            step = (
                self.mission_steps() if self.m.has_beam() else self.lost_beam()
            )
        self.holding = step.beam
        return step

    def lost_beam(self):
        """Advice when no beam is extracted: buncher near 1, field, or
        watch the probe.
        """
        reach = self.m.reach()
        if 'buncher' in self.keys and abs(self.m.knobs['buncher'] - 1) > 0.2:
            return self.turn(
                'buncher',
                1 if self.m.knobs['buncher'] < 1 else -1,
                'SET THE BUNCHER NEAR 1',
                'Bunched beam fits the RF acceptance.',
            )
        if reach < 0.9 and 'b0' in self.keys:
            return self.turn(
                'b0',
                1,
                'RAISE THE FIELD',
                'Flat phase but the beam stops inside: check the field.',
            )
        return Advice(
            'WATCH THE RADIAL PROBE', 'The beam stops inside the machine.'
        )

    def loss_steps(self):
        """Advice to raise the dee voltage or lower the source when
        losses are high, else None.
        """
        m = self.m
        if m.beam.lost_kw > 0.8 * m.machine.max_lost_kw:
            if (
                'dee_kv' in self.keys
                and m.knobs['dee_kv'] < self.keys['dee_kv']['max']
            ):
                return self.turn(
                    'dee_kv',
                    1,
                    'RAISE THE DEE VOLTAGE',
                    (
                        'More gain per turn separates the last turns: '
                        'less beam on the septum.'
                    ),
                )
            if 'source_ua' in self.keys:
                return self.turn(
                    'source_ua',
                    -1,
                    'LOWER THE SOURCE',
                    'Too much beam lost in the machine.',
                )
        return None

    def room_for_more(self):
        """Raising the source must not undo the loss step: a dead band
        between the two.
        """
        return self.m.beam.lost_kw < 0.65 * self.m.machine.max_lost_kw

    def mission_steps(self):
        """Advice on the buncher, else the steps of the mission."""
        if 'buncher' in self.keys and abs(self.m.knobs['buncher'] - 1) > 0.2:
            return self.turn(
                'buncher',
                1 if self.m.knobs['buncher'] < 1 else -1,
                'SET THE BUNCHER NEAR 1',
                (
                    'A bunched beam fits the RF acceptance: more current '
                    'for the same source.'
                ),
            )
        return getattr(self, f'steps_{self.m.key}')()

    def steps_therapy(self):
        """Advice for the therapy mission: peak depth, axis, release."""
        m, p = self.m, self.m.p
        error = m.peak - p['depth_cm']
        if abs(error) > 0.5 * p['half_width_cm']:
            return self.turn(
                'degrader_cm',
                1 if error > 0 else -1,
                'MORE DEGRADER' if error > 0 else 'LESS DEGRADER',
                (
                    f'The Bragg peak is at {m.peak:.1f} cm, the tumour at '
                    f'{p["depth_cm"]:.1f} cm.'
                ),
                fine=abs(error) < 1.5,
            )
        step = self.axis_steps(p['lateral_mm'])
        if step:
            return step
        if m.release_due():
            return Advice('RELEASE SPACE NOW', 'The prescription is reached.')
        return Advice(
            'HOLD SPACE, RELEASE AT THE YELLOW MARK',
            'Peak on the tumour, on axis.',
            beam=True,
        )

    def steps_isotope(self):
        """Advice for the isotope mission: energy, losses, source
        current, axis.
        """
        m, p = self.m, self.m.p
        step = (
            self.energy_step(p['w_ref'], 0.5 * p['w_tolerance'])
            or self.loss_steps()
        )
        if step:
            return step
        if m.target_kw > 0.95 * p['target_kw']:
            return self.turn(
                'source_ua',
                -1,
                'LOWER THE SOURCE',
                'The water is near its power limit.',
            )
        if (
            m.beam.extracted_ua < 1.1 * p['current_min']
            and m.target_kw < 0.85 * p['target_kw']
            and self.room_for_more()
        ):
            return self.turn(
                'source_ua',
                1,
                'RAISE THE SOURCE',
                'More current, more activity.',
            )
        step = self.axis_steps(p['window_mm'])
        if step:
            return step
        return Advice(
            'HOLD SPACE: IRRADIATE',
            'On target, in limits. Activity builds up and saturates.',
            beam=True,
        )

    def steps_experiment(self):
        """Advice for the experiment mission: energy, spread, current,
        position, spot.
        """
        m, p, b = self.m, self.m.p, self.m.beam
        step = self.energy_step(p['w_target'], 0.5 * p['w_tolerance'])
        if step:
            return step
        if b.w_spread > 0.9 * p['spread_max']:
            return self.turn(
                'slit',
                -1,
                'CLOSE THE ENERGY SLIT',
                'The spread is above what the experiment accepts.',
                fine=m.knobs['slit'] < 0.004,
            )
        if b.extracted_ua > 0.9 * p['current_max']:
            return self.turn(
                'source_ua',
                -1,
                'LOWER THE SOURCE',
                'Above the current the detectors accept.',
            )
        if b.extracted_ua < 1.5 * p['current_min']:
            return self.turn(
                'source_ua',
                1,
                'RAISE THE SOURCE',
                'Below the current the experiment needs.',
            )
        step = self.axis_steps(p['position_mm'])
        if step:
            return step
        if b.spot > 0.9 * p['spot_max']:
            return self.quad_hunt(
                False,
                'TURN THE QUADRUPOLE',
                'The spot is too big: find the focus.',
            )
        return Advice(
            'HOLD SPACE: TAKE DATA',
            'All lights green. Keep them green.',
            beam=True,
        )

    def steps_isol(self):
        """Advice for the isol mission: losses, power density, power,
        position.
        """
        m, p = self.m, self.m.p
        step = self.loss_steps()
        if step:
            return step
        if m.density() > 0.8 * p['density_limit']:
            return self.quad_hunt(
                True,
                'TURN THE QUADRUPOLE',
                'The power density is too high: open the spot.',
            )
        if m.power() > 0.95 * p['power_limit_kw']:
            return self.turn(
                'source_ua',
                -1,
                'LOWER THE SOURCE',
                'Near the power limit of the target.',
            )
        if m.power() < 1.1 * p['power_min_kw'] and self.room_for_more():
            return self.turn(
                'source_ua',
                1,
                'RAISE THE SOURCE',
                'Power is what makes the isotopes.',
            )
        step = self.axis_steps(0.5 * p['target_radius_mm'])
        if step:
            return step
        return Advice(
            'HOLD SPACE: IRRADIATE', 'Power and spot within limits.', beam=True
        )

    def steps_qualification(self):
        """Advice for the qualification mission: energy, position,
        uniformity, flux.
        """
        m, p = self.m, self.m.p
        w = m.w_sample()
        if abs(w / p['w_sample'] - 1) > 0.5 * p['w_sample_tolerance']:
            more = w > p['w_sample']
            return self.turn(
                'degrader_cm',
                1 if more else -1,
                'MORE DEGRADER' if more else 'LESS DEGRADER',
                f'{w:.0f} MeV on the sample, {p["w_sample"]:.0f} wanted.',
            )
        step = self.axis_steps(3.0)
        if step:
            return step
        if m.edge_ratio() < p['uniformity'] + 0.02:
            return self.quad_hunt(
                True,
                'TURN THE QUADRUPOLE',
                f'Edge over centre is {m.edge_ratio():.2f}: open the spot.',
            )
        flux = m.flux()
        if flux > 0.8 * p['flux_max']:
            return self.turn(
                'source_ua',
                -1,
                'LOWER THE SOURCE',
                'The flux is above the window.',
                fine=flux < 2 * p['flux_max'],
            )
        if flux < 1.5 * p['flux_min']:
            return self.turn(
                'source_ua',
                1,
                'RAISE THE SOURCE',
                'The flux is below the window.',
                fine=flux > 0.5 * p['flux_min'],
            )
        return Advice(
            'HOLD SPACE: IRRADIATE',
            'Energy, flux and uniformity in spec.',
            beam=True,
        )


class StageCoach:
    """The helper of the level 0 stages: what the screen shows, and
    the key to press.
    """

    def __init__(self, stage):
        self.s = stage

    def advise(self):
        """The advice of the method named by the stage key."""
        return getattr(self, self.s.key)()

    def larmor(self):
        """Advice for the Larmor stage: field, energy, then fire."""
        s = self.s
        if s.flight is not None:
            return Advice(
                'WATCH THE ORBIT',
                (
                    'rho = p / (q B): this shot flies at '
                    f'{s.shot_radius:.2f} m, the ring is {s.target:.2f} m.'
                ),
            )
        wanted = -s.target_rotation * s.p['b_max']  # rotation +1 needs Bz < 0
        needed = (
            physics.momentum(s.w)
            * physics.PROTON_KG
            * physics.C
            / (physics.E_CHARGE * s.target)
        )
        if needed > s.p['b_max']:
            return Advice(
                'HOLD LEFT: LESS ENERGY',
                (
                    'This ring needs more field than the magnet has: '
                    'lower the momentum.'
                ),
                x=-1,
            )
        b = math.copysign(needed, wanted)
        if abs(b - s.b) > 1.5 * s.p['b_rate'] / 60:
            up = b > s.b
            return Advice(
                f"HOLD {'UP' if up else 'DOWN'}: FIELD TO {b:+.2f} T",
                (
                    f'At {s.w:.1f} MeV this ring needs |B| = {needed:.2f} '
                    'T. The arrow on the ring sets the sign.'
                ),
                y=1 if up else -1,
            )
        return Advice(
            'PRESS SPACE: FIRE',
            'Field and energy match the ring.',
            fire=not s.last_action.fire,
        )

    def dee_beat(self):
        """Advice for the dee beat stage: move the RF toward the
        revolution frequency.
        """
        s = self.s
        f_rev = physics.revolution_frequency(s.w, s.p['b'])
        wanted = f_rev * (1 - s.phi / math.pi / 4)
        step = s.p['f_rate'] * s.f_nominal / 60
        if abs(wanted - s.f) < step:
            return Advice(
                'HOLD NOTHING: ON THE BEAT',
                'The needle stands still on the crest: every crossing gains.',
            )
        up = wanted > s.f
        way = (
            'the needle must come back to the green arc'
            if abs(s.phi) > 0.3
            else 'the needle must stop moving'
        )
        return Advice(
            f"HOLD {'RIGHT: MORE' if up else 'LEFT: LESS'} RF",
            (
                f'{way.capitalize()}. A needle that moves means f is not '
                'q B / (2 pi m).'
            ),
            x=1 if up else -1,
        )

    def trim_coils(self):
        """Advice for the trim coils stage: move the field trim to
        follow gamma.
        """
        s = self.s
        wanted = s.needed_delta() + s.phi / 40
        step = s.p['delta_rate'] / 60
        if abs(wanted - s.delta) < step:
            return Advice(
                'HOLD NOTHING: IN STEP',
                'The trim follows gamma: B = gamma B0.',
            )
        up = wanted > s.delta
        return Advice(
            f"HOLD {'UP: MORE' if up else 'DOWN: LESS'} FIELD",
            (
                f'Gamma is {physics.gamma(s.w):.3f}: the field should be '
                f'{100 * s.needed_delta():+.2f}%. The grey needle is what '
                'no trim would do.'
            ),
            y=1 if up else -1,
        )

    def fallout(self):
        """Advice for the fallout stage: wait, or kick the ion back."""
        s = self.s
        rate = math.sqrt(physics.isochronous_index(s.w))
        growing = s.z + s.dz / rate  # the growing mode of z'' = k z
        if abs(growing) <= 1.0 or s.last_action.y:
            return Advice(
                'WAIT',
                (
                    f'Qz^2 = {-rate * rate:+.3f}: the ion drifts away '
                    'from the median plane, faster and faster.'
                ),
            )
        up = growing < 0
        return Advice(
            f"TAP {'UP' if up else 'DOWN'}: KICK BACK",
            (
                'Kick early: a small kick now saves a big one later. '
                'Charges run out.'
            ),
            y=1 if up else -1,
        )

    def thomas(self):
        """Advice for the Thomas stage: flutter and spiral to keep Qz in
        its band.
        """
        s = self.s
        low, high = s.p['band_low'], s.p['band_high']
        if s.qz2 < low**2 and s.flutter < s.p['f_max']:
            return Advice(
                'HOLD UP: MORE FLUTTER',
                (
                    'Energy raises k, so Qz falls: Qz^2 = -k + '
                    'N^2/(N^2-1) F (1 + 2 tan^2 xi).'
                ),
                y=1,
            )
        if s.qz2 < low**2:
            return Advice(
                'HOLD RIGHT: MORE SPIRAL',
                'Flutter is at its limit: the spiral multiplies it.',
                x=1,
            )
        if s.qz2 > high**2:
            return Advice(
                'HOLD DOWN: LESS FLUTTER',
                'Stay below Qz = 1, where field errors drive the beam.',
                y=-1,
            )
        return Advice(
            'STEADY: IN THE GREEN',
            'Qz between the two red lines. Watch it as the energy climbs.',
        )


def coach_for(stage):
    """The Coach for a mission, the StageCoach for a level 0 stage."""
    return Coach(stage) if hasattr(stage, 'specs') else StageCoach(stage)
