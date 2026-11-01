"""Missions: a client request, a machine with knobs, control-room readings. Same reset, step, score as the stages."""
import math

import numpy as np

from lions.constants import E_CHARGE
from lions.game import machine as mc
from lions.game import physics
from lions.game.stages import STAGES, Stage

MACHINE_UNITS = {'radius': 'm', 'sectors': '1', 'flutter': '1', 'spiral_deg': 'deg', 'dees': '1', 'dee_angle_deg': 'deg',
                 'w_injection': 'MeV/u', 'z_injection': 'mm', 'z_aperture': 'mm', 'deflector_scale': 'kV/MV', 'r_size': 'mm',
                 'septum': 'mm', 'deflector_gain': 'mm', 'spot_base': 'mm', 'spot_quad_width': '1', 'drift': 'm',
                 'theta_accept': 'rad', 'space_charge': '1/uA', 'max_lost_kw': 'kW'}
KNOB_UNITS = {'rf_mhz': 'MHz', 'dee_kv': 'kV', 'phase_deg': 'deg', 'buncher': '1', 'source_ua': 'uA', 'b0': 'T',
              'trim': '1', 'deflector_kv': 'kV', 'steer_x': 'mm', 'steer_y': 'mm', 'quad': '1', 'slit': '1', 'degrader_cm': 'cm'}
COMMON_UNITS = {'w_target': 'MeV', 'dee_kv': 'kV', 'source_ua': 'uA', 'time_limit': 's', 'noise': '1', 'drift': '1/s^0.5',
                'field_error': '1', 'steer_error': 'mm', 'quad_error': '1', 'phase_error': 'deg', 'heat_rate': '1/s'}
REFRESH = 0.15    # s of game time between passes of the machine model


def checked(block, units, where):
    out = {}
    for key, unit in units.items():
        item = block[key]
        if item['unit'] != unit:
            raise ValueError(f"{where}.{key}: unit {item['unit']}, expected {unit}")
        out[key] = item['value']
    return out


class Mission(Stage):
    GOAL_UNITS = {}

    def __init__(self, scenario, entry, params):
        super().__init__(scenario, entry, params)
        self.machine = mc.Machine(**checked(entry['machine'], MACHINE_UNITS, f'{self.key}.machine'))
        self.ion = mc.Ion(**entry['ion'])
        self.client, self.mission_id = entry['client'], entry['id']
        self.specs = entry['knobs']
        for spec in self.specs:
            if spec['unit'] != KNOB_UNITS[spec['key']]:
                raise ValueError(f"{self.key}.knob {spec['key']}: unit {spec['unit']}")
        self.harmonic = entry['harmonic']

    def ideal(self):
        """Knobs a perfect operator would reach, hidden errors compensated. Knobs left out of a mission stay here."""
        p, m, e = self.p, self.machine, self.errors
        b0, f = mc.nominal(m, self.ion, p['w_target'])
        return {'rf_mhz': self.harmonic * f / 1e6, 'dee_kv': p['dee_kv'], 'phase_deg': -e.phase_deg, 'buncher': 1.0,
                'source_ua': p['source_ua'], 'b0': b0 / (1 + e.field), 'trim': 1.0,
                'deflector_kv': mc.deflector_nominal(m, self.ion, p['w_target']), 'steer_x': -e.steer_x, 'steer_y': -e.steer_y,
                'quad': -e.quad, 'slit': 0.0, 'degrader_cm': 0.0}

    def setup(self):
        p, r = self.p, self.rng
        generator = np.random.default_rng(self.seed)
        self.draws = generator.standard_normal((mc.MACROPARTICLES, 3))
        self.noise_rng = generator
        self.errors = mc.Errors(field=r.gauss(0, p['field_error']), steer_x=r.gauss(0, p['steer_error']), steer_y=r.gauss(0, p['steer_error']),
                                quad=r.gauss(0, p['quad_error']), phase_deg=r.gauss(0, p['phase_error']))
        self.knobs = self.ideal()
        for spec in self.specs:
            self.knobs[spec['key']] = self.start_value(spec)
        self.selected, self.heat, self.delivered, self.since_refresh = 0, 0.0, 0.0, 0.0
        self.beam_on, self.readings = False, {}
        self.setup_goal()
        self.refresh()

    def start_value(self, spec):
        """Where the shift starts: off the best setting, within the knob range."""
        offset = self.rng.uniform(0.5, 1) * self.rng.choice((-1, 1)) * spec['start_error']
        return self.clip(spec, self.knobs[spec['key']] + offset)

    @staticmethod
    def clip(spec, value):
        return min(max(value, spec['min']), spec['max'])

    def setup_goal(self):
        pass

    def refresh(self):
        self.beam = mc.run(self.machine, self.ion, self.knobs, self.errors, self.draws)
        self.after_beam()
        noise = lambda v: v * (1 + self.noise_rng.normal(0, self.p['noise']))  # noqa: E731
        b = self.beam
        self.readings = {'extracted_uA': noise(b.extracted_ua), 'transmission': noise(b.transmission), 'w_MeV': noise(b.w_mean),
                         'spread': b.w_spread, 'x_mm': b.x + self.noise_rng.normal(0, self.p['noise']),
                         'y_mm': b.y + self.noise_rng.normal(0, self.p['noise']), 'spot_mm': noise(self.spot()), 'lost_kW': noise(b.lost_kw)}

    def after_beam(self):
        pass

    def spot(self):
        return self.beam.spot

    def advance(self, action, dt):
        p = self.p
        if self.y_edge:
            self.selected = (self.selected - int(math.copysign(1, self.y_edge))) % len(self.specs)
        spec = self.specs[self.selected]
        if action.x:
            self.knobs[spec['key']] = self.clip(spec, self.knobs[spec['key']] + action.x * spec['rate'] * dt)
        self.errors.field += self.noise_rng.normal(0, p['drift'] * math.sqrt(dt))
        self.since_refresh += dt
        if self.since_refresh >= REFRESH:
            self.since_refresh = 0.0
            self.refresh()
        self.beam_on = action.fire
        if self.beam_on:
            self.heat = max(self.heat + dt * p['heat_rate'] * (self.beam.lost_kw / self.machine.max_lost_kw - 1), 0.0)
            self.deliver(dt)
            if self.heat >= 1:
                self.finish(False, 'INTERLOCK: TOO MUCH BEAM LOST IN THE MACHINE')
        else:
            self.heat = max(self.heat - dt * p['heat_rate'], 0.0)
        if not self.done and self.time >= p['time_limit']:
            self.finish(False, 'SHIFT OVER')

    def observation(self):
        return {**{f'knob_{k}': float(self.knobs[k]) for k in (s['key'] for s in self.specs)}, **self.readings, 'heat': self.heat,
                'beam_on': float(self.beam_on), 'delivered': self.delivered}

    def metrics(self):
        return {'delivered': self.delivered, 'heat': self.heat, 'time_s': self.time, 'quality': self.quality()}

    def deliver(self, dt):
        raise NotImplementedError

    def quality(self):
        return 0.0


class Therapy(Mission):
    """Water phantom: the Bragg peak on the tumour, the dose as prescribed, nothing behind. Not a clinical model."""
    GOAL_UNITS = {'depth_cm': 'cm', 'half_width_cm': 'cm', 'lateral_mm': 'mm', 'dose_gy': 'Gy', 'dose_tolerance': '1', 'healthy_limit_gy': 'Gy',
                  'field_area_cm2': 'cm2', 'real_s_per_s': '1'}
    UNITS = {**COMMON_UNITS, **GOAL_UNITS}

    def ideal(self):
        out = super().ideal()
        out['degrader_cm'] = physics.proton_range(self.p['w_target']) - self.p['depth_cm']
        return out

    def setup_goal(self):
        self.healthy, self.depth_axis = 0.0, np.linspace(0, 1.15 * physics.proton_range(self.p['w_target']), 240)
        self.bragg = np.zeros_like(self.depth_axis)
        self.peak, self.rate, self.theta, self.deg_transmission = 0.0, 0.0, 0.0, 0.0

    def after_beam(self):
        p, b = self.p, self.beam
        t = self.knobs['degrader_cm']
        self.bragg[:] = 0
        if not len(b.energies):
            self.peak, self.rate, self.theta, self.deg_transmission = 0.0, 0.0, 0.0, 0.0
            return
        w = float(b.w_mean)
        self.theta = physics.scattering_angle(w, t)
        self.deg_transmission = 1 / (1 + (self.theta / self.machine.theta_accept)**2)
        ranges = np.maximum(physics.proton_range(b.energies) - t, 1e-3)
        self.peak = float(ranges.mean())
        sigma = float(np.hypot(ranges.std(), physics.straggling(self.peak) + (physics.straggling(max(t, 1e-3)) if t > 0 else 0)))
        mean = self.peak
        # Depth dose of a proton: dE/dx = E^(1-p) / (p alpha), folded with the range spread.
        e = np.maximum(mean - self.depth_axis, 1e-3)
        stopping = e**(1 / physics.BRAGG_P - 1) / (physics.BRAGG_P * physics.BRAGG_ALPHA**(1 / physics.BRAGG_P))
        stopping[self.depth_axis > mean] = 0
        step = self.depth_axis[1] - self.depth_axis[0]
        width = max(int(4 * sigma / step), 1)
        shape = np.exp(-0.5 * (np.arange(-width, width + 1) * step / sigma)**2)
        self.bragg[:] = np.convolve(stopping, shape / shape.sum(), mode='same')
        protons = b.extracted_ua * 1e-6 * self.deg_transmission / E_CHARGE
        at_peak = float(self.bragg.max())
        self.rate = physics.GY_PER_FLUENCE * protons / p['field_area_cm2'] * at_peak    # Gy/s at the peak

    def spot(self):
        drift_mm = self.machine.drift * 1e3 * self.theta
        return float(np.hypot(self.beam.spot, drift_mm))

    def on_tumour(self):
        p = self.p
        return abs(self.peak - p['depth_cm']) <= p['half_width_cm'] and math.hypot(self.beam.x, self.beam.y) <= p['lateral_mm']

    def deliver(self, dt):
        p = self.p
        dose = self.rate * dt * p['real_s_per_s']
        if self.on_tumour():
            self.delivered += dose
        else:
            self.healthy += dose
        if self.healthy >= p['healthy_limit_gy']:
            self.finish(False, 'TOO MUCH DOSE OUTSIDE THE TUMOUR')
        elif self.delivered >= p['dose_gy'] * (1 + p['dose_tolerance']):
            self.finish(False, 'OVERDOSE: RELEASE THE BEAM IN TIME')
        self.points = int(3000 * min(self.delivered / p['dose_gy'], 1))

    def advance(self, action, dt):
        was_on = self.beam_on
        super().advance(action, dt)
        p = self.p
        if was_on and not self.beam_on and not self.done and abs(self.delivered / p['dose_gy'] - 1) <= p['dose_tolerance']:
            self.points = int(3000 + 3000 * self.quality() + 20 * max(p['time_limit'] - self.time, 0))
            self.finish(True, 'PRESCRIPTION DELIVERED')

    def quality(self):
        total = self.delivered + self.healthy
        return self.delivered / total if total else 0.0


class Isotope(Mission):
    """18O(p,n)18F in enriched water: A' = lambda (Y I - A), Y = 8 GBq/uA at 18 MeV (Silva, Niroomand-Rad; see the scenario)."""
    GOAL_UNITS = {'activity_gbq': 'GBq', 'yield_gbq_ua': 'GBq/uA', 'w_ref': 'MeV', 'w_threshold': 'MeV', 'half_life_min': 'min',
                  'real_min_per_s': 'min/s', 'target_kw': 'kW', 'window_mm': 'mm'}
    UNITS = {**COMMON_UNITS, **GOAL_UNITS}

    def setup_goal(self):
        self.activity, self.target_kw = 0.0, 0.0

    def after_beam(self):
        self.target_kw = self.beam.extracted_ua * self.beam.w_mean * 1e-3

    def on_window(self):
        return math.hypot(self.beam.x, self.beam.y) <= self.p['window_mm']

    def deliver(self, dt):
        p, b = self.p, self.beam
        minutes = dt * p['real_min_per_s']
        lam = math.log(2) / p['half_life_min']
        rate = 0.0
        if self.on_window():
            fraction = max((b.w_mean - p['w_threshold']) / (p['w_ref'] - p['w_threshold']), 0.0)
            rate = p['yield_gbq_ua'] * fraction * b.extracted_ua
        self.activity += minutes * lam * (rate - self.activity)
        self.delivered = self.activity
        if self.target_kw > p['target_kw']:
            self.finish(False, 'TARGET BOILED: TOO MUCH POWER ON THE WATER')
        elif not self.on_window() and b.extracted_ua > 1:
            self.finish(False, 'BEAM ON THE TARGET WINDOW: WINDOW BROKEN')
        elif self.activity >= p['activity_gbq']:
            self.points = int(3000 + 30 * max(p['time_limit'] - self.time, 0))
            self.finish(True, 'BATCH READY FOR THE HOSPITAL')
        else:
            self.points = int(3000 * self.activity / p['activity_gbq'])

    def quality(self):
        return self.activity / self.p['activity_gbq']


class Experiment(Mission):
    """A nuclear physics run: energy, spread and current in spec, held for the requested time."""
    GOAL_UNITS = {'w_tolerance': '1', 'spread_max': '1', 'current_min': 'uA', 'current_max': 'uA', 'position_mm': 'mm', 'spot_max': 'mm', 'hold_s': 's'}
    UNITS = {**COMMON_UNITS, **GOAL_UNITS}

    def in_spec(self):
        p, b = self.p, self.beam
        return (b.w_mean > 0 and abs(b.w_mean / p['w_target'] - 1) <= p['w_tolerance'] and b.w_spread <= p['spread_max']
                and p['current_min'] <= b.extracted_ua <= p['current_max'] and math.hypot(b.x, b.y) <= p['position_mm'] and b.spot <= p['spot_max'])

    def deliver(self, dt):
        p = self.p
        if self.in_spec():
            self.delivered += dt
            self.points = int(5000 * self.delivered / p['hold_s'])
        if self.delivered >= p['hold_s']:
            self.points = int(5000 + 30 * max(p['time_limit'] - self.time, 0))
            self.finish(True, 'THE EXPERIMENT HAS ITS DATA')

    def quality(self):
        return self.delivered / self.p['hold_s']


class IsolTarget(Mission):
    """A production target: as much beam power as it survives, spread over the spot."""
    GOAL_UNITS = {'power_limit_kw': 'kW', 'density_limit': 'kW/mm2', 'target_radius_mm': 'mm', 'energy_goal_kws': 'kW s'}
    UNITS = {**COMMON_UNITS, **GOAL_UNITS}

    def power(self):
        return self.beam.extracted_ua * self.beam.w_mean * 1e-3

    def density(self):
        return self.power() / (2 * math.pi * max(self.beam.spot, 1e-3)**2)

    def deliver(self, dt):
        p, b = self.p, self.beam
        centred = math.hypot(b.x, b.y) + b.spot <= p['target_radius_mm']
        if self.power() > p['power_limit_kw']:
            self.finish(False, 'TARGET OVER ITS POWER LIMIT')
        elif self.density() > p['density_limit']:
            self.finish(False, 'HOT SPOT: THE TARGET MELTED')
        elif b.extracted_ua > 1 and not centred:
            self.finish(False, 'BEAM OFF THE TARGET')
        else:
            self.delivered += self.power() * dt
            self.points = int(5000 * self.delivered / p['energy_goal_kws'])
            if self.delivered >= p['energy_goal_kws']:
                self.points += int(30 * max(p['time_limit'] - self.time, 0))
                self.finish(True, 'TARGET IRRADIATED')

    def quality(self):
        return self.delivered / self.p['energy_goal_kws']


STAGES.update({'therapy': Therapy, 'isotope': Isotope, 'experiment': Experiment, 'isol': IsolTarget})
