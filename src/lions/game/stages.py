"""Level 0 stages: plain models, no Qt. A human at the arcade, a script or an RL agent plays the same reset, step, score."""
import math
import random
from collections import namedtuple
from importlib import resources

import yaml

from lions.game import physics

Action = namedtuple('Action', 'x y fire', defaults=(0, 0, False))   # x, y in [-1, 1]


def load_scenario(name='level0'):
    """Stage objects from a scenario file. A parameter whose unit differs from the one the stage expects is an error."""
    import lions.game.missions  # noqa: F401, E402   registers the mission stages
    text = resources.files('lions.game').joinpath('scenarios', f'{name}.yaml').read_text()
    scenario = yaml.safe_load(text)
    stages = []
    for entry in scenario['stages']:
        cls = STAGES[entry['stage']]
        params = {}
        for key, unit in cls.UNITS.items():
            item = entry['params'][key]
            if item['unit'] != unit:
                raise ValueError(f"{entry['stage']}.{key}: unit {item['unit']}, expected {unit}")
            params[key] = item['value']
        stages.append(cls(scenario['scenario'], entry, params))
    return stages


class Stage:
    UNITS = {}

    def __init__(self, scenario, entry, params):
        self.scenario, self.key, self.p = scenario, entry['stage'], params
        self.title, self.brief, self.lesson = entry['title'], entry['brief'], entry['lesson']
        self.mandatory = entry.get('mandatory', True)    # False: failing costs no life
        self.tutorial = False     # no clock and no shot limit: the arcade pairs it with a coach
        self.switches = entry.get('switches', {})

    def reset(self, seed):
        self.seed, self.rng = seed, random.Random(seed)
        self.time, self.points, self.steps, self.knob_changes = 0.0, 0, 0, 0
        self.done, self.success, self.message = False, False, ''
        self.last_action = Action()
        self.setup()
        return self.observation()

    def step(self, action, dt):
        """Advance by dt seconds of game time. Returns observation, metrics, done."""
        if self.done:
            return self.observation(), self.metrics(), True
        last = self.last_action
        self.knob_changes += (action.x != 0 and last.x == 0) + (action.y != 0 and last.y == 0)
        self.fire_edge = action.fire and not last.fire
        self.y_edge = action.y if action.y != 0 and last.y == 0 else 0
        self.last_action = action
        self.time += dt
        self.steps += 1
        self.advance(action, dt)
        return self.observation(), self.metrics(), self.done

    def finish(self, success, message):
        self.done, self.success, self.message = True, success, message

    def score(self):
        """The record of a run: scenario, seed, switches, budget, metrics (SPEC section 7). No player identity."""
        return {'scenario': self.scenario, 'stage': self.key, 'seed': self.seed, 'switches': self.switches,
                'budget': {'knob_changes': self.knob_changes, 'steps': self.steps, 'time_s': self.time},
                'success': self.success, 'points': self.points, 'metrics': self.metrics()}

    def setup(self):
        raise NotImplementedError

    def advance(self, action, dt):
        raise NotImplementedError

    def observation(self):
        raise NotImplementedError

    def metrics(self):
        return {}


class Larmor(Stage):
    UNITS = {'b_max': 'T', 'b_min': 'T', 'b_start': 'T', 'b_rate': 'T/s', 'w_min': 'MeV', 'w_max': 'MeV', 'w_start': 'MeV',
             'w_rate': 'MeV/s', 'target_min': 'm', 'target_max': 'm', 'tolerance': '1', 'targets': '1', 'shots': '1',
             'flight_time': 's'}

    def setup(self):
        p = self.p
        self.b, self.w = p['b_start'], p['w_start']
        self.hits, self.shots = 0, 0
        self.flight, self.shot_radius, self.shot_rotation = None, None, None
        self.new_target()

    def new_target(self):
        self.target = self.rng.uniform(self.p['target_min'], self.p['target_max'])
        self.target_rotation = self.rng.choice((-1, 1))

    def advance(self, action, dt):
        p = self.p
        if self.flight is not None:
            self.flight += dt / p['flight_time']
            if self.flight >= 1:
                self.flight = None
                self.judge()
            return
        self.b = min(max(self.b + action.y * p['b_rate'] * dt, -p['b_max']), p['b_max'])
        self.w = min(max(self.w + action.x * p['w_rate'] * dt, p['w_min']), p['w_max'])
        if self.fire_edge:
            self.shots += 1
            if abs(self.b) < p['b_min']:
                self.message = 'NO FIELD: THE ION FLIES STRAIGHT'
                self.shot_radius = None
                self.after_shot()
                return
            self.shot_radius, self.shot_rotation = physics.radius(self.w, self.b), physics.rotation(self.b)
            self.flight = 0.0

    def judge(self):
        error = self.shot_radius / self.target - 1
        if self.shot_rotation != self.target_rotation:
            self.message = 'WRONG WAY: FLIP THE SIGN OF B'
        elif abs(error) > self.p['tolerance']:
            self.message = f'RADIUS {100 * error:+.0f}%: TOO BIG' if error > 0 else f'RADIUS {100 * error:+.0f}%: TOO SMALL'
        else:
            self.hits += 1
            self.points += int(500 * (2 - abs(error) / self.p['tolerance']))
            self.message = 'LOCKED'
            if self.hits < self.p['targets']:
                self.new_target()
        self.after_shot()

    def after_shot(self):
        if self.hits >= self.p['targets']:
            self.finish(True, 'ALL RINGS LOCKED')
        elif self.shots >= self.p['shots'] and not self.tutorial:
            self.finish(False, 'OUT OF SHOTS')

    def observation(self):
        return {'b_T': self.b, 'w_MeV': self.w, 'target_m': self.target, 'target_rotation': self.target_rotation,
                'hits': self.hits, 'shots': self.shots}

    def metrics(self):
        return {'hits': self.hits, 'shots': self.shots}


class GapStage(Stage):
    """An ion crossing a gap twice per turn: gain q V cos(phi), phase slip pi (f_rf / f_rev - 1) per crossing."""

    def setup_gap(self, phi):
        self.w, self.phi, self.crossings, self.combo, self.crest = self.p['w_start'], phi, 0, 0, 0
        self.frac, self.arcs = 0.0, []
        self.arcs = [self.arc(0.0)]

    def arc(self, x_start):
        """The half circle after a crossing, on the gap line: (x_start, x_end). Odd crossings go to +x, under the line."""
        direction = 1 if len(self.arcs) % 2 else -1
        return x_start, x_start + direction * 2 * physics.radius(self.w, self.field())

    def cross(self):
        p = self.p
        self.w = max(self.w + p['v_gap'] * math.cos(self.phi), 0.0)
        self.crossings += 1
        if math.cos(self.phi) >= p['crest_cos']:
            self.combo += 1
            self.crest += 1
            self.points += 5 * min(self.combo, 20)
        else:
            self.combo = 0
        if self.w <= p['w_lost']:
            self.finish(False, 'DECELERATED: BEAM LOST')
            return
        self.update_phase()
        self.phi = math.remainder(self.phi, 2 * math.pi)
        self.arcs.append(self.arc(self.arcs[-1][1]))

    def advance_gap(self, dt):
        self.frac += dt / self.p['crossing_time']
        while self.frac >= 1 and not self.done:
            self.frac -= 1
            self.cross()
        if not self.done and self.w >= self.p['w_target']:
            self.finish(True, 'TARGET ENERGY')
        elif not self.done and self.time >= self.p['time_limit'] and not self.tutorial:
            self.finish(False, 'TIME UP')

    def observation(self):
        return {'w_MeV': self.w, 'phi_rad': self.phi, 'crossings': self.crossings, 'combo': self.combo}

    def metrics(self):
        return {'w_MeV': self.w, 'crossings': self.crossings, 'crest_fraction': self.crest / max(self.crossings, 1)}


class DeeBeat(GapStage):
    UNITS = {'b': 'T', 'v_gap': 'MV', 'w_start': 'MeV', 'w_target': 'MeV', 'w_lost': 'MeV', 'f_offset_min': '1',
             'f_offset_max': '1', 'f_rate': '1/s', 'crossing_time': 's', 'crest_cos': '1', 'time_limit': 's'}

    def field(self):
        return self.p['b']

    def setup(self):
        p = self.p
        self.f_nominal = physics.revolution_frequency(p['w_start'], p['b'])
        sign = self.rng.choice((-1, 1))
        self.f = self.f_nominal * (1 + sign * self.rng.uniform(p['f_offset_min'], p['f_offset_max']))
        self.setup_gap(self.rng.uniform(-math.pi / 2, math.pi / 2))

    def update_phase(self):
        self.phi += math.pi * (self.f / physics.revolution_frequency(self.w, self.p['b']) - 1)

    def advance(self, action, dt):
        self.f += action.x * self.p['f_rate'] * self.f_nominal * dt
        self.advance_gap(dt)

    def observation(self):
        return {**super().observation(), 'f_rf_Hz': self.f}


class TrimCoils(GapStage):
    UNITS = {'b': 'T', 'v_gap': 'MV', 'w_start': 'MeV', 'w_target': 'MeV', 'w_lost': 'MeV', 'delta_rate': '1/s',
             'crossing_time': 's', 'crest_cos': '1', 'time_limit': 's'}

    def field(self):
        return self.p['b'] * (1 + self.delta)

    def setup(self):
        self.delta, self.ghost_phi = 0.0, 0.0
        self.f = physics.revolution_frequency(self.p['w_start'], self.p['b'])
        self.setup_gap(0.0)

    def update_phase(self):
        self.phi += math.pi * (self.f / physics.revolution_frequency(self.w, self.field()) - 1)
        self.ghost_phi += math.pi * (self.f / physics.revolution_frequency(self.w, self.p['b']) - 1)

    def needed_delta(self):
        """B = gamma B0, relative to the injection energy, for which the RF is set."""
        return physics.gamma(self.w) / physics.gamma(self.p['w_start']) - 1

    def advance(self, action, dt):
        self.delta += action.y * self.p['delta_rate'] * dt
        self.advance_gap(dt)

    def observation(self):
        return {**super().observation(), 'delta': self.delta, 'delta_needed': self.needed_delta()}


class Fallout(Stage):
    UNITS = {'w_start': 'MeV', 'w_per_turn': 'MeV', 'turn_time': 's', 'gap': 'mm', 'z_start': 'mm', 'kick': 'mm/rad',
             'charges': '1', 'turns_target': '1'}

    def setup(self):
        p = self.p
        self.z = p['z_start'] * self.rng.choice((-1, 1))
        self.dz, self.turns, self.charges = 0.0, 0.0, p['charges']
        self.w = p['w_start']
        self.trail = []

    def advance(self, action, dt):
        p = self.p
        if self.y_edge and self.charges > 0:
            self.dz += math.copysign(p['kick'], self.y_edge)
            self.charges -= 1
        dtheta = 2 * math.pi * dt / p['turn_time']
        self.w = p['w_start'] + p['w_per_turn'] * self.turns
        s = math.sqrt(physics.isochronous_index(self.w))    # z'' = k z, exact over the step
        c, h = math.cosh(s * dtheta), math.sinh(s * dtheta)
        self.z, self.dz = c * self.z + h * self.dz / s, s * h * self.z + c * self.dz
        self.turns += dt / p['turn_time']
        self.points = int(1000 * self.turns)
        self.trail.append(self.z)
        del self.trail[:-600]
        if abs(self.z) >= p['gap']:
            self.finish(False, f'HIT THE POLE AFTER {self.turns:.1f} TURNS')
        elif self.turns >= p['turns_target']:
            self.finish(True, 'HELD ON. NOT FOR LONG')

    def observation(self):
        return {'z_mm': self.z, 'dz_mm_rad': self.dz, 'w_MeV': self.w, 'turns': self.turns, 'charges': self.charges}

    def metrics(self):
        return {'turns': self.turns, 'charges_left': self.charges}


class Thomas(Stage):
    UNITS = {'n_sectors': '1', 'w_start': 'MeV', 'w_target': 'MeV', 'w_per_turn': 'MeV', 'turn_time': 's', 'gap': 'mm',
             'z_start': 'mm', 'f_max': '1', 'f_start': '1', 'f_rate': '1/s', 'xi_max': 'deg', 'xi_rate': 'deg/s',
             'gradient_error': '1', 'dipole_error': 'mm', 'band_low': '1', 'band_high': '1'}
    DTHETA = 0.02    # rad, leapfrog step of the vertical motion

    def setup(self):
        p = self.p
        self.flutter, self.xi_deg, self.w, self.turns = p['f_start'], 0.0, p['w_start'], 0.0
        self.z, self.dz, self.theta = p['z_start'] * self.rng.choice((-1, 1)), 0.0, 0.0
        self.trail, self.in_band = [], 0.0
        self.update_tunes()

    def update_tunes(self):
        k = physics.isochronous_index(self.w)
        spiral = math.radians(self.xi_deg)
        self.qz2 = physics.qz_squared(k, self.p['n_sectors'], self.flutter, spiral)
        self.qr2 = physics.qr_squared(k, self.p['n_sectors'], self.flutter, spiral)
        self.qz = math.sqrt(self.qz2) if self.qz2 > 0 else 0.0

    def advance(self, action, dt):
        p = self.p
        self.flutter = min(max(self.flutter + action.y * p['f_rate'] * dt, 0.0), p['f_max'])
        self.xi_deg = min(max(self.xi_deg + action.x * p['xi_rate'] * dt, 0.0), p['xi_max'])
        self.w = p['w_start'] + p['w_per_turn'] * self.turns
        self.update_tunes()
        total = 2 * math.pi * dt / p['turn_time']
        n = max(1, math.ceil(total / self.DTHETA))
        h = total / n
        for _ in range(n):    # z'' + (Qz^2 + e cos theta) z = d cos theta: the errors act at 2 Qz = 1 and at Qz = 1
            self.z += 0.5 * h * self.dz
            self.theta += h
            c = math.cos(self.theta)
            self.dz += h * (-(self.qz2 + p['gradient_error'] * c) * self.z + p['dipole_error'] * c)
            self.z += 0.5 * h * self.dz
        self.turns += dt / p['turn_time']
        if p['band_low'] <= self.qz <= p['band_high']:
            self.in_band += dt / p['turn_time']
        self.points = int(500 * self.in_band + 10 * (self.w - p['w_start']))
        self.trail.append(self.z)
        del self.trail[:-600]
        if abs(self.z) >= p['gap']:
            self.finish(False, f'HIT THE POLE AT {self.w:.0f} MEV')
        elif self.w >= p['w_target']:
            self.finish(True, 'TARGET ENERGY')

    def observation(self):
        return {'z_mm': self.z, 'dz_mm_rad': self.dz, 'w_MeV': self.w, 'flutter': self.flutter, 'spiral_deg': self.xi_deg,
                'qz2': self.qz2, 'qr2': self.qr2}

    def metrics(self):
        return {'w_MeV': self.w, 'turns_in_band': self.in_band, 'turns': self.turns}


STAGES = {'larmor': Larmor, 'dee_beat': DeeBeat, 'trim_coils': TrimCoils, 'fallout': Fallout, 'thomas': Thomas}
