"""Gymnasium environment on the missions: tune a machine from control-room readings until the client's spec is met.

Observation: knob positions and readings, never the hidden errors. Action: bounded increments of every knob, at most one step each,
as an operator turns them. Reward: progress of Mission.tune_score, a cost per action, a bonus when in spec.
"""
import math

import gymnasium as gym
import numpy as np

from lions.game.missions import BUDGET, STEP_S
from lions.game.stages import load_scenario

ACTION_COST = 0.002
SPEC_BONUS = 1.0


def mission_by_id(mission_id):
    return next(m for m in load_scenario('missions') if m.mission_id == mission_id)


class MissionEnv(gym.Env):
    metadata = {'render_modes': []}

    def __init__(self, mission_id, seeds=None):
        """seeds: the conditions drawn from at reset. Training and evaluation use separate sets."""
        self.m = mission_by_id(mission_id)
        self.seeds = seeds
        self.steps = 0
        n = len(self.m.specs)
        self.action_space = gym.spaces.Box(-1.0, 1.0, (n,), np.float32)
        self.observation_space = gym.spaces.Box(-10.0, 10.0, (n + 9,), np.float32)

    def observe(self):
        m, r = self.m, self.m.readings
        knobs = [2 * (m.knobs[s['key']] - s['min']) / (s['max'] - s['min']) - 1 for s in m.specs]
        beam = m.has_beam()
        seen = np.flatnonzero(m.beam.probe > 0)
        phase = m.beam.phase[seen] if len(seen) else np.zeros(1)
        readings = [float(beam), r['reach'], r['transmission'], r['w_MeV'] / m.p['w_target'] - 1 if beam else 0.0,
                    math.tanh(r['x_mm'] / 20) if beam else 0.0, math.tanh(r['y_mm'] / 20) if beam else 0.0,
                    r['lost_kW'] / m.machine.max_lost_kw, float(phase.mean()) / 90, float(phase[-1] - phase[0]) / 90]
        return np.clip(np.array(knobs + readings, np.float32), -10, 10)

    def reset(self, *, seed=None, options=None):
        super().reset(seed=seed)
        drawn = int(self.np_random.choice(self.seeds)) if self.seeds is not None else int(self.np_random.integers(2**31))
        self.m.reset(drawn)
        self.steps, self.score = 0, self.m.tune_score()
        return self.observe(), {'condition': drawn}

    def act(self, action):
        m = self.m
        for spec, a in zip(m.specs, np.clip(action, -1, 1)):
            m.knobs[spec['key']] = m.clip(spec, m.knobs[spec['key']] + float(a) * spec['rate'] * STEP_S)

    def step(self, action):
        m = self.m
        action = np.clip(action, -1, 1)
        self.act(action)
        m.wait(STEP_S)
        m.refresh()
        self.steps += 1
        score = m.tune_score()
        reward = score - self.score - ACTION_COST * float(np.abs(action).mean())
        self.score = score
        terminated = bool(m.ready())
        if terminated:
            reward += SPEC_BONUS
        return self.observe(), float(reward), terminated, self.steps >= BUDGET, {'in_spec': terminated, 'score': score}
