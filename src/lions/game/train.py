"""Trains the PPO policy of a mission for the arcade and tests it on the held-out conditions. Seed 0, never the best of several.

usage: uv run --extra game --extra rl python -m lions.game.train F18 [EYE ...]
"""

import json
import sys

from lions.game.env import mission_by_id
from lions.game.tuners import (
    HELD_OUT,
    RL_STEPS,
    policy,
    policy_path,
    train_policy,
    tune,
)


def main(missions, steps=None):
    for mission_id in missions:
        path = policy_path(mission_id)
        path.parent.mkdir(parents=True, exist_ok=True)
        progress = path.with_suffix('.progress')
        total = steps or RL_STEPS['ppo']
        model = train_policy(
            mission_id,
            'ppo',
            0,
            path,
            steps,
            lambda n: progress.write_text(f'{n / total:.3f}'),
        )
        runs = []
        for condition in HELD_OUT:
            m = mission_by_id(mission_id)
            m.reset(condition)
            runs.append(tune(m, policy, condition, model=model).reached)
        path.with_suffix('.json').write_text(
            json.dumps(
                {
                    'in_spec': sum(r is not None for r in runs),
                    'conditions': len(runs),
                    'steps': total,
                    'seed': 0,
                }
            )
        )
        progress.unlink(missing_ok=True)
        print(
            mission_id,
            sum(r is not None for r in runs),
            '/',
            len(runs),
            flush=True,
        )


if __name__ == '__main__':
    main(sys.argv[1:])
