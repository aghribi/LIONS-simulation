"""Post-hoc diagnostic of case T, not the registered metric. After the first read in spec: hold the knobs 50 passes, read without noise, replay the actions blind.

usage: uv run --extra game --extra rl --group bench python benchmarks/case_t_check.py
"""
import multiprocessing, sys
from pathlib import Path
import numpy as np, pyarrow as pa, pyarrow.parquet as pq
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'benchmarks'))
import case_t as T
from lions.game.env import MissionEnv, mission_by_id
OUT = ROOT/'benchmarks/runs/case_t'
HOLD = 50

def verdict(m):
    noise = m.p['noise']
    hold = 0
    for _ in range(HOLD):
        m.refresh(); hold += bool(m.ready())
    m.p['noise'] = 0.0
    m.refresh(); truth = bool(m.ready()); terms = {k: float(v) for k, v in m.spec_terms().items()}
    m.p['noise'] = noise
    worst = max(terms, key=terms.get)
    return hold/HOLD, truth, worst, terms[worst]

def classical(job):
    mission, method, seed = job
    m = mission_by_id(mission); m.reset(seed); rng = np.random.default_rng(seed)
    if method == 'coach':
        reached = T.coach_run(m)
    else:
        f = T.Counter(m)
        try: {'nelder_mead': T.nelder_mead, 'random': T.random_search, 'bayesian': T.bayesian}[method](f, rng)
        except T.Counter.Done: pass
        reached = f.reached
    row = dict(mission=mission, method=method, training_seed=-1, condition=seed, in_spec=reached is not None, evals=float(reached or 'nan'))
    row.update(zip(('hold','truth','worst','worst_term'), verdict(m)) if reached else dict(hold=float('nan'), truth=False, worst='', worst_term=float('nan')))
    return row

def rl(job):
    from stable_baselines3 import PPO, SAC
    import torch; torch.set_num_threads(1)
    mission, algo, seed, conds = job
    model = {'ppo': PPO, 'sac': SAC}[algo].load(OUT/f'{mission}_{algo}_{seed}.zip', device='cpu')
    rows, replay = [], None
    env = MissionEnv(mission)
    for c in conds:
        env.seeds = [c]; obs, _ = env.reset(seed=c); done = False; reached = None; acts = []
        start = {s['key']: env.m.knobs[s['key']] for s in env.m.specs}
        while not done:
            a, _ = model.predict(obs, deterministic=True); acts.append(a.copy())
            obs, _, term, trunc, _ = env.step(a); done = term or trunc
            if term: reached = env.steps
        row = dict(mission=mission, method=algo, training_seed=seed, condition=c, in_spec=reached is not None, evals=float(reached or 'nan'))
        if reached:
            ideal = env.m.ideal()
            row['moved'] = float(np.mean([abs(env.m.knobs[s['key']]-start[s['key']])/s['start_error'] for s in env.m.specs]))
            row['off_ideal'] = float(np.max([abs(env.m.knobs[s['key']]-ideal[s['key']])/s['start_error'] for s in env.m.specs]))
            row.update(zip(('hold','truth','worst','worst_term'), verdict(env.m)))
            if replay is None: replay = acts
        else:
            row.update(moved=float('nan'), off_ideal=float('nan'), hold=float('nan'), truth=False, worst='', worst_term=float('nan'))
        rows.append(row)
    # open loop: the first successful action sequence, blind, on every condition
    if replay is not None:
        for r, c in zip(rows, conds):
            env.seeds = [c]; env.reset(seed=c); ok = False
            for a in replay:
                _, _, term, _, _ = env.step(a)
                if term: ok = True; break
            r['open_loop'] = ok
    else:
        for r in rows: r['open_loop'] = False
    return rows

if __name__ == '__main__':
    fresh = list(range(20000, 20020))
    with multiprocessing.get_context('spawn').Pool(6) as pool:
        rows = pool.map(classical, [(m, a, s) for m in T.ALL for a in ('coach','nelder_mead','random','bayesian') for s in T.HELD_OUT], chunksize=1)
        pq.write_table(pa.Table.from_pylist(rows), OUT/'check_classical.parquet')
        print('classical done', flush=True)
        jobs = [(m, a, s, T.HELD_OUT) for m in ('F18','EYE') for a in ('ppo',) for s in T.RL_SEEDS] + [('F18','ppo',s,fresh) for s in T.RL_SEEDS]
        parts = pool.map(rl, jobs, chunksize=1)
        for i, p in enumerate(parts):
            for r in p: r['set'] = 'fresh' if r['condition'] >= 20000 else 'held_out'
        pq.write_table(pa.Table.from_pylist([r for p in parts for r in p]), OUT/'check_rl.parquet')
    print('done')
