"""Case E: differentiability spike, LIONS magnetic tracking ported to PyTorch (DESIGN.md section 17). Throwaway.

usage: uv run --group bench --with torch python benchmarks/case_e.py [outdir]
       uv run --group bench --with torch python benchmarks/case_e.py --cost N PARTICLES CHECKPOINT  (one cost point, JSON)
Writes case_e.parquet.
"""
import contextlib
import io
import json
import math
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import torch
from torch.utils.checkpoint import checkpoint

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'benchmarks'))

import lions  # noqa: E402
import synthetic as syn  # noqa: E402
from provenance import provenance  # noqa: E402
from fields_b import B0, MAP, write_field  # noqa: E402

torch.set_default_dtype(torch.float64)
W, Z0 = 1.0, 1e-3                       # MeV, m
STEPS = [360, 3600]
FD = {'s': 1e-6, 'r0': 1e-7, 'pr0': 1e-7}   # s relative, r0 in m, pr0 in beta gamma
SWEEP = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
COST_STEPS, COST_PARTICLES = [1_000, 10_000, 100_000], [1, 1000]
OUTPUTS = ('r', 'pr', 'phase')


def write_case(case):
    if case == 'A':
        syn.write_pot_central('pot1.bin')
        syn.write_pot_outer('pot2.dat')
        syn.write_bmap('bmap.dat', bz=B0, ntet=361, **MAP)
    else:
        write_field(case, 361)


def load(n_steps, turns=1):
    d = lions.State()
    d.run.swap_bytes = 'n'
    with contextlib.redirect_stdout(io.StringIO()):
        syn.write_deck('deck.dat', fhf=syn.revolution_frequency(W, B0), r0=syn.larmor_radius(W, B0), w0=W,
                       z0=Z0, turns=turns, npaspt=n_steps)
        lions.load(d, 'deck.dat')
    return d


class Port:
    """LIONS compute_fields (magnetic part, map region) and boris_push with E = 0."""

    def __init__(self, d):
        f = d.field
        self.maps = [torch.tensor(a) for a in (f.bz, f.d2bz, f.dbzdx, f.dbzdy)]
        self.rmin, self.dr, self.dtheta = f.rmin, f.dr, f.dtheta
        self.bcoef, self.c_dt = d.bunch.boris_b_coef[0], d.run.c_dt

    def fields(self, x, y, z, s):
        r = torch.sqrt(x * x + y * y)
        theta = torch.remainder(torch.atan2(y, x), 2 * math.pi)
        i = ((r - self.rmin) / self.dr).detach().long()
        j = (theta / self.dtheta).detach().long()
        a = r - self.rmin - i * self.dr
        b = theta - j * self.dtheta
        d1, d2 = a - self.dr, b - self.dtheta
        c1, c2, c3, c4 = d1 * d2, a * d2, a * b, d1 * b
        vz = z * z * 0.5
        cell = [[m[i + p, j + q] for p in (0, 1) for q in (0, 1)] for m in self.maps]
        (b00, b01, b10, b11), (e00, e01, e10, e11) = cell[0], cell[1]
        h = 1.0 / (self.dr * self.dtheta)
        bz = (c1 * (b00 - vz * e00) - c2 * (b10 - vz * e10) + c3 * (b11 - vz * e11) - c4 * (b01 - vz * e01)) * h
        bx, by = ((c1 * g00 - c2 * g10 + c3 * g11 - c4 * g01) * z * h for g00, g01, g10, g11 in cell[2:])
        return s * bx, s * by, s * bz

    def step(self, x, y, z, px, py, pz, s):
        bx, by, bz = self.fields(x, y, z, s)
        gamma = torch.sqrt(1.0 + px**2 + py**2 + pz**2)
        vec = self.bcoef / gamma
        a1, a2, a3 = vec * bx, vec * by, vec * bz
        v1, v2, v3 = a1 * a1, a2 * a2, a3 * a3
        vec = 1.0 / (1.0 + v1 + v2 + v3)
        s12 = a1 * a2
        s21, s12 = s12 - a3, s12 + a3
        s13 = a1 * a3
        s31, s13 = s13 + a2, s13 - a2
        s23 = a2 * a3
        s32, s23 = s23 - a1, s23 + a1
        v4, v5, v6 = vec * (0.0 + 2.0 * px), vec * (0.0 + 2.0 * py), vec * (0.0 + 2.0 * pz)
        px, py, pz = (px - (v2 + v3) * v4 + s12 * v5 + s13 * v6,
                      py + s21 * v4 - (v1 + v3) * v5 + s23 * v6,
                      pz + s31 * v4 + s32 * v5 - (v1 + v2) * v6)
        vec = self.c_dt / torch.sqrt(1.0 + px**2 + py**2 + pz**2)
        return x + vec * px, y + vec * py, z + vec * pz, px, py, pz

    def segment(self, n, s, phase, *state):
        for _ in range(n):
            new = self.step(*state, s)
            x0, y0, x1, y1 = state[0], state[1], new[0], new[1]
            phase = phase + torch.atan2(x0 * y1 - y0 * x1, x0 * x1 + y0 * y1)
            state = new
        return (phase, *state)

    def run(self, state, s, n, segment=None):
        """(r, pr, phase) after n steps; segment = steps per checkpointed block, None = no checkpointing."""
        phase = torch.zeros_like(state[0])
        if segment is None:
            phase, *state = self.segment(n, s, phase, *state)
        else:
            done = 0
            while done < n:
                k = min(segment, n - done)
                phase, *state = checkpoint(lambda *a, k=k: self.segment(k, s, *a), phase, *state,
                                           use_reentrant=False)
                done += k
        x, y = state[0], state[1]
        r = torch.sqrt(x * x + y * y)
        return r, (x * state[3] + y * state[4]) / r, phase


def initial(r0, pr0, p, n=1):
    zero = torch.zeros(n)
    pr0 = pr0 + zero
    return (r0 + zero, zero, Z0 + zero, pr0, torch.sqrt(p**2 - pr0**2), zero)


def validate(d):
    """Max relative difference between LIONS and the port over one turn, positions and momenta."""
    n = d.run.n_steps
    rows = []

    def record(st):
        b = st.bunch
        rows.append([b.x[0], b.y[0], b.z[0], b.px[0], b.py[0], b.pz[0]])
        return st.run.step < n

    with contextlib.redirect_stdout(io.StringIO()):
        lions.run('deck.dat', on_step=record)
    ref = np.array(rows)
    port = Port(d)
    state = tuple(torch.tensor([v]) for v in ref[0])
    out = []
    with torch.no_grad():
        for _ in range(n):
            state = port.step(*state, 1.0)
            out.append([v.item() for v in state])
    out = np.array(out)
    scale = np.repeat([np.abs(ref[1:, :3]).max(), np.abs(ref[1:, 3:]).max()], 3)   # positions, momenta
    return np.max(np.abs(out - ref[1:]) / scale)


def gradients(d, n):
    """Autograd and central finite differences of (r, pr, phase) after n steps w.r.t. (s, r0, pr0)."""
    port, p = Port(d), d.bunch.p[0]
    r0 = syn.larmor_radius(W, B0)
    params = {'s': 1.0, 'r0': r0, 'pr0': 0.0}

    def forward(v):
        return port.run(initial(v['r0'], v['pr0'], p), v['s'], n)

    leaf = {k: torch.tensor(v, requires_grad=True) for k, v in params.items()}
    out = forward(leaf)
    ad = {o: torch.autograd.grad(out[i].sum(), list(leaf.values()), retain_graph=True)
          for i, o in enumerate(OUTPUTS)}
    rows = []
    with torch.no_grad():
        for j, k in enumerate(params):
            for h_rel in sorted({FD[k], *SWEEP}, reverse=True):
                h = h_rel * (params[k] if k == 's' else 1.0)
                hi, lo = forward({**params, k: params[k] + h}), forward({**params, k: params[k] - h})
                for i, o in enumerate(OUTPUTS):
                    fd = ((hi[i] - lo[i]) / (2 * h)).item()
                    g = ad[o][j].item()
                    rows.append({'output': o, 'param': k, 'fd_step': h_rel, 'registered': h_rel == FD[k],
                                 'autograd': g, 'fd': fd, 'rel_diff': abs(g - fd) / abs(fd) if fd else math.inf})
    return rows


def cost_point(n, particles, ckpt):
    """Forward plus backward wall time and peak RSS for one configuration (run in its own process)."""
    with contextlib.chdir(ROOT / 'benchmarks' / 'runs' / 'case_e' / 'B'):
        d = load(360)
    port, p = Port(d), d.bunch.p[0]
    r0 = torch.linspace(0.14, 0.15, particles)
    s = torch.tensor(1.0, requires_grad=True)
    t0 = time.perf_counter()
    r, pr, phase = port.run(initial(r0, 0.0, p, particles), s, n, int(math.isqrt(n)) if ckpt else None)
    phase.sum().backward()
    wall = time.perf_counter() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {'steps': n, 'particles': particles, 'checkpoint': ckpt, 'wall_s': wall,
            'peak_rss_mb': rss / (1 if platform.system() == 'Linux' else 1024) / 1024}


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    records, cost = [], []
    for case in ('A', 'B0', 'B'):
        work = outdir / case
        work.mkdir(exist_ok=True)
        with contextlib.chdir(work):
            write_case(case)
            err = validate(load(STEPS[0]))
            records.append({'case': case, 'steps': STEPS[0], 'output': 'state', 'param': 'validation',
                            'fd_step': None, 'registered': True, 'autograd': None, 'fd': None, 'rel_diff': err})
            if err > 1e-13:
                raise RuntimeError(f'case {case}: port differs from LIONS by {err} relative')
            for n in STEPS:
                records += [{'case': case, 'steps': n, **row} for row in gradients(load(STEPS[0]), n)]
        print(f'case {case} done', file=sys.stderr)
    meta = {'case': f'E: PyTorch port of the LIONS Boris step, {W} MeV proton, z0 = {Z0} m',
            'torch': torch.__version__,
            **provenance(), 'machine': platform.platform(), 'python': platform.python_version()}
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / 'case_e.parquet')
    ram_mb = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / 2**20
    for n in COST_STEPS:
        for particles in COST_PARTICLES:
            for ckpt in (False, True):
                # memory grows linearly in steps without checkpointing: skip a point extrapolated past half the RAM
                prev = [c for c in cost if c['particles'] == particles and not c['checkpoint'] and c['peak_rss_mb']]
                if not ckpt and prev and prev[-1]['peak_rss_mb'] * n / prev[-1]['steps'] > ram_mb / 2:
                    cost.append({'steps': n, 'particles': particles, 'checkpoint': ckpt, 'wall_s': None,
                                 'peak_rss_mb': None,
                                 'error': f"skipped: extrapolated {prev[-1]['peak_rss_mb'] * n / prev[-1]['steps']:.0f} MB"
                                          f" > half of {ram_mb:.0f} MB RAM"})
                    continue
                res = subprocess.run([sys.executable, __file__, '--cost', str(n), str(particles), str(int(ckpt))],
                                     capture_output=True, text=True)
                # the RF-Track import prints a banner on stdout: the result is the last line
                cost.append({'error': None, **json.loads(res.stdout.strip().splitlines()[-1])} if res.returncode == 0 else
                            {'steps': n, 'particles': particles, 'checkpoint': ckpt, 'wall_s': None,
                             'peak_rss_mb': None, 'error': res.stderr[-300:]})
                print(cost[-1], file=sys.stderr)
    pq.write_table(pa.Table.from_pylist(cost).replace_schema_metadata(meta), outdir / 'case_e_cost.parquet')


if __name__ == '__main__':
    if sys.argv[1:2] == ['--cost']:
        print(json.dumps(cost_point(int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4])))))
    else:
        out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'benchmarks' / 'runs' / 'case_e'
        main(out.resolve())
