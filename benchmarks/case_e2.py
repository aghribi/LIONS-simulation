"""Case E2: speed spike of the case E PyTorch port against plain LIONS (DESIGN.md section 20). Throwaway.

usage: uv run --group bench --with torch python benchmarks/case_e2.py [outdir]
Writes case_e2.parquet.
"""
import contextlib
import importlib.metadata
import io
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'benchmarks'))

import lions  # noqa: E402
import synthetic as syn  # noqa: E402
from case_b import B0  # noqa: E402
from case_e import W, Port, initial, load, write_case  # noqa: E402

STEPS = 3600
PARTICLES = [1, 1000, 100_000]
BACKWARD_PARTICLES = 1000
REPEATS = 3
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'


def make_port(d, dtype, compiled):
    port = Port(d)
    port.maps = [m.to(DEVICE, dtype) for m in port.maps]
    if compiled:
        port.step = torch.compile(port.step, dynamic=False)
    return port


def start(d, n, dtype):
    r0 = torch.linspace(0.14, 0.15, n, dtype=torch.float64)
    return tuple(v.to(DEVICE, dtype) for v in initial(r0, 0.0, d.bunch.p[0], n))


def timed(fn):
    """Median wall time of REPEATS calls after one warm-up (the warm-up time is returned too)."""
    def sync():
        if DEVICE == 'cuda':
            torch.cuda.synchronize()
    t0 = time.perf_counter()
    fn()
    sync()
    first = time.perf_counter() - t0
    times = []
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        fn()
        sync()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), first


def lions_time(n):
    with contextlib.redirect_stdout(io.StringIO()):
        syn.write_deck('deck_speed.dat', fhf=syn.revolution_frequency(W, B0), r0=syn.larmor_radius(W, B0), w0=W,
                       turns=STEPS // 360, npaspt=360, nbrion=n, rmax1=0.001)
        return timed(lambda: lions.run('deck_speed.dat'))


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True,
                            text=True).stdout.strip()
    records = []

    def add(code, dtype, compiled, mode, n, wall, first):
        records.append({'code': code, 'dtype': dtype, 'compiled': compiled, 'mode': mode, 'particles': n,
                        'steps': STEPS, 'wall_s': wall, 'first_call_s': first,
                        'particle_steps_per_s': n * STEPS / wall, 'metric': None, 'value': None})

    work = outdir / 'B'
    work.mkdir(exist_ok=True)
    with contextlib.chdir(work):
        write_case('B')
        d = load(360)
        for n in PARTICLES:
            add('lions', 'float64', False, 'forward', n, *lions_time(n))
            print(f'lions {n} done', file=sys.stderr)
        for dtype in (torch.float64, torch.float32):
            for compiled in (False, True):
                port = make_port(d, dtype, compiled)
                name = str(dtype).removeprefix('torch.')
                for n in PARTICLES:
                    state = start(d, n, dtype)
                    with torch.no_grad():
                        add('torch', name, compiled, 'forward', n,
                            *timed(lambda: port.run(state, 1.0, STEPS)))
                    print(f'torch {name} compiled={compiled} {n} done', file=sys.stderr)
                state = start(d, BACKWARD_PARTICLES, dtype)

                def fwd_bwd():
                    s = torch.tensor(1.0, dtype=dtype, device=DEVICE, requires_grad=True)
                    port.run(state, s, STEPS)[2].sum().backward()
                add('torch', name, compiled, 'forward+backward', BACKWARD_PARTICLES,
                    *timed(fwd_bwd))

    # float32 and compiled results against eager float64, cases B0 and B
    for case in ('B0', 'B'):
        cw = outdir / case
        cw.mkdir(exist_ok=True)
        with contextlib.chdir(cw):
            write_case(case)
            d = load(360)
            out = {}
            for dtype, compiled in ((torch.float64, False), (torch.float64, True), (torch.float32, False)):
                with torch.no_grad():
                    r, _, phase = make_port(d, dtype, compiled).run(start(d, BACKWARD_PARTICLES, dtype), 1.0, STEPS)
                out[(dtype, compiled)] = r.double().cpu(), phase.double().cpu()
            ref_r, ref_phase = out[(torch.float64, False)]
            for key, label in (((torch.float64, True), 'compiled_float64'), ((torch.float32, False), 'float32')):
                r, phase = out[key]
                for metric, value in ((f'{label}_r_max_rel', ((r - ref_r).abs() / ref_r).max().item()),
                                      (f'{label}_phase_max_abs_rad', (phase - ref_phase).abs().max().item())):
                    records.append({'code': 'torch', 'dtype': None, 'compiled': None, 'mode': case,
                                    'particles': BACKWARD_PARTICLES, 'steps': STEPS, 'wall_s': None,
                                    'first_call_s': None, 'particle_steps_per_s': None, 'metric': metric,
                                    'value': value})

    meta = {'case': f'E2: speed spike, {W} MeV protons, {STEPS} steps, device {DEVICE}',
            'torch': torch.__version__, 'lions': importlib.metadata.version('lions'),
            'cpu': platform.processor() or platform.machine(), 'repeats': str(REPEATS),
            'git_commit': commit, 'machine': platform.platform(), 'python': platform.python_version()}
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / 'case_e2.parquet')


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'benchmarks' / 'runs' / 'case_e2'
    main(out.resolve())
