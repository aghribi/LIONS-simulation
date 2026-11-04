"""Case E2: speed spike of the case E PyTorch port against plain LIONS
(DESIGN.md section 20). Throwaway.

usage: uv run --group gpu python benchmarks/case_e2.py [outdir]
[--device cpu|mps|cuda] Writes case_e2.parquet. MPS has no float64: only
the float32 variants run there.
"""

import argparse
import contextlib
import io
import platform
import statistics
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
from provenance import provenance  # noqa: E402
from fields_b import B0  # noqa: E402
from case_e import W, Port, initial, load, write_case  # noqa: E402

STEPS = 3600
PARTICLES = [1, 1000, 100_000]
BACKWARD_PARTICLES = 1000
REPEATS = 3
DEVICE = 'cpu'


def make_port(d, dtype, compiled):
    # int x float promotes to the default dtype; case_e sets float64,
    # which MPS cannot hold
    torch.set_default_dtype(dtype)
    port = Port(d)
    port.maps = [m.to(DEVICE, dtype) for m in port.maps]
    if compiled:
        # a fresh cache per variant: past the recompile limit dynamo
        # falls back to eager without saying so
        torch._dynamo.reset()
        port.step = torch.compile(port.step, dynamic=False)
    return port


def start(d, n, dtype):
    r0 = torch.linspace(0.14, 0.15, n, dtype=torch.float64)
    return tuple(
        v.to(DEVICE, dtype) for v in initial(r0, 0.0, d.bunch.p[0], n)
    )


def timed(fn):
    """Median wall time of REPEATS calls after one warm-up (the warm-up
    time is returned too).
    """

    def sync():
        if DEVICE == 'cuda':
            torch.cuda.synchronize()
        elif DEVICE == 'mps':
            torch.mps.synchronize()

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
        syn.write_deck(
            'deck_speed.dat',
            fhf=syn.revolution_frequency(W, B0),
            r0=syn.larmor_radius(W, B0),
            w0=W,
            turns=STEPS // 360,
            npaspt=360,
            nbrion=n,
            rmax1=0.001,
        )
        return timed(lambda: lions.run('deck_speed.dat'))


def main(outdir):
    # the accuracy section runs the CPU reference, then restores the
    # device
    global DEVICE
    outdir.mkdir(parents=True, exist_ok=True)
    records, inputs = [], []
    dtypes = (
        (torch.float32,) if DEVICE == 'mps' else (torch.float64, torch.float32)
    )

    def add(code, dtype, compiled, mode, n, run):
        """Time run() and record it; a variant that fails (a backend that
        cannot compile it) is recorded too.
        """
        try:
            wall, first = run()
            error = None
        except Exception as exc:
            wall = first = None
            error = f'{type(exc).__name__}: {str(exc).strip().splitlines()[0][:300]}'
            print(
                f'{code} {dtype} compiled={compiled} {mode} {n}: {error}',
                file=sys.stderr,
            )
        records.append(
            {
                'code': code,
                'dtype': dtype,
                'compiled': compiled,
                'mode': mode,
                'particles': n,
                'steps': STEPS,
                'wall_s': wall,
                'first_call_s': first,
                'particle_steps_per_s': n * STEPS / wall if wall else None,
                'metric': None,
                'value': None,
                'error': error,
            }
        )

    work = outdir / 'B'
    work.mkdir(exist_ok=True)
    with contextlib.chdir(work):
        write_case('B')
        inputs.append(work / 'bmap.dat')
        d = load(360)
        for n in PARTICLES:
            add('lions', 'float64', False, 'forward', n, lambda: lions_time(n))
            print(f'lions {n} done', file=sys.stderr)
        for dtype in dtypes:
            for compiled in (False, True):
                port = make_port(d, dtype, compiled)
                name = str(dtype).removeprefix('torch.')
                for n in PARTICLES:
                    state = start(d, n, dtype)
                    with torch.no_grad():
                        add(
                            'torch',
                            name,
                            compiled,
                            'forward',
                            n,
                            lambda: timed(lambda: port.run(state, 1.0, STEPS)),
                        )
                    print(
                        f'torch {name} compiled={compiled} {n} done',
                        file=sys.stderr,
                    )
                state = start(d, BACKWARD_PARTICLES, dtype)

                def fwd_bwd():
                    s = torch.tensor(
                        1.0, dtype=dtype, device=DEVICE, requires_grad=True
                    )
                    port.run(state, s, STEPS)[2].sum().backward()

                add(
                    'torch',
                    name,
                    compiled,
                    'forward+backward',
                    BACKWARD_PARTICLES,
                    lambda: timed(fwd_bwd),
                )

    # each device variant against eager float64 on the CPU, cases B0 and
    # B
    device = DEVICE
    variants = (
        [(torch.float32, False)]
        if device == 'mps'
        else [
            (torch.float64, False),
            (torch.float64, True),
            (torch.float32, False),
        ]
    )
    for case in ('B0', 'B'):
        cw = outdir / case
        cw.mkdir(exist_ok=True)
        with contextlib.chdir(cw):
            write_case(case)
            inputs.append(cw / 'bmap.dat')
            d = load(360)
            out = {}
            for dev, dtype, compiled in [('cpu', torch.float64, False)] + [
                (device, *v) for v in variants
            ]:
                DEVICE = dev
                with torch.no_grad():
                    r, _, phase = make_port(d, dtype, compiled).run(
                        start(d, BACKWARD_PARTICLES, dtype), 1.0, STEPS
                    )
                out[(dev, dtype, compiled)] = (
                    r.cpu().double(),
                    phase.cpu().double(),
                )
            DEVICE = device
            ref_r, ref_phase = out[('cpu', torch.float64, False)]
            for (dev, dtype, compiled), (r, phase) in out.items():
                if (dev, dtype, compiled) == ('cpu', torch.float64, False):
                    continue
                label = f"{dev}_{str(dtype).removeprefix('torch.')}{'_compiled' if compiled else ''}"
                for metric, value in (
                    (
                        f'{label}_r_max_rel',
                        ((r - ref_r).abs() / ref_r).max().item(),
                    ),
                    (
                        f'{label}_phase_max_abs_rad',
                        (phase - ref_phase).abs().max().item(),
                    ),
                ):
                    records.append(
                        {
                            'code': 'torch',
                            'dtype': None,
                            'compiled': None,
                            'mode': case,
                            'particles': BACKWARD_PARTICLES,
                            'steps': STEPS,
                            'wall_s': None,
                            'first_call_s': None,
                            'particle_steps_per_s': None,
                            'metric': metric,
                            'value': value,
                            'error': None,
                        }
                    )

    device_name = {
        'cuda': lambda: torch.cuda.get_device_name(),
        'mps': lambda: 'Apple GPU (MPS)',
        'cpu': lambda: platform.processor() or platform.machine(),
    }[DEVICE]()
    meta = {
        'case': f'E2: speed spike, {W} MeV protons, {STEPS} steps, device {DEVICE}',
        'device_name': device_name,
        'torch': torch.__version__,
        'cpu': platform.processor() or platform.machine(),
        'repeats': str(REPEATS),
        **provenance(*inputs),
        'machine': platform.platform(),
        'python': platform.python_version(),
    }
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta),
        outdir / 'case_e2.parquet',
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'outdir',
        nargs='?',
        type=Path,
        default=ROOT / 'benchmarks' / 'runs' / 'case_e2',
    )
    parser.add_argument(
        '--device', choices=('cpu', 'mps', 'cuda'), default='cpu'
    )
    args = parser.parse_args()
    DEVICE = args.device
    main(args.outdir.resolve())
