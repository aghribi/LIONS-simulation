"""Case G3: OPALX 1.0.0 on the PSI Ring map (DESIGN.md section 25).

usage: uv run --group bench [--with pandas] python benchmarks/case_g3.py --step a|b|s|report [--outroot DIR]
a: closed orbits and tunes at the 24 G2 step 3 energies; b: det(M_r) against steps and FDSTEP at 3 energies;
s: SpectralTunes at the 24 energies; 2: LIONS and spline multi-turn tunes against amplitude (step 2); 2b: linear range of
the LIONS Jacobian; 2c: the same on a commensurate turn; report: the DESIGN.md section 25 numbers from the tables (needs pandas).
OPALX is run as is from its image. The PSI files (no licence, internal use) stay in benchmarks/external/ (git-ignored).
"""
import argparse
import contextlib
import json
import math
import multiprocessing
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

import case_g1  # noqa: E402
import case_g2  # noqa: E402
import orbits  # noqa: E402
from case_g2 import F_RF, HARMONIC, PSI, STEPS, SYMMETRY  # noqa: E402
from provenance import provenance  # noqa: E402

IMAGE = 'ghcr.io/opalx-project/opalx@sha256:994394a56b0d8ddeeae26285703fe9fea5e51408ebc1c811c1a3e966b4ed82b4'
INTEGRATORS = ('BORIS', 'RK4', 'DOP853')
R_FRAME = 2.1314    # m, radius of the RING frame origin; x = r - R_FRAME
R_SECTOR = 3.3      # m, element pose radius (RMIN + RMAX) / 2
RESIDUAL_MAX = 1e-9   # m, as G1 step 2
R_ERROR_MAX = 1e-3    # relative, a converged orbit further from ic.dat counts as a non-solution
G3B_ROWS = (40, 110, 190)    # ic.dat rows of G2 step 2b: 152, 292, 452 MeV
G3B_STEPS = (720, 1440, 2880)
FDSTEPS = (1e-7, 1e-6, 1e-5)
TURNS, TUNESAMPLE = 100, 50
WORKERS = 6
AMPLITUDES = (1e-5, 1e-4, 1e-3, 5e-3)    # m
ENVELOPE_TURNS = 10
LINEAR_AMPLITUDES = (1e-7, 1e-6, 1e-5, 1e-4)    # m
LINEAR_CHECK = (1, 2, 5, 10, 20)
LINEAR_TURNS = 20


def lattice():
    sectors = []
    for k in range(SYMMETRY):
        a = 2 * math.pi * k / SYMMETRY
        sectors.append(f'SM{k}: CYCLOTRONSECTOR, FMAPFN="bfield.dat", SYM={SYMMETRY}, RMIN=1.9, RMAX=4.7, VMIN=-0.05, '
                       f'VMAX=0.05, X={R_SECTOR * math.cos(a) - R_FRAME!r}, Y=0, Z={R_SECTOR * math.sin(a)!r}, '
                       f'THETA={-a!r}, PHI=0, PSI=0;')
    names = ','.join(f'SM{k}' for k in range(SYMMETRY))
    return '\n'.join(sectors) + f'\nMYCYCL: RING=({names}), X={R_FRAME!r};\n'


def run_opalx(work, name, body):
    work.mkdir(parents=True, exist_ok=True)
    if not (work / 'bfield.dat').exists():
        shutil.copy(PSI / 'bfield.dat', work / 'bfield.dat')
    (work / 'launch.dat').write_text('1\nx px y py z pz\n0 0 0 0 0 0.4\n')
    (work / f'{name}.in').write_text(f'''OPTION, VERSION=10900, ECHO=FALSE;
TITLE, STRING="G3 {name}";
{lattice()}DIST0: DISTRIBUTION, TYPE=FROMFILE, FNAME="launch.dat", NPARTDIST=1;
ES0: EMISSIONSOURCE, DISTRIBUTION=DIST0;
SOURCES0: EMISSIONSOURCELIST=(ES0);
{body}
QUIT;
''')
    with open(work / f'{name}.log', 'w') as log:
        return subprocess.run(['docker', 'run', '--rm', '--user', f'{os.getuid()}:{os.getgid()}', '-v', f'{work}:/data',
                               IMAGE, f'{name}.in', '--info', '2'], stdout=log, stderr=subprocess.STDOUT).returncode


def cof(work, energy, r_ic, pr_ic, integrator, steps, fd):
    """One COF at the ic.dat energy and start: the parsed JSON."""
    name = f'cof_{energy:.3f}_{integrator}_{steps}_{fd:g}'
    out = work / f'{name}.json'
    out.unlink(missing_ok=True)
    body = f'''BEAM0: BEAM, PARTICLE=PROTON, NALLOC=1, BCHARGE=1.602176634e-19, SOURCES=SOURCES0, CHARGE=1,
    ENERGY=PMASS+{energy / 1000!r};
ic: COF, LINE=MYCYCL, BEAM=BEAM0, DT={HARMONIC / F_RF / steps!r}, MAXSTEPS={3 * steps}, MAXPATH=40,
    TIMEINTEGRATOR="{integrator}", X={r_ic / 1000 - R_FRAME!r}, PX={pr_ic / 1000!r}, Y=0, PY=0,
    FDSTEP={{{fd!r},{fd!r},{fd!r},{fd!r}}}, OUTPUT="{out.name}";'''
    code = run_opalx(work, name, body)
    if not out.exists():
        raise RuntimeError(f'no COF result (exit {code}, see {name}.log)')
    return json.loads(out.read_text())


def metrics(result, energy, r_ic, pr_ic, qr_ref, qz_ref):
    p, _, q_ref = case_g2.reference(energy, r_ic, pr_ic, qr_ref, qz_ref)
    u, m = np.array(result['coordinates']), np.array(result['matrix'])
    r_err = (u[0] + R_FRAME) / (r_ic / 1000) - 1
    residual = float(np.max(np.abs(np.array(result['residual'])[[0, 2]])))
    if residual > RESIDUAL_MAX or abs(r_err) > R_ERROR_MAX:
        raise RuntimeError(f'non-solution: residual {residual:.1e} m, r error {r_err:.1e}')
    qr, qz = orbits.linear_tunes(m)
    return {'r_eo_rel_error': (r_err, '1'), 'pr_eo_error_over_p': ((u[1] - pr_ic / 1000) / p, '1'),
            'qr': (qr, '1'), 'qz': (qz, '1'), 'qr_error': (qr - q_ref[0], '1'), 'qz_error': (qz - q_ref[1], '1'),
            'det_r_error': (np.linalg.det(m[:2, :2]) - 1, '1'), 'det_z_error': (np.linalg.det(m[2:, 2:]) - 1, '1'),
            'det_error': (np.linalg.det(m) - 1, '1'), 'residual': (residual, 'm'),
            'energy_drift': (result['relative_energy_drift'], '1')}


def fixpo():
    """FIXPO rows: E [MeV], r [mm], pr [1e-3 m0c], Qr, Qz."""
    ic, ref = np.loadtxt(PSI / 'ic.dat'), np.loadtxt(PSI / 'refsol.dat')
    return np.hstack([ic, ref[:, 1:]])


def write(records, outdir, name, case):
    meta = {'case': case, 'opalx_image': IMAGE, 'opalx_tag': 'OPALX-1.0.0', 'rf_hz': repr(F_RF),
            'harmonic': str(HARMONIC), 'trim_coil': 'none', 'mass': 'OPALX PMASS 938.27208816 MeV',
            **provenance(PSI / 'bfield.dat', PSI / 'ic.dat', PSI / 'refsol.dat'),
            'machine': platform.platform(), 'python': platform.python_version()}
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / name)


def sweep(outdir, runs, name, case):
    """runs: (row, integrator, steps, fd). Failures are recorded, as in G2."""
    records = []
    for row, integrator, steps, fd in runs:
        energy, r_ic, pr_ic, qr_ref, qz_ref = map(float, row)
        base = {'code': 'opalx', 'integrator': integrator, 'steps_per_turn': steps, 'fdstep': fd,
                'energy_mev': float(energy), 'qr_half_integer_band': bool(386 <= energy <= 430)}
        try:
            values = metrics(cof(outdir, energy, r_ic, pr_ic, integrator, steps, fd), energy, r_ic, pr_ic, qr_ref, qz_ref)
        except Exception as exc:
            records.append({**base, 'metric': 'failed', 'value': None, 'unit': f'{type(exc).__name__}: {exc}'[:200]})
            print(f'{integrator} {steps} {energy} MeV failed: {exc}', file=sys.stderr)
            continue
        records += [{**base, 'metric': k, 'value': float(v), 'unit': u} for k, (v, u) in values.items()]
        print(f"{integrator} {steps} {fd:g} {energy} MeV qr={values['qr_error'][0]:+.2e} qz={values['qz_error'][0]:+.2e} "
              f"det_r={values['det_r_error'][0]:+.2e}", file=sys.stderr)
        write(records, outdir, name, case)
    write(records, outdir, name, case)


def step_a(outdir):
    rows = fixpo()
    sweep(outdir, [(row, i, STEPS, 1e-6) for i in INTEGRATORS for row in [*rows[::10], rows[-1]]], 'case_g3_a.parquet',
          'G3a: OPALX COF on the PSI Ring map, 24 energies, 1440 steps')


def step_b(outdir):
    rows = fixpo()[list(G3B_ROWS)]
    runs = [(row, i, n, 1e-6) for i in INTEGRATORS for n in G3B_STEPS for row in rows]
    runs += [(row, 'BORIS', STEPS, fd) for fd in FDSTEPS if fd != 1e-6 for row in rows]
    sweep(outdir, runs, 'case_g3_b.parquet', 'G3b: OPALX det(M_r) against steps and FDSTEP')


def step_s(outdir):
    rows = fixpo()
    rows = np.vstack([rows[::10], rows[-1]])
    triples = ','.join(f'{e / 1000!r},{r / 1000!r},{pr / 1000!r}' for e, r, pr, _, _ in rows.tolist())
    code = run_opalx(outdir, 'tunes', f'''FS0: FIELDSOLVER, TYPE=NONE, NX=16, NY=16, NZ=16;
BEAM0: BEAM, PARTICLE=PROTON, NALLOC=1, BCHARGE=1.602176634e-19, SOURCES=SOURCES0, CHARGE=1, TUNEINITIAL={{{triples}}};
TRACK, LINE=MYCYCL, BEAM=BEAM0, MAXSTEPS={TURNS * STEPS}, DT={HARMONIC / F_RF / STEPS!r}, STEPSPERTURN={STEPS};
RUN, METHOD="PARALLEL", FIELDSOLVER=FS0, SPECTRALTUNES=TRUE, TURNS={TURNS}, TUNESAMPLE={TUNESAMPLE}, TUNEINTEGRATOR="RK4", TUNESECTOR="SM0";
ENDTRACK;''')
    print(f'exit {code}; see {outdir}/tunes-tunes.csv', file=sys.stderr)


def amplitude_run(args):
    """All codes at one energy: tunes and amplitude ratios against the launch amplitude."""
    row, outdir = args
    grids, attrs, species, mass_kg, spline = SHARED
    energy, r_ic, pr_ic, qr_ref, qz_ref = row
    p, t, q_ref = case_g2.reference(energy, r_ic, pr_ic, qr_ref, qz_ref)
    work = outdir / f'lions_{energy:.3f}'
    work.mkdir(exist_ok=True)
    with contextlib.chdir(work):
        case_g1.write_lions_map(grids, attrs, 1)
    f_rev = F_RF / HARMONIC
    trackers = {'lions': orbits.LionsTracker(work, energy, f_rev, STEPS, species).track,
                'lions_hermite': orbits.LionsTracker(work, energy, f_rev, STEPS, species, {'interpolation': 'hermite'}).track,
                'spline': orbits.BorisTracker(spline, mass_kg, case_g1.lc.E_CHARGE, p, f_rev, STEPS).track}
    records, tol = [], (r_ic / 1000, pr_ic / 1000)
    for code, track in trackers.items():
        base = {'code': code, 'energy_mev': float(energy), 'qr_half_integer_band': bool(386 <= energy <= 430)}
        try:
            co, res = orbits.closed_orbit(track, tol)
            if res > RESIDUAL_MAX:
                raise RuntimeError(f'closed orbit residual {res:.1e} m')
            states = np.array([[co[0] + a, co[1], a, 0.0] for a in AMPLITUDES])
            series = orbits.section_series(*track(states, TURNS + 1.2))
        except Exception as exc:
            records.append({**base, 'amplitude_m': None, 'metric': 'failed', 'value': None,
                            'unit': f'{type(exc).__name__}: {exc}'[:200]})
            print(f'{code} {energy} MeV failed: {exc}', file=sys.stderr)
            continue
        for a, tbt in zip(AMPLITUDES, series):
            r, z = tbt[:TURNS, 1] - co[0], tbt[:TURNS, 3]
            rms = lambda x: np.sqrt(np.mean((x - np.mean(x))**2))
            qr, qz = orbits.tune(r), orbits.tune(z)
            values = {'qr': qr, 'qz': qz, 'qr_error': qr - q_ref[0], 'qz_error': qz - q_ref[1], 'turns': float(len(tbt)),
                      'amp_ratio_r': rms(r[-ENVELOPE_TURNS:]) / rms(r[:ENVELOPE_TURNS]),
                      'amp_ratio_z': rms(z[-ENVELOPE_TURNS:]) / rms(z[:ENVELOPE_TURNS])}
            records += [{**base, 'amplitude_m': a, 'metric': k, 'value': float(v), 'unit': '1'} for k, v in values.items()]
    print(f'{energy} MeV done', file=sys.stderr)
    return records


def step_2(outdir):
    global SHARED
    grids, attrs, species, mass_kg, spline, _ = case_g2.setup()
    SHARED = grids, attrs, species, mass_kg, spline    # inherited by fork: the spline field does not pickle
    rows = fixpo()
    runs = [(row.tolist(), outdir) for row in [*rows[::10], rows[-1]]]
    with multiprocessing.get_context('fork').Pool(WORKERS) as pool:
        records = [r for part in pool.imap_unordered(amplitude_run, runs) for r in part]
    write(records, outdir, 'case_g3_2.parquet', 'G3 step 2: LIONS and spline multi-turn tunes against launch amplitude')


def step_2b(outdir, commensurate=False):
    grids, attrs, species, _, _, _ = case_g2.setup()
    work = outdir / 'lions'
    work.mkdir(exist_ok=True)
    with contextlib.chdir(work):
        case_g1.write_lions_map(grids, attrs, 1)
    records = []
    amplitudes = LINEAR_AMPLITUDES[:2] if commensurate else LINEAR_AMPLITUDES
    for energy, r_ic, pr_ic, _, _ in fixpo()[list(G3B_ROWS)].tolist():
        p = np.sqrt((1 + energy / case_g2.MP_MEV)**2 - 1)
        f_rev, turn_error = F_RF / HARMONIC, None
        track = orbits.LionsTracker(work, energy, f_rev, STEPS, species).track
        co, _ = orbits.closed_orbit(track, (r_ic / 1000, pr_ic / 1000))
        for _ in range(6 if commensurate else 0):
            t_turn = orbits.section_series(*track(co[None], 1.2))[0][0, 0]
            turn_error = t_turn * f_rev - 1
            f_rev = 1 / t_turn
            track = orbits.LionsTracker(work, energy, f_rev, STEPS, species).track
            co, _ = orbits.closed_orbit(track, co[:2])
        if commensurate:
            print(f'{energy} MeV: last turn-time correction {turn_error:.1e}, f_rev {f_rev!r}', file=sys.stderr)
        scale = np.array([1.0, co[0] / p])    # pr to m, as the Jacobian probes
        m = orbits.jacobian(track, co, 1e-7, 1e-7 * p / co[0])[:2, :2]
        if commensurate:
            records.append({'energy_mev': float(energy), 'amplitude_m': None, 'turn': None, 'unit': '1',
                            'metric': 'det_r_error', 'value': float(np.linalg.det(m) - 1)})
        states = np.array([[co[0] + a, co[1], 0.0, 0.0] for a in amplitudes])
        series = orbits.section_series(*track(states, LINEAR_TURNS + 1.2))
        for a, tbt in zip(amplitudes, series):
            u0 = np.array([a, 0.0])
            for n in LINEAR_CHECK:
                u = tbt[n - 1, 1:3] - co[:2]
                predicted = np.linalg.matrix_power(m, n) @ u0
                row = {'energy_mev': float(energy), 'amplitude_m': a, 'turn': n, 'unit': '1'}
                records += [{**row, 'metric': 'deviation', 'value': float(np.linalg.norm((u - predicted) * scale) / np.linalg.norm(u * scale))},
                            {**row, 'metric': 'growth', 'value': float(np.linalg.norm(u * scale) / a)},
                            {**row, 'metric': 'growth_predicted', 'value': float(np.linalg.norm(predicted * scale) / a)}]
    name, case = ('case_g3_2c.parquet', 'G3 step 2c: commensurate turn') if commensurate else \
        ('case_g3_2b.parquet', 'G3 step 2b: linear range of the LIONS one-turn Jacobian')
    write(records, outdir, name, case)


def report(outroot):
    import pandas as pd
    fold, g2 = case_g2.fold, outroot.parent / 'case_g2'
    interior = lambda e: (e >= 92) & (e <= 492) & ~((e >= 386) & (e <= 430))   # as the G2 table
    pivot = lambda t: t.pivot_table(index='energy_mev', columns='metric', values='value')
    a = pd.read_parquet(outroot / 'step1a' / 'case_g3_a.parquet')
    step3, step2 = pd.read_parquet(g2 / 'step3' / 'case_g2_step3.parquet'), pd.read_parquet(g2 / 'step2' / 'case_g2_step2.parquet')
    opal = pivot(step3[step3.code == 'opal'])
    opal = opal[abs(opal.index - 512) > 1]    # OPAL non-solution, G2 step 3
    lions = pivot(step2[step2.code == 'lions'])
    for name, d in a.groupby('integrator'):
        d = pivot(d)
        i = d[interior(d.index)]
        j = i.join(opal, rsuffix='_o', how='inner')
        print(name, 'failed', int(d.get('failed', pd.Series(dtype=float)).notna().sum()), 'interior', len(i),
              'Qr clipped at', d.index[d.qr >= 0.5 - 1e-12].tolist())
        for m in ('r_eo_rel_error', 'qr_error', 'qz_error', 'det_r_error', 'det_z_error'):
            print(f'  |{m}| median {i[m].abs().median():.2e} max {i[m].abs().max():.2e}')
        for m, dq in (('r', j.r_eo_rel_error - j.r_eo_rel_error_o), ('Qr', j.qr - fold(j.qr_o)), ('Qz', j.qz - fold(j.qz_o))):
            print(f'  minus OPAL 2024.1 {m} median {dq.abs().median():.2e} max {dq.abs().max():.2e} (n={len(j)})')
        k = d.join(lions, rsuffix='_l', how='inner')
        print(f'  det_r correlation with LIONS det, 24 energies: {np.corrcoef(k.det_r_error, k.det_error_l)[0, 1]:+.2f}')
    b = pd.read_parquet(outroot / 'step1b' / 'case_g3_b.parquet')
    print(b[b.metric == 'det_r_error'].pivot_table(index=['integrator', 'steps_per_turn', 'fdstep'], columns='energy_mev',
                                                   values='value').to_string(float_format='{:+.3e}'.format))
    s = pd.read_csv(outroot / 'step1s' / 'tunes-tunes.csv')
    ref = np.loadtxt(PSI / 'refsol.dat')
    row = [int(np.argmin(abs(ref[:, 0] - e))) for e in s.energy_MeV]
    scale = (STEPS * TURNS / TUNESAMPLE - 1) * TUNESAMPLE / STEPS / TURNS    # sampled span / 100 nominal turns
    dr, dz = s.nu_r / scale - ref[row, 1], s.nu_y / scale - ref[row, 2]
    print(f'spectral minus FIXPO (not folded): Qr median {dr.abs().median():.1e} max {dr.abs().max():.1e}, '
          f'{(dr.abs() <= 5e-3).sum()} of {len(s)} within 5e-3; Qz median {dz.abs().median():.1e}, '
          f'{(dz.abs() <= 5e-3).sum()} within 5e-3, misses at {s.energy_MeV[dz.abs() > 5e-3].tolist()}: {dz[dz.abs() > 5e-3].round(4).tolist()}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--step', choices=('a', 'b', 's', '2', '2b', '2c', 'report'), required=True)
    parser.add_argument('--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_g3')
    args = parser.parse_args()
    if args.step == 'report':
        sys.exit(report(args.outroot.resolve()))
    outdir = args.outroot.resolve() / (f'step{args.step}' if args.step.startswith('2') else f'step1{args.step}')
    outdir.mkdir(parents=True, exist_ok=True)
    {'a': step_a, 'b': step_b, 's': step_s, '2': step_2, '2b': step_2b, '2c': lambda d: step_2b(d, True)}[args.step](outdir)
