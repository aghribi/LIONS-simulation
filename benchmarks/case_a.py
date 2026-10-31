"""Case A: uniform Bz, no E, one turn of a 1 MeV proton in LIONS, OPAL-cycl and RF-Track.

usage: uv run --group bench python benchmarks/case_a.py [outdir]
Writes case_a.parquet (one row per code, integrator, resolution and metric).
OPAL runs in the `lions-bench/opal:2024.1.0` image (benchmarks/opal/Dockerfile); outdir must be
visible to Docker (under $HOME with colima).
"""
import contextlib
import importlib.metadata
import io
import platform
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'benchmarks'))

import lions  # noqa: E402
import synthetic as syn  # noqa: E402
from provenance import provenance  # noqa: E402
import convert  # noqa: E402

OPAL_IMAGE = 'lions-bench/opal:2024.1.0'
W0 = 1.0          # MeV/u
BZ = -1.0         # T, LIONS sign for a proton
STEPS = [90, 180, 360, 720, 1440]
RFT_TOLERANCES = [1e-6, 1e-9, 1e-12]


def fit_circle(xy):
    """Least-squares circle through points (n, 2): centre and radius."""
    a = np.column_stack([xy, np.ones(len(xy))])
    sol, *_ = np.linalg.lstsq(a, (xy**2).sum(axis=1), rcond=None)
    centre = sol[:2] / 2
    return centre, np.sqrt(sol[2] + centre @ centre)


def write_inputs(n):
    syn.write_pot_central('pot1.bin')
    syn.write_pot_outer('pot2.dat')
    syn.write_bmap('bmap.dat', BZ)
    syn.write_deck('deck.dat', fhf=syn.revolution_frequency(W0, BZ), r0=syn.larmor_radius(W0, BZ),
                   w0=W0, turns=2, npaspt=n)


def initial_state():
    d = lions.State()
    d.run.swap_bytes = 'n'
    with contextlib.redirect_stdout(io.StringIO()):
        lions.load(d, 'deck.dat')
    b = d.bunch
    pos = np.array([[b.x[0], b.y[0], b.z[0]]])
    mom = np.array([[b.px[0], b.py[0], b.pz[0]]])
    return d, pos, mom, b.mass[0], b.charge[0]


def metrics(xy, p_first, p_last, rho):
    """n positions one step apart: radius, |p| and phase errors over the n - 1 steps."""
    centre, radius = fit_circle(xy)
    n = len(xy)
    angle = np.unwrap(np.arctan2(*(xy - centre).T[::-1]))
    swept = abs(angle[-1] - angle[0])
    return {'radius_rel_error': (radius / rho - 1, '1'),
            'p_rel_change': (p_last / p_first - 1, '1'),
            'phase_rel_error': (swept / (2 * np.pi * (n - 1) / n) - 1, '1')}


def run_lions(n, rho):
    rows = []

    def record(d):
        b = d.bunch
        rows.append((b.x[0], b.y[0], b.p[0]))
        return d.run.step < n - 1  # on_step also runs at step 0: steps 0 .. n - 1

    with contextlib.redirect_stdout(io.StringIO()):
        lions.run('deck.dat', on_step=record)
    a = np.array(rows)
    m = metrics(a[:, :2], a[0, 2], a[-1, 2], rho)
    m['radius_vs_boris_polygon'] = (fit_circle(a[:, :2])[1] / (rho * np.hypot(1, np.pi / n)) - 1, '1')
    return m


def opal_version():
    """OPAL version and git revision, as printed by the image."""
    log = subprocess.run(['docker', 'run', '--rm', OPAL_IMAGE, 'opal', '--version-full'],
                         capture_output=True, text=True).stdout
    version = re.search(r'Version (\S+)', log)
    rev = re.search(r'git rev\. (\w+)', log)
    if not (version and rev):
        raise RuntimeError(f'cannot read the OPAL version from {OPAL_IMAGE}')
    return f'{version.group(1)} ({rev.group(1)})'


def run_opal(n, d, pos, mom, mass_kg, workdir):
    mass_gev = mass_kg * syn.C**2 / syn.QE / 1e9
    pc_gev = np.linalg.norm(mom[0]) * mass_gev
    f_rev_mhz = d.rf.frequency / d.rf.harmonic / 1e6
    convert.write_opal_carbon(workdir / 'bmap_carbon.dat', d.field)
    (workdir / 'part.dat').write_text('1\n0 0 0 0 0 0\n')
    name = f'a_{{}}_{n}'
    out = {}
    for integrator in ('LF2', 'RK4'):
        stem = name.format(integrator)
        (workdir / f'{stem}.in').write_text(f'''OPTION, ECHO=FALSE;
OPTION, PSDUMPFREQ=100000000;
OPTION, SPTDUMPFREQ=1;
OPTION, VERSION=20240;
Title, string="LIONS benchmark case A";
cyc: CYCLOTRON, TYPE=CARBONCYCL, CYHARMON=1, PHIINIT=0.0, PRINIT=0.0,
     RINIT={np.hypot(*pos[0, :2]) * 1e3:.17g}, SYMMETRY=1.0, RFFREQ={f_rev_mhz:.17g},
     FMAPFN="bmap_carbon.dat", MINZ=-1000, MAXZ=1000, MINR=0, MAXR=10000;
l1: LINE=(cyc);
dist: DISTRIBUTION, TYPE=FROMFILE, FNAME="part.dat";
fs: FIELDSOLVER, FSTYPE=NONE, MX=8, MY=8, MT=8, PARFFTX=true, PARFFTY=true, PARFFTT=false,
    BCFFTX=open, BCFFTY=open, BCFFTT=open;
beam1: BEAM, MASS={mass_gev:.17g}, CHARGE=1.0, PC={pc_gev:.17g}, NPART=1, BCURRENT=0,
       BFREQ={f_rev_mhz:.17g};
TRACK, LINE=l1, BEAM=beam1, MAXSTEPS={n}, STEPSPERTURN={n}, TIMEINTEGRATOR="{integrator}";
RUN, METHOD="CYCLOTRON-T", BEAM=beam1, FIELDSOLVER=fs, DISTRIBUTION=dist;
ENDTRACK;
STOP;
''')
        subprocess.run(['docker', 'run', '--rm', '-v', f'{workdir}:/data', OPAL_IMAGE,
                        'bash', '-c', f'opal {stem}.in > {stem}.log 2>&1'], check=True)
        r, p = convert.read_opal_orbit(workdir / f'{stem}-trackOrbit.dat')
        r, p = r[:n], p[:n]
        rho = np.linalg.norm(p[0]) * mass_kg * syn.C / (syn.QE * abs(BZ))
        out[integrator] = metrics(r[:, :2], np.linalg.norm(p[0]), np.linalg.norm(p[-1]), rho)
    return out


def run_rftrack(n, pos, mom, mass_kg):
    import RF_Track as rft
    mass_mev = mass_kg * syn.C**2 / syn.QE / 1e6
    length = 1.0
    field = rft.ConstantField(length)
    field.set_static_Bfield(*convert.to_rftrack_field((0.0, 0.0, BZ)))
    volume = rft.Volume()
    volume.add(field, 0.0, 0.0, -length / 2)
    x0 = convert.to_rftrack(pos, mom, mass_mev)
    p_mev = np.linalg.norm(mom[0]) * mass_mev
    rho = p_mev * 1e6 / (syn.C * abs(BZ))
    period_mm = 2 * np.pi * rho * np.hypot(p_mev, mass_mev) / p_mev * 1e3
    runs = [('analytic', None)] + [(alg, tol) for alg in ('rkf45', 'rk2') for tol in RFT_TOLERANCES]
    out = {}
    for alg, tol in runs:
        opt = rft.TrackingOptions()
        opt.dt_mm = period_mm / n
        opt.tt_dt_mm = period_mm / n
        opt.t_max_mm = period_mm * (n - 1) / n
        opt.odeint_algorithm = alg
        if tol is not None:
            opt.odeint_epsabs = tol
            opt.odeint_epsrel = tol
        bunch = volume.track(rft.Bunch6dT(mass_mev, 1.0, 1.0, x0), opt)
        table = volume.get_transport_table('%mean_X %mean_Y %mean_Z')
        p_end = bunch.get_phase_space('%P')[0, 0]
        pos_l, _ = convert.from_rftrack(np.column_stack([table[:, 0], np.zeros(len(table)), table[:, 1],
                                                         np.zeros(len(table)), table[:, 2],
                                                         np.zeros(len(table))]), mass_mev)
        out[(alg, tol)] = metrics(pos_l[:n, :2], p_mev, p_end, rho)
    return out


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    import RF_Track as rft
    versions = {'lions': importlib.metadata.version('lions'),
                'opal': opal_version(), 'rf_track': str(rft.RF_TRACK_VERSION)}
    records = []

    def add(code, integrator, n, tol, m):
        for metric, (value, unit) in m.items():
            records.append({'code': code, 'version': versions[code], 'integrator': integrator,
                            'steps_per_turn': n, 'tolerance': tol, 'metric': metric,
                            'value': float(value), 'unit': unit})

    for n in STEPS:
        work = outdir / f'n{n}'
        work.mkdir(exist_ok=True)
        with contextlib.chdir(work):
            write_inputs(n)
            d, pos, mom, mass, charge = initial_state()
            rho = np.linalg.norm(mom[0]) * mass * syn.C / (charge * abs(BZ))
            add('lions', 'boris', n, None, run_lions(n, rho))
            for integ, m in run_opal(n, d, pos, mom, mass, work).items():
                add('opal', integ, n, None, m)
            for (alg, tol), m in run_rftrack(n, pos, mom, mass).items():
                add('rf_track', alg, n, tol, m)
        print(f'n={n} done', file=sys.stderr)

    table = pa.Table.from_pylist(records).replace_schema_metadata({
        'case': 'A: uniform Bz = -1 T, 1 MeV proton, one turn, no E',
        'units': 'value: dimensionless relative quantities (unit column)',
        **provenance(), 'machine': platform.platform(), 'python': platform.python_version(),
        **{f'version_{k}': v for k, v in versions.items()}})
    pq.write_table(table, outdir / 'case_a.parquet')
    return table


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'benchmarks' / 'runs' / 'case_a'
    for row in main(out.resolve()).to_pylist():
        print(row)
