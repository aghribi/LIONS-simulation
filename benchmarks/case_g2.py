"""Case G2: PSI Ring (map S03AV) against FIXPO (DESIGN.md section 21).

usage: uv run --group bench python benchmarks/case_g2.py [--step 1] [--outroot DIR]
Step 1: conventions of the files. 1a: units of the file theta derivatives; 1b: closed orbits at theta = 0 with the
numpy Boris tracker, both directions, against ic.dat under both pr readings; 1c: OPAL closed-orbit finder with and
without trim coil TC1.
The PSI files (no licence, internal use) stay in benchmarks/external/ (git-ignored); OPAL gets a copy of the map in
the run directory.
"""
import argparse
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.constants import e, m_p, physical_constants
from scipy.interpolate import RectBivariateSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

import orbits  # noqa: E402
from provenance import provenance  # noqa: E402

PSI = ROOT / 'benchmarks' / 'external' / 'psi_ring_github_io_7ed8895'
MP_MEV = physical_constants['proton mass energy equivalent in MeV'][0]
F_RF, HARMONIC, SYMMETRY = 50.65e6, 6, 8
ENERGIES = [72, 150, 300, 450, 530]     # MeV, nearest ic.dat rows
STEPS = 1440
PAD = 5
FLOAT = re.compile(r'[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?')
TC1 = 'tc1: TRIMCOIL, TYPE="PSI-BFIELD-MIRRORED", RMIN=4350.0, RMAX=4470.0, BMAX=1.4e-3, SLPTC=0.6;'


def read_ring(path):
    """RING map as OPAL reads it (getFieldFromFile_Ring): r [m], theta step [deg], Bz and three theta derivatives.

    Arrays (nrad, ntet) in file units (kG, per file azimuth unit); the first data record is skipped as in OPAL.
    """
    text, pos = path.read_text(), 0

    def scan(width=None):
        # fscanf: skip whitespace, then %s (up to whitespace) or %16lE (the longest float within 16 characters).
        # The numbers are fixed width, so a negative value touches its neighbour.
        nonlocal pos
        while text[pos].isspace():
            pos += 1
        if width is None:
            end = pos
            while end < len(text) and not text[end].isspace():
                end += 1
        else:
            end = FLOAT.match(text, pos, pos + width).end()
        out, pos = text[pos:end], end
        return out

    rmin, dr, tmin, dt = (float(scan()) for _ in range(4))
    dr = 1 / -dr if dr < 0 else dr
    dt = 1 / -dt if dt < 0 else dt
    if tmin != 0:
        raise ValueError(f'{path}: theta_min {tmin}, OPAL needs 0')
    for _ in range(13):
        scan()
    nrad, ntet = int(scan()), int(scan())
    for _ in range(5):
        scan()
    lpar = int(scan())
    for _ in range(4):
        scan()
    for _ in range(lpar):
        scan(16)
    for _ in range(6):
        scan()
    while scan() != 'LREC=':
        pass
    for _ in range(5):
        scan()
    blocks = np.empty((4, nrad, ntet))
    for ir in range(nrad):
        if ir:
            for _ in range(6):
                scan()
        for k in range(4):
            blocks[k, ir] = [float(scan(16)) for _ in range(ntet)]
    if text[pos:].strip():
        raise ValueError(f'{path}: data left after the last radius')
    return (rmin + dr * np.arange(nrad)) * 1e-3, dt, blocks


def periodic_d5(f, h):
    """5-point central first derivative along axis 1, periodic."""
    return (np.roll(f, 2, 1) - 8 * np.roll(f, 1, 1) + 8 * np.roll(f, -1, 1) - np.roll(f, -2, 1)) / (12 * h)


def ring_field(r, dt_deg, bz_kg, direction):
    """Median-plane field b(x, y, z) in T for protons travelling counterclockwise (LIONS frame, Bz < 0).

    direction +1: the map azimuth increases along the motion; -1: it decreases. Quintic spline, periodic per sector.
    """
    ntet = bz_kg.shape[1]
    jj = np.arange(-PAD, ntet + PAD)
    spline = RectBivariateSpline(r, np.radians(dt_deg) * jj, -0.1 * bz_kg[:, jj % ntet], kx=5, ky=5)
    sector = 2 * np.pi / SYMMETRY

    def field(x, y, z):
        t = (direction * np.arctan2(y, x)) % sector
        bz = spline.ev(np.hypot(x, y), t)
        return np.zeros_like(bz), np.zeros_like(bz), bz

    return field


def ic_rows():
    rows = np.loadtxt(PSI / 'ic.dat')
    return {e_: rows[np.argmin(np.abs(rows[:, 0] - e_))] for e_ in ENERGIES}


def opal_closed_orbit(work, energy_mev, r_mm, pr_ic, trim):
    """OPAL 2024.1 closed-orbit finder (tune mode) at one energy: (r [m], pr [beta gamma], Qr, Qz)."""
    work.mkdir(parents=True, exist_ok=True)
    (work / 'data').mkdir(exist_ok=True)
    shutil.copy(PSI / 'bfield.dat', work / 'bfield.dat')
    gev = energy_mev / 1000
    # FMLOWE, FMHIGHE, DENERGY in GeV in OPAL 2024.1; the element cannot be named RING (protected keyword).
    (work / 'tune.in').write_text(f'''OPTION, ECHO=FALSE;
OPTION, CLOTUNEONLY=TRUE;
TITLE, STRING="PSI Ring closed orbit";
{TC1 if trim else ''}
psiring: CYCLOTRON, TYPE="RING", CYHARMON={HARMONIC}, PHIINIT=0.0, PRINIT={pr_ic / 1000:.17g}, RINIT={r_mm:.17g},
    SYMMETRY={SYMMETRY}.0, RFFREQ={F_RF / 1e6:.17g}, FMAPFN="bfield.dat", FMLOWE={gev:.17g}, FMHIGHE={gev + 0.002:.17g}
    {', TRIMCOIL=tc1' if trim else ''};
l1: LINE=(psiring);
dist: DISTRIBUTION, TYPE=GAUSSMATCHED, LINE=l1, NSTEPS={STEPS}, DENERGY=0.002, MAXSTEPSCO=100, NSECTORS={SYMMETRY},
    SECTOR=FALSE;
fs: FIELDSOLVER, FSTYPE=NONE, MX=8, MY=8, MT=8, PARFFTX=true, PARFFTY=true, PARFFTT=false,
    BCFFTX=open, BCFFTY=open, BCFFTT=open;
beam1: BEAM, PARTICLE=PROTON, PC={np.sqrt((gev + MP_MEV / 1000)**2 - (MP_MEV / 1000)**2):.17g}, NPART=1,
    BCURRENT=1.0E-6, CHARGE=1.0, BFREQ={F_RF / 1e6:.17g};
TRACK, LINE=l1, BEAM=beam1, MAXSTEPS=720, STEPSPERTURN=720;
RUN, METHOD="CYCLOTRON-T", BEAM=beam1, FIELDSOLVER=fs, DISTRIBUTION=dist;
ENDTRACK;
STOP;
''')
    subprocess.run(['docker', 'run', '--rm', '-v', f'{work.resolve()}:/data', orbits.OPAL_IMAGE,
                    'bash', '-c', 'opal tune.in > tune.log 2>&1'], check=True)
    rows = np.atleast_2d(np.loadtxt(work / 'data' / 'tunes.dat'))
    if not len(rows) or abs(rows[0, 0] - energy_mev) > 1e-6:
        raise RuntimeError(f'{work}: no OPAL closed orbit at {energy_mev} MeV (see tune.log)')
    _, r, pr, _, qr, qz = rows[0]
    return r, pr, qr, qz


def step1(outdir):
    r, dt_deg, (bz, dbt, _, _) = read_ring(PSI / 'bfield.dat')
    records = []

    # 1a: units of the file theta derivatives
    for unit, h in (('per_deg', dt_deg), ('per_rad', np.radians(dt_deg))):
        diff = np.abs(dbt - periodic_d5(bz, h))
        for stat, v in (('max', diff.max()), ('median', np.median(diff))):
            records.append({'run': '1a', 'energy_mev': None, 'case': unit, 'metric': f'dbz_dtheta_diff_{stat}',
                            'value': float(v), 'unit': 'kG per file unit'})
        print(f'1a {unit}: max {diff.max():.3e}, median {np.median(diff):.3e} (|dBz/dtheta| max {np.abs(dbt).max():.3e})',
              file=sys.stderr)

    # 1b: closed orbits with the numpy Boris tracker, both directions
    ic = ic_rows()
    for direction in (1, -1):
        field = ring_field(r, dt_deg, bz, direction)
        for energy, (e_ic, r_ic, pr_ic) in ic.items():
            gamma = 1 + e_ic / MP_MEV
            p = np.sqrt(gamma**2 - 1)
            readings = {'pr_1e-3_m0c': pr_ic / 1000, 'pr_mrad': pr_ic / 1000 * p}
            guess = (r_ic / 1000, direction * np.mean(list(readings.values())))
            tracker = orbits.BorisTracker(field, m_p, e, p, F_RF / HARMONIC, STEPS)
            case = f'direction{direction:+d}'
            try:
                s, res = orbits.closed_orbit(tracker.track, guess)
            except Exception as exc:
                records.append({'run': '1b', 'energy_mev': e_ic, 'case': case, 'metric': 'failed', 'value': None,
                                'unit': f'{type(exc).__name__}: {exc}'[:200]})
                print(f'1b {case} {e_ic} MeV failed: {exc}', file=sys.stderr)
                continue
            values = {'r_minus_ic_mm': (s[0] * 1000 - r_ic, 'mm'), 'residual': (res, 'm'),
                      'pr_beta_gamma': (s[1], '1'),
                      **{f'pr_rel_error_{k}': (s[1] / v - 1, '1') for k, v in readings.items()}}
            for k, (v, u) in values.items():
                records.append({'run': '1b', 'energy_mev': e_ic, 'case': case, 'metric': k, 'value': float(v),
                                'unit': u})
            print(f'1b {case} {e_ic} MeV ' + ' '.join(f'{k}={v[0]:+.3e}' for k, v in values.items()), file=sys.stderr)

    # 1c: OPAL closed-orbit finder, with and without TC1
    for energy in (72, 530):
        e_ic, r_ic, pr_ic = ic[energy]
        for trim in (False, True):
            case = 'tc1' if trim else 'no_trim'
            r_o, pr_o, qr, qz = opal_closed_orbit(outdir / f'opal_{energy}_{case}', e_ic, r_ic, pr_ic, trim)
            values = {'r_m': (r_o, 'm'), 'pr_beta_gamma': (pr_o, '1'), 'qr': (qr, '1'), 'qz': (qz, '1')}
            for k, (v, u) in values.items():
                records.append({'run': '1c', 'energy_mev': e_ic, 'case': case, 'metric': k, 'value': float(v),
                                'unit': u})
            print(f'1c {case} {e_ic} MeV ' + ' '.join(f'{k}={v[0]:.7g}' for k, v in values.items()), file=sys.stderr)

    meta = {'case': 'G2 step 1: PSI Ring file conventions', 'proton_mass_mev': repr(MP_MEV),
            'rf_hz': repr(F_RF), 'harmonic': str(HARMONIC), 'steps_per_turn': str(STEPS),
            'opal_image': orbits.OPAL_IMAGE,
            **provenance(PSI / 'bfield.dat', PSI / 'ic.dat'),
            'machine': platform.platform(), 'python': platform.python_version()}
    outdir.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(records).replace_schema_metadata(meta), outdir / 'case_g2_step1.parquet')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--step', type=int, choices=(1,), default=1)
    parser.add_argument('--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_g2')
    args = parser.parse_args()
    {1: step1}[args.step](args.outroot.resolve() / f'step{args.step}')
