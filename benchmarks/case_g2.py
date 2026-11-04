"""Case G2: PSI Ring (map S03AV) against FIXPO (DESIGN.md section 21).

usage: uv run --group bench python benchmarks/case_g2.py [--step 1]
[--outroot DIR] Step 1: conventions of the files. 1a: units of the file
theta derivatives; 1b: closed orbits at theta = 0 with the numpy Boris
tracker, both directions, against ic.dat under both pr readings; 1c:
OPAL closed-orbit finder with and without trim coil TC1. The PSI files
(no licence, internal use) stay in benchmarks/external/ (git-ignored);
OPAL gets a copy of the map in the run directory.
"""

import argparse
import contextlib
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

import case_g1  # noqa: E402
import orbits  # noqa: E402
from provenance import provenance  # noqa: E402

PSI = ROOT / 'benchmarks' / 'external' / 'psi_ring_github_io_7ed8895'
MP_MEV = physical_constants['proton mass energy equivalent in MeV'][0]
F_RF, HARMONIC, SYMMETRY = 50.65e6, 6, 8
ENERGIES = [72, 150, 300, 450, 530]  # MeV, nearest ic.dat rows
STEPS = 1440
RFT_SAMPLES, RFT_TOL = 2880, 1e-12  # as G1 step 1
PAD = 5
FLOAT = re.compile(r'[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?')
TC1 = (
    'tc1: TRIMCOIL, TYPE="PSI-BFIELD-MIRRORED", RMIN=4350.0, '
    'RMAX=4470.0, BMAX=1.4e-3, SLPTC=0.6;'
)


def read_ring(path):
    """RING map as OPAL reads it (getFieldFromFile_Ring): r [m], theta
    step [deg], Bz and three theta derivatives.

    Arrays (nrad, ntet) in file units (kG, per file azimuth unit); the
    first data record is skipped as in OPAL.
    """
    text, pos = path.read_text(), 0

    def scan(width=None):
        # fscanf: skip whitespace, then %s (up to whitespace) or %16lE
        # (the longest float within 16 characters). The numbers are
        # fixed width, so a negative value touches its neighbour.
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
    return (
        np.roll(f, 2, 1)
        - 8 * np.roll(f, 1, 1)
        + 8 * np.roll(f, -1, 1)
        - np.roll(f, -2, 1)
    ) / (12 * h)


def ring_field(r, dt_deg, bz_kg, direction):
    """Median-plane field b(x, y, z) in T for protons travelling
    counterclockwise (LIONS frame, Bz < 0).

    direction +1: the map azimuth increases along the motion; -1: it
    decreases. Quintic spline, periodic per sector.
    """
    ntet = bz_kg.shape[1]
    jj = np.arange(-PAD, ntet + PAD)
    spline = RectBivariateSpline(
        r, np.radians(dt_deg) * jj, -0.1 * bz_kg[:, jj % ntet], kx=5, ky=5
    )
    sector = 2 * np.pi / SYMMETRY

    def field(x, y, z):
        t = (direction * np.arctan2(y, x)) % sector
        bz = spline.ev(np.hypot(x, y), t)
        return np.zeros_like(bz), np.zeros_like(bz), bz

    return field


def ic_rows():
    """Return the ic.dat row nearest each energy in ENERGIES."""
    rows = np.loadtxt(PSI / 'ic.dat')
    return {e_: rows[np.argmin(np.abs(rows[:, 0] - e_))] for e_ in ENERGIES}


def opal_closed_orbit(work, energy_mev, r_mm, pr_ic, trim):
    """OPAL 2024.1 closed-orbit finder (tune mode) at one energy: (r
    [m], pr [beta gamma], Qr, Qz).
    """
    work.mkdir(parents=True, exist_ok=True)
    (work / 'data').mkdir(exist_ok=True)
    shutil.copy(PSI / 'bfield.dat', work / 'bfield.dat')
    gev = energy_mev / 1000
    # FMLOWE, FMHIGHE, DENERGY in GeV in OPAL 2024.1; the element cannot
    # be named RING (protected keyword).
    (work / 'tune.in').write_text(f'''OPTION, ECHO=FALSE;
OPTION, VERSION=20240;
OPTION, CLOTUNEONLY=TRUE;
TITLE, STRING="PSI Ring closed orbit";
{TC1 if trim else ''}
psiring: CYCLOTRON, TYPE="RING", CYHARMON={HARMONIC},\
 PHIINIT=0.0,\
 PRINIT={pr_ic / 1000:.17g},\
 RINIT={r_mm:.17g},
    SYMMETRY={SYMMETRY}.0, RFFREQ={F_RF / 1e6:.17g},\
 FMAPFN="bfield.dat", FMLOWE={gev:.17g},\
 FMHIGHE={gev + 0.002:.17g}
    {', TRIMCOIL=tc1' if trim else ''};
l1: LINE=(psiring);
dist: DISTRIBUTION, TYPE=GAUSSMATCHED, LINE=l1, NSTEPS={STEPS}, DENERGY=0.002,\
 MAXSTEPSCO=100, NSECTORS={SYMMETRY},
    SECTOR=FALSE;
fs: FIELDSOLVER, FSTYPE=NONE, MX=8, MY=8, MT=8, PARFFTX=true,\
 PARFFTY=true, PARFFTT=false,
    BCFFTX=open, BCFFTY=open, BCFFTT=open;
beam1: BEAM, PARTICLE=PROTON,\
 PC={np.sqrt((gev + MP_MEV / 1000) ** 2 - (MP_MEV / 1000) ** 2):.17g},\
 NPART=1,
    BCURRENT=1.0E-6, CHARGE=1.0, BFREQ={F_RF / 1e6:.17g};
TRACK, LINE=l1, BEAM=beam1, MAXSTEPS=720, STEPSPERTURN=720;
RUN, METHOD="CYCLOTRON-T", BEAM=beam1, FIELDSOLVER=fs, DISTRIBUTION=dist;
ENDTRACK;
STOP;
''')
    subprocess.run(
        [
            'docker',
            'run',
            '--rm',
            '-v',
            f'{work.resolve()}:/data',
            orbits.OPAL_IMAGE,
            'bash',
            '-c',
            'opal tune.in > tune.log 2>&1',
        ],
        check=True,
    )
    rows = np.atleast_2d(np.loadtxt(work / 'data' / 'tunes.dat'))
    if not len(rows) or abs(rows[0, 0] - energy_mev) > 1e-6:
        raise RuntimeError(
            f'{work}: no OPAL closed orbit at {energy_mev} MeV (see tune.log)'
        )
    _, r, pr, _, qr, qz = rows[0]
    return r, pr, qr, qz


def step1(outdir):
    """Step 1: file conventions (1a units, 1b Boris orbits, 1c OPAL)."""
    r, dt_deg, (bz, dbt, _, _) = read_ring(PSI / 'bfield.dat')
    records = []

    # 1a: units of the file theta derivatives
    for unit, h in (('per_deg', dt_deg), ('per_rad', np.radians(dt_deg))):
        diff = np.abs(dbt - periodic_d5(bz, h))
        for stat, v in (('max', diff.max()), ('median', np.median(diff))):
            records.append(
                {
                    'run': '1a',
                    'energy_mev': None,
                    'case': unit,
                    'metric': f'dbz_dtheta_diff_{stat}',
                    'value': float(v),
                    'unit': 'kG per file unit',
                }
            )
        print(
            f'1a {unit}: max {diff.max():.3e}, '
            f'median {np.median(diff):.3e} '
            f'(|dBz/dtheta| max {np.abs(dbt).max():.3e})',
            file=sys.stderr,
        )

    # 1b: closed orbits with the numpy Boris tracker, both directions
    ic = ic_rows()
    for direction in (1, -1):
        field = ring_field(r, dt_deg, bz, direction)
        for energy, (e_ic, r_ic, pr_ic) in ic.items():
            gamma = 1 + e_ic / MP_MEV
            p = np.sqrt(gamma**2 - 1)
            readings = {
                'pr_1e-3_m0c': pr_ic / 1000,
                'pr_mrad': pr_ic / 1000 * p,
            }
            guess = (r_ic / 1000, direction * np.mean(list(readings.values())))
            tracker = orbits.BorisTracker(
                field, m_p, e, p, F_RF / HARMONIC, STEPS
            )
            case = f'direction{direction:+d}'
            try:
                s, res = orbits.closed_orbit(tracker.track, guess)
            except Exception as exc:
                records.append(
                    {
                        'run': '1b',
                        'energy_mev': e_ic,
                        'case': case,
                        'metric': 'failed',
                        'value': None,
                        'unit': f'{type(exc).__name__}: {exc}'[:200],
                    }
                )
                print(f'1b {case} {e_ic} MeV failed: {exc}', file=sys.stderr)
                continue
            values = {
                'r_minus_ic_mm': (s[0] * 1000 - r_ic, 'mm'),
                'residual': (res, 'm'),
                'pr_beta_gamma': (s[1], '1'),
                **{
                    f'pr_rel_error_{k}': (s[1] / v - 1, '1')
                    for k, v in readings.items()
                },
            }
            for k, (v, u) in values.items():
                records.append(
                    {
                        'run': '1b',
                        'energy_mev': e_ic,
                        'case': case,
                        'metric': k,
                        'value': float(v),
                        'unit': u,
                    }
                )
            print(
                f'1b {case} {e_ic} MeV '
                + ' '.join(f'{k}={v[0]:+.3e}' for k, v in values.items()),
                file=sys.stderr,
            )

    # 1c: OPAL closed-orbit finder, with and without TC1
    for energy in (72, 530):
        e_ic, r_ic, pr_ic = ic[energy]
        for trim in (False, True):
            case = 'tc1' if trim else 'no_trim'
            r_o, pr_o, qr, qz = opal_closed_orbit(
                outdir / f'opal_{energy}_{case}', e_ic, r_ic, pr_ic, trim
            )
            values = {
                'r_m': (r_o, 'm'),
                'pr_beta_gamma': (pr_o, '1'),
                'qr': (qr, '1'),
                'qz': (qz, '1'),
            }
            for k, (v, u) in values.items():
                records.append(
                    {
                        'run': '1c',
                        'energy_mev': e_ic,
                        'case': case,
                        'metric': k,
                        'value': float(v),
                        'unit': u,
                    }
                )
            print(
                f'1c {case} {e_ic} MeV '
                + ' '.join(f'{k}={v[0]:.7g}' for k, v in values.items()),
                file=sys.stderr,
            )

    meta = {
        'case': 'G2 step 1: PSI Ring file conventions',
        'proton_mass_mev': repr(MP_MEV),
        'rf_hz': repr(F_RF),
        'harmonic': str(HARMONIC),
        'steps_per_turn': str(STEPS),
        'opal_image': orbits.OPAL_IMAGE,
        **provenance(PSI / 'bfield.dat', PSI / 'ic.dat'),
        'machine': platform.platform(),
        'python': platform.python_version(),
    }
    outdir.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta),
        outdir / 'case_g2_step1.parquet',
    )


def g1_grids(r, dt_deg, bz, dbt):
    """The map in the case_g1 grids/attrs form (T, T/m, T/rad; sign of
    the file), dBz/dr by OPAL's 5-point scheme.
    """
    n = len(r)
    h = r[1] - r[0]
    dbr = np.empty_like(bz)
    for i in range(n):
        lo = min(
            max(i - 2, 0), n - 5
        )  # stencil clamped at the edges, as OPAL getdiffs
        x = (np.arange(lo, lo + 5) - i) * h
        w = np.linalg.solve(
            np.vander(x, increasing=True).T, [0, 1, 0, 0, 0]
        )  # d/dr at x = 0
        dbr[i] = w @ bz[lo : lo + 5]
    grids = {'bz': 0.1 * bz, 'dbz_dr': 0.1 * dbr, 'dbz_dtheta': 0.1 * dbt}
    attrs = {
        'sector_count': SYMMETRY,
        'ri_m': r[0],
        'rstep_m': h,
        'nrad': n,
        'theta_step_rad': np.radians(dt_deg),
    }
    return grids, attrs


def fold(q):
    """Fold a tune into [0, 0.5] by 1/(2 pi) arccos(cos(2 pi q))."""
    return np.arccos(np.cos(2 * np.pi * q)) / (2 * np.pi)


def setup():
    """Map in the G1 form, LIONS species and its mass, the quintic
    spline field, FIXPO rows (E, r, pr, Qr, Qz).
    """
    r, dt_deg, (bz, dbt, _, _) = read_ring(PSI / 'bfield.dat')
    grids, attrs = g1_grids(r, dt_deg, bz, dbt)
    species = case_g1.lions_species(MP_MEV * 1e6, 1)
    mass_kg = case_g1.lc.AMU * (species[0] - case_g1.lc.ELECTRON_MASS_U)
    ic, ref = np.loadtxt(PSI / 'ic.dat'), np.loadtxt(PSI / 'refsol.dat')
    if not np.array_equal(ic[:, 0], ref[:, 0]):
        raise ValueError('ic.dat and refsol.dat energies differ')
    return (
        grids,
        attrs,
        species,
        mass_kg,
        case_g1.smooth_field(grids, attrs, 1, 5),
        np.hstack([ic, ref[:, 1:]]),
    )


def reference(energy, r_ic, pr_ic, qr_ref, qz_ref):
    """p, the orbit_metrics reference dict and the folded FIXPO
    tunes."""
    p = np.sqrt((1 + energy / MP_MEV) ** 2 - 1)
    return (
        p,
        {
            'r_theta0_m': r_ic / 1000,
            'pr_over_p_theta0': pr_ic / 1000 / p,
            'f_rev_hz': F_RF / HARMONIC,
        },
        [fold(qr_ref), fold(qz_ref)],
    )


def step2(outdir):
    """Step 2: LIONS and the spline tracker on the map, no trim coil."""
    grids, attrs, species, mass_kg, spline, fixpo = setup()
    f_rev = F_RF / HARMONIC
    work = outdir / 'lions_map'
    work.mkdir(parents=True, exist_ok=True)
    with contextlib.chdir(work):
        ntet = case_g1.write_lions_map(grids, attrs, 1)
    records = []
    for energy, r_ic, pr_ic, qr_ref, qz_ref in fixpo:
        p, t, q_ref = reference(energy, r_ic, pr_ic, qr_ref, qz_ref)
        flagged = 386 <= energy <= 430
        for code, track in (
            (
                'lions',
                orbits.LionsTracker(work, energy, f_rev, STEPS, species).track,
            ),
            (
                'spline',
                orbits.BorisTracker(
                    spline, mass_kg, case_g1.lc.E_CHARGE, p, f_rev, STEPS
                ).track,
            ),
        ):
            base = {
                'code': code,
                'energy_mev': float(energy),
                'qr_half_integer_band': flagged,
            }
            try:
                values = case_g1.orbit_metrics(track, t, p, q_ref)
            except Exception as exc:
                records.append(
                    {
                        **base,
                        'metric': 'failed',
                        'value': None,
                        'unit': f'{type(exc).__name__}: {exc}'[:200],
                    }
                )
                print(f'{code} {energy} MeV failed: {exc}', file=sys.stderr)
                continue
            for k, (v, u) in values.items():
                records.append(
                    {**base, 'metric': k, 'value': float(v), 'unit': u}
                )
            print(
                f'{code} {energy} MeV '
                + ' '.join(f'{k}={v[0]:+.2e}' for k, v in values.items()),
                file=sys.stderr,
            )
    meta = {
        'case': 'G2 step 2: LIONS and the spline tracker on the PSI '
        'Ring map, no trim coil',
        'proton_mass_mev': repr(MP_MEV),
        'species': repr(species),
        'rf_hz': repr(F_RF),
        'harmonic': str(HARMONIC),
        'steps_per_turn': str(STEPS),
        'lions_map_columns': str(ntet),
        'dbz_dr': '5-point differences in r, clamped stencils (OPAL getdiffs)',
        'dbz_dtheta': 'file, per rad',
        'tunes': 'folded into [0, 0.5], refsol.dat folded the same way',
        **provenance(PSI / 'bfield.dat', PSI / 'ic.dat', PSI / 'refsol.dat'),
        'machine': platform.platform(),
        'python': platform.python_version(),
    }
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta),
        outdir / 'case_g2_step2.parquet',
    )


def block_determinants(outdir, name, case, make_trackers):
    """det - 1 of the r block, z block and M, and the tunes, at 152,
    292, 452 MeV for each tracker.

    make_trackers(energy, p, grids, attrs, species, mass_kg, spline,
    work) -> {code: (track, jacobian steps)}.
    """
    grids, attrs, species, mass_kg, spline, fixpo = setup()
    work = outdir / 'lions_map'
    work.mkdir(parents=True, exist_ok=True)
    with contextlib.chdir(work):
        case_g1.write_lions_map(grids, attrs, 1)
    records = []
    for energy, r_ic, pr_ic, _, _ in fixpo[[40, 110, 190]]:
        p, t, _ = reference(energy, r_ic, pr_ic, 0, 0)
        for code, (track, steps) in make_trackers(
            energy, p, grids, attrs, species, mass_kg, spline, work
        ).items():
            s, _ = orbits.closed_orbit(
                track, (t['r_theta0_m'], p * t['pr_over_p_theta0'])
            )
            for h in steps:
                m = orbits.jacobian(track, s, h, h * p / s[0])
                qr, qz = orbits.linear_tunes(m)
                values = {
                    'det_r_block': np.linalg.det(m[:2, :2]) - 1,
                    'det_z_block': np.linalg.det(m[2:, 2:]) - 1,
                    'det_m': np.linalg.det(m) - 1,
                    'qr': qr,
                    'qz': qz,
                }
                for k, v in values.items():
                    records.append(
                        {
                            'code': code,
                            'energy_mev': float(energy),
                            'jacobian_step_m': h,
                            'metric': k,
                            'value': float(v),
                            'unit': '1',
                        }
                    )
                print(
                    f'{code} {energy} MeV h={h:.0e} '
                    + ' '.join(f'{k}={v:+.3e}' for k, v in values.items()),
                    file=sys.stderr,
                )
    meta = {
        'case': case,
        'steps_per_turn': str(STEPS),
        'rftrack_samples': str(RFT_SAMPLES),
        'rftrack_tolerance': str(RFT_TOL),
        **provenance(PSI / 'bfield.dat', PSI / 'ic.dat'),
        'machine': platform.platform(),
        'python': platform.python_version(),
    }
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta),
        outdir / f'case_g2_{name}.parquet',
    )


def step2b(outdir):
    """LIONS (three Jacobian steps), numpy Boris on the LIONS bilinear
    field, spline tracker.
    """

    def trackers(energy, p, grids, attrs, species, mass_kg, spline, work):
        f_rev, q = F_RF / HARMONIC, case_g1.lc.E_CHARGE
        bilinear = case_g1.bilinear_field(grids, attrs, 1)
        return {
            'lions': (
                orbits.LionsTracker(work, energy, f_rev, STEPS, species).track,
                (1e-7, 1e-6, 1e-5),
            ),
            'bilinear_numpy': (
                orbits.BorisTracker(
                    bilinear, mass_kg, q, p, f_rev, STEPS
                ).track,
                (case_g1.H,),
            ),
            'spline': (
                orbits.BorisTracker(spline, mass_kg, q, p, f_rev, STEPS).track,
                (case_g1.H,),
            ),
        }

    block_determinants(
        outdir,
        'step2b',
        'G2 step 2b: block determinants of the one-turn map',
        trackers,
    )


def step2c(outdir):
    """RF-Track (adaptive) on the LIONS bilinear field."""

    def trackers(energy, p, grids, attrs, species, mass_kg, spline, work):
        r_max = attrs['ri_m'] + attrs['rstep_m'] * (attrs['nrad'] - 1)
        rt = orbits.RFTrackTracker(
            case_g1.bilinear_field(grids, attrs, 1),
            r_max,
            mass_kg,
            p,
            HARMONIC / F_RF,
            RFT_SAMPLES,
            RFT_TOL,
        )
        return {'rftrack_bilinear': (rt.track, (case_g1.H,))}

    block_determinants(
        outdir,
        'step2c',
        'G2 step 2c: RF-Track on the LIONS bilinear field',
        trackers,
    )


def step3(outdir):
    """OPAL closed-orbit finder (RING reader, no trim coil) and RF-Track
    on the spline, every 20 MeV and 530 MeV.
    """
    _, attrs, _, mass_kg, spline, fixpo = setup()
    r_max = attrs['ri_m'] + attrs['rstep_m'] * (attrs['nrad'] - 1)
    records = []
    for energy, r_ic, pr_ic, qr_ref, qz_ref in [*fixpo[::10], fixpo[-1]]:
        p, t, q_ref = reference(energy, r_ic, pr_ic, qr_ref, qz_ref)

        def opal():
            r_o, pr_o, qr, qz = opal_closed_orbit(
                outdir / f'opal_{energy:.3f}', energy, r_ic, pr_ic, False
            )
            return {
                'r_eo_rel_error': (r_o / t['r_theta0_m'] - 1, '1'),
                'pr_eo_error_over_p': (
                    (pr_o - p * t['pr_over_p_theta0']) / p,
                    '1',
                ),
                'qr_error': (fold(qr) - q_ref[0], '1'),
                'qz_error': (fold(qz) - q_ref[1], '1'),
                'qr': (qr, '1'),
                'qz': (qz, '1'),
            }

        def rftrack():
            tracker = orbits.RFTrackTracker(
                spline,
                r_max,
                mass_kg,
                p,
                HARMONIC / F_RF,
                RFT_SAMPLES,
                RFT_TOL,
            )
            return case_g1.orbit_metrics(tracker.track, t, p, q_ref)

        for code, run in (('opal', opal), ('rftrack', rftrack)):
            base = {
                'code': code,
                'energy_mev': float(energy),
                'qr_half_integer_band': bool(386 <= energy <= 430),
            }
            try:
                values = run()
            except Exception as exc:
                records.append(
                    {
                        **base,
                        'metric': 'failed',
                        'value': None,
                        'unit': f'{type(exc).__name__}: {exc}'[:200],
                    }
                )
                print(f'{code} {energy} MeV failed: {exc}', file=sys.stderr)
                continue
            for k, (v, u) in values.items():
                records.append(
                    {**base, 'metric': k, 'value': float(v), 'unit': u}
                )
            print(
                f'{code} {energy} MeV '
                + ' '.join(f'{k}={v[0]:+.2e}' for k, v in values.items()),
                file=sys.stderr,
            )
    meta = {
        'case': 'G2 step 3: OPAL (RING reader, closed-orbit finder) and '
        'RF-Track (spline) on the PSI Ring map',
        'proton_mass_mev': repr(MP_MEV),
        'rf_hz': repr(F_RF),
        'harmonic': str(HARMONIC),
        'opal_image': orbits.OPAL_IMAGE,
        'opal_steps_per_turn': str(STEPS),
        'rftrack_samples': str(RFT_SAMPLES),
        'rftrack_tolerance': str(RFT_TOL),
        'trim_coil': 'none',
        **provenance(PSI / 'bfield.dat', PSI / 'ic.dat', PSI / 'refsol.dat'),
        'machine': platform.platform(),
        'python': platform.python_version(),
    }
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta),
        outdir / 'case_g2_step3.parquet',
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--step', choices=('1', '2', '2b', '2c', '3'), default='1'
    )
    parser.add_argument(
        '--outroot',
        type=Path,
        default=ROOT / 'benchmarks' / 'runs' / 'case_g2',
    )
    args = parser.parse_args()
    {'1': step1, '2': step2, '2b': step2b, '2c': step2c, '3': step3}[
        args.step
    ](args.outroot.resolve() / f'step{args.step}')
