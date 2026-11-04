"""Case V: LIONS judges the tunes of a mission,
docs/training/PHYSICAL.md (F18) and PHYSICAL_NEON.md (NEON). Synthetic
map: says nothing about a real machine.

usage: CASE_P_MISSION=NEON uv run --extra game --extra rl --extra bo
--group bench --group dev python benchmarks/case_v.py
{tunes,check,verdict,limit} [--outroot DIR] [--mapdir DIR] [--reference
DIR] [--workers N] [--start-mev W] [--extraction {mean,energy}] tunes:
the tunes of case T3 that read in spec, regenerated and compared with
its table. The F18 tunes of case T need the code at a39faa0: the score
has changed since. check: time one pass, sign of the phase conversion,
field scaling. verdict: every tune through LIONS, one pass per tune;
finished tunes are kept in parts/ and skipped on a rerun. limit: every
tune to 4 times T0's turn limit, read at 1, 2 and 4 times it, with
particles still circulating reported apart. --start-mev: start at the
closed-orbit energy nearest above W, for the verdict without the central
region; case_v_tunes.parquet is then copied into --outroot.
--phase-rms-deg: Gaussian start phases of this rms around the tune's
phase, in place of T0's buncher law. Deterministic quantiles.
--extraction: mean, the orbit of mean radius R, as F18. energy: the
orbit at T0's energy for R, where beta c / omega = R.
"""

import argparse
import contextlib
import io
import json
import multiprocessing
import os
import pickle
import platform
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))

import case_p  # noqa: E402 (after sys.path)
import lions  # noqa: E402 (after sys.path)
import lions.cli  # noqa: E402 (after sys.path)
import orbits  # noqa: E402 (after sys.path)
import synthetic as syn  # noqa: E402 (after sys.path)
from lions.game import machine as mc  # noqa: E402 (after sys.path)
from lions.game.reduced import gaussian_phases, install_from_env  # noqa: E402 (after sys.path)
from lions.game import tuners as t  # noqa: E402 (after sys.path)
from lions.game.env import mission_by_id  # noqa: E402 (after sys.path)
from provenance import provenance  # noqa: E402 (after sys.path)

MISSION = case_p.MISSION
MODEL = install_from_env()
TUNERS = {
    'coach': t.coach,
    'turbo': t.turbo,
    'constrained_bo': t.constrained_bo,
    'mobo': t.mobo,
}
N = mc.MACROPARTICLES
PPO_SEED = 0
MAP_FILES = ('bmap.dat', 'pot1.bin', 'pot2.dat')
REST_MEV = case_p.MASS_KG * syn.C**2 / syn.MEV
DEFLECTOR_DEG, DEFLECTOR_R = (
    250.0,
    1.513,
    # entry of the Mo deck: another ion on another map, a labelled
    # reading, decides nothing
)
EXTRACTED, Z_APERTURE, OUT_OF_MAP, TURN_LIMIT = (
    'extracted',
    'z_aperture',
    'out_of_map',
    'turn_limit',
)
LIMIT_FACTORS = (1, 2, 4)
BEAM_FIELDS = (
    'transmission',
    'extracted_ua',
    'w_mean',
    'w_spread',
    'lost_kw',
    'turns',
    'x',
    'y',
)


def wrap(phi):
    """Wrap phi [rad] into [-pi, pi)."""
    return np.remainder(phi + np.pi, 2 * np.pi) - np.pi


def lions_phi0(phi, crest_deg):
    """T0 phase [rad, 0 at the crest] to the RF phase omega t [deg] of a
    LIONS particle crossing theta = 0. Later is larger in both.
    """
    return crest_deg + np.degrees(phi)


def t0_phase(omega_t, crest_deg):
    """Inverse of lions_phi0, in [-pi, pi). omega_t [rad]."""
    return wrap(omega_t - np.radians(crest_deg))


class Map:
    """What the verdict takes from case P: start orbit, crest,
    extraction radius.
    """

    def __init__(self, mapdir, start_mev=None, extraction='mean'):
        self.dir, self.extraction = Path(mapdir).resolve(), extraction
        rec = [
            r
            for r in pq.read_table(
                self.dir / 'case_p_orbits.parquet'
            ).to_pylist()
            if r['code'] == 'boris_analytic'
        ]
        gains = pq.read_table(self.dir / 'case_p_gain.parquet').to_pylist()
        by_energy = [(r['w_mev'], r['r_co_m']) for r in rec]
        table = dict(by_energy)
        self.w0 = (
            case_p.start_energy(gains, table)
            if start_mev is None
            else min(a[0] for a in by_energy if a[0] >= start_mev)
        )
        self.r0 = table[self.w0]
        self.crest_deg = float(
            np.median(
                [g['crest_phi0_deg'] for g in gains if g['w_mev'] >= self.w0]
            )
        )
        self.r_ext = float(
            np.interp(
                case_p.MACHINE.radius,
                [r['rbar_m'] for r in rec],
                [r['r_co_m'] for r in rec],
            )
            if extraction == 'mean'
            else np.interp(
                case_p.momentum(case_p.W_EXT),
                case_p.momentum(np.array([r['w_mev'] for r in rec])),
                [r['r_co_m'] for r in rec],
            )
        )
        self.inputs = [
            self.dir / name
            for name in (
                'bmap.dat',
                'pot2.dat',
                'case_p_gain.parquet',
                'case_p_orbits.parquet',
            )
        ]

    def workdir(self, path):
        """Make path hold symlinks to the map files; return it."""
        path.mkdir(parents=True, exist_ok=True)
        for name in MAP_FILES:
            if not (path / name).exists():
                os.symlink(self.dir / name, path / name)
        return path


def load(deck_args):
    """Write deck.dat from deck_args and case P settings; load it."""
    syn.write_deck(
        'deck.dat',
        npaspt=case_p.STEPS,
        teta1=case_p.MACHINE.dee_angle_deg / 2,
        species=case_p.SPECIES,
        tail=orbits.deck_tail(case_p.OPTIONS),
        **deck_args,
    )
    d = lions.State()
    d.run.swap_bytes = 'n'
    lions.load(d, 'deck.dat')
    return d


def reference_orbit(workdir, mp):
    """Closed orbit at the start energy, no RF: x, y, z, px, py, pz at
    each step over one turn, in LIONS' staggered convention.
    """
    rows = []

    def record(d):
        b = d.bunch
        rows.append([b.x[0], b.y[0], b.z[0], b.px[0], b.py[0], b.pz[0]])

    with contextlib.chdir(workdir), contextlib.redirect_stdout(io.StringIO()):
        d = load(
            dict(
                fhf=case_p.HARMONIC * case_p.F0,
                r0=mp.r0,
                w0=mp.w0,
                turns=1,
                v_dee=0.0,
                harm=float(case_p.HARMONIC),
            )
        )
        lions.cli.check_box(d.bunch, d.field, d.box)
        lions.cli.compute_fields(d)
        orbits.half_kick(d, -0.5)  # the deck momentum is p(-1/2)
        lions.cli.track(d, on_step=record)
    return np.array(rows)


class Pass:
    """One pass of 120 macroparticles through LIONS. The RF is one clock
    for the whole bunch, so a particle's RF phase is its arrival time at
    theta = 0.

    Each one starts before theta = 0 on the closed orbit, with its
    electric kick off until it arrives there: it arrives at the start
    energy, at its own phase.
    """

    def __init__(
        self,
        workdir,
        mp,
        ref,
        knobs,
        errors,
        draws,
        max_turns=None,
        scale_in_map=False,
        factors=(),
        phase_rms_deg=None,
    ):
        machine, ion = case_p.MACHINE, case_p.ION
        self.mp, self.knobs, self.errors = mp, knobs, errors
        self.b0 = knobs['b0'] * (1 + errors.field)
        self.s = self.b0 / mc.nominal(machine, ion, case_p.W_EXT * case_p.A)[0]
        w_inj = machine.w_injection * ion.mass_number
        self.g0 = ion.gamma(w_inj)
        f_rf = knobs['rf_mhz'] * 1e6
        self.h = max(round(f_rf / ion.f_rev(w_inj, self.b0)), 1)
        if self.h != case_p.HARMONIC:
            raise ValueError(
                f'harmonic {self.h}: the map is built for {case_p.HARMONIC}'
            )
        self.limit = (
            int(np.ceil(mc.turn_limit(machine, ion, knobs, self.b0, self.h)))
            if max_turns is None
            else max_turns
        )
        self.base, self.snaps = self.limit, {}
        self.read_at = {self.base * f * case_p.STEPS: f for f in factors}
        if factors:
            self.limit = self.base * max(factors)
        self.phi = (
            mc.start_phases(knobs, errors)
            if phase_rms_deg is None
            else gaussian_phases(knobs, errors, phase_rms_deg)
        )
        self.z0 = mc.start_heights(machine, draws) * 1e-3
        hold = (
            (self.phi + np.pi) / (2 * np.pi) * case_p.STEPS / self.h
        )  # steps to theta = 0, under half a turn at harmonic 2
        at = case_p.STEPS - hold
        start = np.array(
            [
                np.interp(at, np.arange(case_p.STEPS + 1), ref[:, c])
                for c in range(6)
            ]
        )
        self.release = np.ceil(hold).astype(int)
        self.fate = np.full(N, TURN_LIMIT, dtype=object)
        self.w_end, self.w_last, self.w_prev = (
            np.full(N, mp.w0 * case_p.A) for _ in range(3)
        )
        self.r_deflector = np.full(N, np.nan)
        self.turns = np.zeros(N, dtype=int)
        self.phi_first = np.full(N, np.nan)
        self.flagged, self.w_flag = np.zeros(N, dtype=bool), np.full(N, np.nan)
        self.phase_sum, self.counts = (
            np.zeros((N, mc.PROBE_BINS)),
            np.zeros((N, mc.PROBE_BINS)),
        )
        phi0 = lions_phi0(
            -np.pi, mp.crest_deg
        )  # the clock starts when hold = 0 has phase -pi
        with (
            contextlib.chdir(workdir),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            d = load(
                dict(
                    fhf=f_rf,
                    r0=mp.r0,
                    w0=mp.w0,
                    turns=self.limit + 1,
                    v_dee=knobs['dee_kv'] * 1e3,
                    harm=float(self.h),
                    nbrion=N,
                    phi0=phi0,
                )
            )
            b = d.bunch
            d.traps.follow = 0  # dead particles stop costing steps
            b.x[:], b.y[:], b.z[:] = (
                start[0] / self.s,
                start[1] / self.s,
                self.z0,
            )
            b.px[:], b.py[:], b.pz[:] = start[3], start[4], 0.0
            b.gamma[:] = np.sqrt(1 + b.px**2 + b.py**2 + b.pz**2)
            if not scale_in_map:
                # scales Bx, By, Bz together; b_factor in the deck is
                # only printed
                b.boris_b_coef *= self.s
            self.e_coef = b.boris_e_coef.copy()
            b.boris_e_coef[self.release > 0] = 0.0
            self.prev = b.x.copy(), b.y.copy(), b.gamma.copy(), d.run.t
            self.rmin, self.rmax = d.field.rmin, d.field.rmax - 0.02
            lions.cli.track(d, on_step=self.step)
        self.turn_max = int(self.turns.max())
        for f in factors:
            self.snaps.setdefault(f, self.state(np.hypot(b.x, b.y)))

    def step(self, d):
        """Per-step hook: release, crossings, fates. True to go on."""
        k, b = d.run.step, d.bunch
        if k == 0:
            return True
        released = self.release == k
        b.boris_e_coef[released] = self.e_coef[released]
        xp, yp, gp, tp = self.prev
        live = b.alive.copy()
        c = np.flatnonzero(live & (yp < 0) & (b.y >= 0) & (b.x > 0))
        if len(c):
            self.crossing(c, d, xp, yp, gp, tp)
        w = (b.gamma - 1) * REST_MEV
        self.w_end[live] = w[live]
        r = np.hypot(b.x, b.y)
        line = np.radians(DEFLECTOR_DEG)

        def side(x, y):
            return y * np.cos(line) - x * np.sin(line)

        at = (
            live
            & (side(xp, yp) < 0)
            & (side(b.x, b.y) >= 0)
            & (b.x * np.cos(line) + b.y * np.sin(line) > 0)
        )
        self.r_deflector[at] = r[at]
        over_z = live & (np.abs(b.z) > case_p.MACHINE.z_aperture * 1e-3)
        out = live & ~over_z & ((r < self.rmin) | (r >= self.rmax))
        self.fate[over_z], self.fate[out] = Z_APERTURE, OUT_OF_MAP
        b.alive[over_z | out] = False
        self.prev = b.x.copy(), b.y.copy(), b.gamma.copy(), d.run.t
        if k in self.read_at:
            self.snaps[self.read_at[k]] = self.state(r)
        return bool(b.alive.any()) and k < self.limit * case_p.STEPS

    def state(self, r):
        """Copy of the per-particle arrays at radius r [m]."""
        return {
            'fate': self.fate.copy(),
            'w_end': self.w_end.copy(),
            'w_last': self.w_last.copy(),
            'w_prev': self.w_prev.copy(),
            'turns': self.turns.copy(),
            'r': r.copy(),
            'counts': self.counts.sum(0),
            'phase_sum': self.phase_sum.sum(0),
        }

    def crossing(self, c, d, xp, yp, gp, tp):
        """Update particles c that crossed y = 0, interpolated."""
        b, ion = d.bunch, case_p.ION
        f = -yp[c] / (b.y[c] - yp[c])
        when = tp + f * (d.run.t - tp)
        w = (gp[c] + f * (b.gamma[c] - gp[c]) - 1) * REST_MEV
        r = xp[c] + f * (b.x[c] - xp[c])
        phi = t0_phase(d.rf.omega * when, self.mp.crest_deg)
        first = self.turns[c] == 0
        self.phi_first[c[first]] = phi[first]
        self.turns[c] += 1
        self.w_prev[c] = self.w_last[c]
        self.w_last[c] = w
        g = ion.gamma(w)
        rho = ion.rigidity(w) / (
            self.b0 * (1 + self.knobs['trim'] * (g / self.g0 - 1))
        )  # as machine.run: the radial probe of both models
        i = np.minimum(
            (rho / case_p.MACHINE.radius * mc.PROBE_BINS).astype(int),
            mc.PROBE_BINS - 1,
        )
        self.phase_sum[c, i] += np.degrees(phi)
        self.counts[c, i] += 1
        behind = (np.cos(phi) < mc.BEHIND_RF) & ~self.flagged[c]
        self.flagged[c[behind]] = True
        self.w_flag[c[behind]] = w[behind]
        done = c[r >= self.mp.r_ext]
        self.fate[done] = EXTRACTED
        b.alive[done] = False

    def beam(self, behind_rule=False):
        """Beam through the closed forms of machine.extract.
        behind_rule: T0's loss rule applied on top of the dynamics.
        """
        out = self.fate == EXTRACTED
        if behind_rule:
            ok = out & ~self.flagged
            power = (
                self.w_flag[self.flagged].sum()
                + self.w_end[~out & ~self.flagged].sum()
            )
        else:
            ok, power = out, self.w_end[~out].sum()
        return mc.extract(
            case_p.MACHINE,
            case_p.ION,
            self.knobs,
            self.errors,
            self.w_last,
            self.w_prev,
            ok,
            power,
            self.turn_max,
            self.h,
            self.counts.sum(0),
            self.phase_sum.sum(0),
            self.counts.sum(0),
        )

    def beam_at(self, s, count_circulating):
        """Beam read at a snapshot. A particle still circulating has no
        fate: count_circulating says whether its energy enters the lost
        power.
        """
        out = s['fate'] == EXTRACTED
        lost = ~out if count_circulating else ~out & (s['fate'] != TURN_LIMIT)
        return mc.extract(
            case_p.MACHINE,
            case_p.ION,
            self.knobs,
            self.errors,
            s['w_last'],
            s['w_prev'],
            out,
            s['w_end'][lost].sum(),
            int(s['turns'].max()),
            self.h,
            s['counts'],
            s['phase_sum'],
            s['counts'],
        )


def judged(condition, knobs, field, beam=None):
    """A mission in the state of the tune: spec terms and verdict for a
    beam, T0's own if none.
    """
    m = mission_by_id(MISSION)
    m.reset(condition)
    m.knobs, m.errors.field = dict(knobs), field
    if beam is None:
        m.refresh()
    else:
        m.beam = beam
        m.after_beam()
    return m


def row_of(prefix, m):
    """Beam fields and spec terms of mission m, keys prefixed."""
    b = m.beam
    return {
        **{f'{prefix}_{k}': float(getattr(b, k)) for k in BEAM_FIELDS},
        f'{prefix}_terms': json.dumps(m.spec_terms()),
        f'{prefix}_in_spec': bool(m.ready()),
    }


def write(records, path, what, mp, units, phase_rms_deg=None):
    """Write records to Parquet at path, with map metadata."""
    meta = {
        'case': f'synthetic {MISSION}, {what}',
        'registration': case_p.REGISTRATION,
        'label': 'SYNTHETIC',
        **provenance(*mp.inputs),
        'machine': platform.platform(),
        'python': platform.python_version(),
        'steps_per_turn': str(case_p.STEPS),
        'interpolation': 'hermite',
        'start_mev': str(mp.w0),
        'crest_deg': str(mp.crest_deg),
        'r_ext_m': str(mp.r_ext),
        'tuning_model': MODEL or 'T0',
        'extraction': mp.extraction,
        'start_phase_law': 'buncher'
        if phase_rms_deg is None
        else f'gaussian, {phase_rms_deg} deg rms, buncher unused',
        'units': units,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist(records).replace_schema_metadata(meta), path
    )


def tune_job(job):
    """One tune of case T3, regenerated: the Parquet holds no field
    error. torch on one thread, as there.
    """
    import torch

    torch.set_num_threads(1)
    was = job
    m = mission_by_id(MISSION)
    m.reset(was['condition'])
    run = t.tune(m, TUNERS[was['method']], was['condition'])
    knobs = {k: float(m.knobs[k]) for k in run.keys}
    same = (
        run.reached is not None
        and float(run.reached) == was['evals_to_spec']
        and knobs == json.loads(was['final_knobs'])
        and bool(m.ready()) == was['true_at_stop']
    )
    print(
        was['method'],
        was['condition'],
        run.reached,
        was['evals_to_spec'],
        '' if same else 'DIFFERENT',
        flush=True,
    )
    return same, {
        'method': was['method'],
        'training_seed': -1,
        'condition': was['condition'],
        'evals': run.reached,
        'field_error': float(m.errors.field),
        'knobs': json.dumps(m.knobs),
        't0_true_at_stop': bool(m.ready()),
        'harmonic': harmonic_of(m.knobs, m.errors),
    }


def harmonic_of(knobs, errors):
    """RF harmonic number at injection for the knobs and field error."""
    w_inj = case_p.MACHINE.w_injection * case_p.ION.mass_number
    return max(
        round(
            knobs['rf_mhz']
            * 1e6
            / case_p.ION.f_rev(w_inj, knobs['b0'] * (1 + errors.field))
        ),
        1,
    )


def tunes(outdir, reference, mp, workers):
    """Every tune that read in spec on the original box in case T3. Each
    must end as in its table, or the run stops.
    """
    runs = [
        r
        for r in pq.read_table(
            reference / 'case_t3_classical.parquet', use_threads=False
        ).to_pylist()
        if r['mission'] == MISSION
        and r['box'] == 'original'
        and r['in_spec_measured']
    ]
    with multiprocessing.get_context('spawn').Pool(workers) as pool:
        done = pool.map(tune_job, runs, chunksize=1)
    rows = [{'tune': i, **d[1]} for i, d in enumerate(done)]
    write(
        rows,
        outdir / 'case_v_tunes.parquet',
        'tunes',
        mp,
        'field_error 1, knobs as the scenario',
    )
    wrong = sum(not d[0] for d in done)
    print(
        f'{len(rows)} tunes; true at stop '
        f'{sum(r["t0_true_at_stop"] for r in rows)}; harmonics '
        f'{sorted({r["harmonic"] for r in rows})}; '
        f'{wrong} differ from case T3'
    )
    return wrong


def tunes_from_table(outdir, reference, mp):
    """Tunes read in spec on the tuning model, from their table: PPO
    seed 0 and the other methods. The stored state must give the stored
    true read, or the run stops.
    """
    runs = [
        r
        for sub in ('classical', 'rl')
        for r in pq.read_table(
            reference / f'case_t3_{sub}.parquet', use_threads=False
        ).to_pylist()
        if r['mission'] == MISSION
        and r['in_spec_measured']
        and r['training_seed'] in (-1, PPO_SEED)
    ]
    rows = [
        {
            'tune': i,
            'method': r['method'],
            'training_seed': r['training_seed'],
            'condition': r['condition'],
            'evals': r['evals_to_spec'],
            'field_error': r['field_error_at_stop'],
            'knobs': r['knobs_at_stop'],
            't0_true_at_stop': r['true_at_stop'],
        }
        for i, r in enumerate(runs)
    ]
    for r in rows:
        m = of_tune(r)
        r['harmonic'] = harmonic_of(m.knobs, m.errors)
    wrong = sum(bool(of_tune(r).ready()) != r['t0_true_at_stop'] for r in rows)
    write(
        rows,
        outdir / 'case_v_tunes.parquet',
        'tunes',
        mp,
        'field_error 1, knobs as the scenario',
    )
    print(
        f'{len(rows)} tunes; true at stop '
        f'{sum(r["t0_true_at_stop"] for r in rows)}; harmonics '
        f'{sorted({r["harmonic"] for r in rows})}; '
        f'{wrong} differ from the table'
    )
    return wrong


def of_tune(row):
    """Mission in the state of a tune row."""
    return judged(
        row['condition'], json.loads(row['knobs']), row['field_error']
    )


def check(outdir, mp):
    """Time one pass and check the phase and field-scale conversions."""
    work = mp.workdir(outdir / 'check' / 'base')
    ref = reference_orbit(work, mp)
    print(
        'map: start',
        mp.w0,
        'MeV, r0',
        mp.r0,
        'm; crest',
        mp.crest_deg,
        'deg; r_ext',
        mp.r_ext,
        'm',
    )
    print(
        'closed orbit closes: |dr| after one turn',
        float(np.hypot(*(ref[-1, :2] - ref[0, :2]))),
        'm',
    )
    row = pq.read_table(outdir / 'case_v_tunes.parquet').to_pylist()[0]
    m = of_tune(row)
    began = time.perf_counter()
    p = Pass(work, mp, ref, m.knobs, m.errors, m.draws)
    print(
        f'pass: {time.perf_counter() - began:.1f} s wall, '
        f'limit {p.limit} turns, longest {p.turn_max}, fates',
        dict(zip(*np.unique(p.fate.astype(str), return_counts=True))),
    )
    print(
        'first crossing against the set phase: max |diff|',
        float(np.nanmax(np.abs(np.degrees(wrap(p.phi_first - p.phi))))),
        'deg; never crossed',
        int(np.isnan(p.phi_first).sum()),
    )
    # Sign and size of the phase conversion against T0's update: phi +=
    # 2 pi (f_rf / f_rev - h) per turn.
    calm = {**m.knobs, 'phase_deg': -m.errors.phase_deg, 'buncher': 0.0}
    turns, nominal_b0 = (
        60,
        mc.nominal(case_p.MACHINE, case_p.ION, case_p.W_EXT * case_p.A)[0],
    )
    calm['b0'] = nominal_b0 / (1 + m.errors.field)
    calm['rf_mhz'] = case_p.HARMONIC * case_p.F0 / 1e6
    for name, change in (
        ('rf +1e-4', {'rf_mhz': calm['rf_mhz'] * (1 + 1e-4)}),
        ('b0 +1e-4', {'b0': calm['b0'] * (1 + 1e-4)}),
    ):
        slip = []
        for knobs in (calm, {**calm, **change}):
            q = Pass(work, mp, ref, knobs, m.errors, m.draws, max_turns=turns)
            on = (np.abs(q.phi) < np.radians(20)) & (q.turns > 0)
            last = (
                q.phase_sum.sum(1)[on] / q.counts.sum(1)[on]
            ).mean()  # mean over the turns of each particle near the crest
            slip.append(last)
        predicted = (
            360
            * case_p.HARMONIC
            * (1e-4 if 'rf' in name else -1e-4)
            * turns
            / 2
        )
        print(
            f'{name}: mean phase over {turns} turns moves '
            f'{slip[1] - slip[0]:+.3f} deg, '
            f'T0 update predicts {predicted:+.3f} deg'
        )
    # boris_b_coef against a map rewritten with the same scale.
    s = 1.002
    knobs = {**calm, 'b0': nominal_b0 * s / (1 + m.errors.field)}
    field = case_p.load_field(mp.dir)
    scaled = mp.workdir(outdir / 'check' / 'scaled')
    (scaled / 'bmap.dat').unlink()

    def nodes(i):
        return lambda r, th: s * field.nodes(r, th)[i]

    with contextlib.chdir(scaled):
        syn.write_bmap(
            'bmap.dat',
            nodes(0),
            nodes(1),
            nodes(2),
            nodes(3),
            rmin=0.005,
            rmax=case_p.MACHINE.radius * 1.15,
            nr=case_p.MAP_NR,
            ntet=361,
        )
    a = Pass(work, mp, ref, knobs, m.errors, m.draws, max_turns=20)
    c = Pass(
        scaled,
        mp,
        ref,
        knobs,
        m.errors,
        m.draws,
        max_turns=20,
        scale_in_map=True,
    )
    print(
        'scale 1.002 by coefficient against rewritten map, 20 turns: max |dw|',
        float(np.max(np.abs(a.w_end - c.w_end))),
        'MeV; of',
        float(a.w_end.max()),
        'MeV',
    )
    shutil.rmtree(outdir / 'check')


def one(job):
    """Run one tune through LIONS; keep its rows in parts/."""
    row, mapdir, parts, start_mev, extraction, phase_rms_deg = job
    path = parts / f"{row['tune']:03d}.pkl"
    if path.exists():
        return row['tune']
    mp = Map(mapdir, start_mev, extraction)
    work = mp.workdir(parts / f"work_{row['tune']:03d}")
    ref = reference_orbit(work, mp)
    m = of_tune(row)
    t0 = row_of('t0', m)
    t0_beam = m.beam
    began = time.perf_counter()
    p = Pass(
        work, mp, ref, m.knobs, m.errors, m.draws, phase_rms_deg=phase_rms_deg
    )
    wall = time.perf_counter() - began
    base = {
        k: row[k] for k in ('tune', 'method', 'training_seed', 'condition')
    }
    inputs = {
        'rf_mhz': float(m.knobs['rf_mhz']),
        'b0_t': float(p.b0),
        'dee_kv': float(m.knobs['dee_kv']),
        'phase_deg': float(m.knobs['phase_deg'] + m.errors.phase_deg),
        'buncher': float(m.knobs['buncher']),
    }
    lions_m = judged(row['condition'], m.knobs, row['field_error'], p.beam())
    rule_m = judged(
        row['condition'], m.knobs, row['field_error'], p.beam(behind_rule=True)
    )
    reading = p.fate == EXTRACTED
    verdict = {
        **base,
        **inputs,
        'field_scale': float(p.s),
        'wall_s': wall,
        'limit_turns': p.limit,
        'extracted_past_deflector_entry': float(
            np.mean(p.r_deflector[reading] >= DEFLECTOR_R)
        )
        if reading.any()
        else float('nan'),
        **t0,
        **row_of('lions', lions_m),
        **row_of('rule', rule_m),
    }
    with np.errstate(invalid='ignore', divide='ignore'):
        by_bin = p.phase_sum / p.counts
    particles = [
        {
            **base,
            **inputs,
            'particle': i,
            'phi_start_deg': float(np.degrees(p.phi[i])),
            'phi_first_deg': float(np.degrees(p.phi_first[i])),
            'z0_m': float(p.z0[i]),
            'survives': p.fate[i] == EXTRACTED,
            'fate': str(p.fate[i]),
            'turns': int(p.turns[i]),
            'w_final_mev': float(p.w_end[i]),
            'w_extraction_mev': float(p.w_last[i])
            if p.fate[i] == EXTRACTED
            else float('nan'),
            'behind_rf': bool(p.flagged[i]),
            'r_at_deflector_azimuth_m': float(p.r_deflector[i]),
            'phase_by_bin_deg': [float(v) for v in by_bin[i]],
        }
        for i in range(N)
    ]
    bl = lions_m.beam
    curves = [
        {
            **base,
            'bin': j,
            'r_m': float(t0_beam.probe_r[j]),
            't0_probe': float(t0_beam.probe[j]),
            't0_phase_deg': float(t0_beam.phase[j]),
            'lions_probe': float(bl.probe[j]),
            'lions_phase_deg': float(bl.phase[j]),
        }
        for j in range(mc.PROBE_BINS)
    ]
    with open(path, 'wb') as f:
        pickle.dump((verdict, particles, curves), f)
    shutil.rmtree(work)
    print(
        row['tune'],
        row['method'],
        row['condition'],
        f'{wall:.0f} s',
        'lions in spec',
        verdict['lions_in_spec'],
        'transmission',
        verdict['t0_transmission'],
        verdict['lions_transmission'],
        flush=True,
    )
    return row['tune']


def classes(s, fate):
    """Count, summed and median end energy [MeV] of a fate class."""
    w = s['w_end'][s['fate'] == fate]
    return {
        f'n_{fate}': len(w),
        f'w_sum_{fate}_mev': float(w.sum()),
        f'w_median_{fate}_mev': float(np.median(w))
        if len(w)
        else float('nan'),
    }


def limit_one(job):
    """Run one tune to the longest turn limit; keep rows in parts/."""
    row, mapdir, parts = job
    path = parts / f"{row['tune']:03d}.pkl"
    if path.exists():
        return row['tune']
    mp = Map(mapdir)
    work = mp.workdir(parts / f"work_{row['tune']:03d}")
    ref = reference_orbit(work, mp)
    m = of_tune(row)
    began = time.perf_counter()
    p = Pass(work, mp, ref, m.knobs, m.errors, m.draws, factors=LIMIT_FACTORS)
    wall = time.perf_counter() - began
    base = {
        k: row[k] for k in ('tune', 'method', 'training_seed', 'condition')
    }
    tunes_, particles = [], []
    for f in LIMIT_FACTORS:
        s = p.snaps[f]
        with_c = judged(
            row['condition'], m.knobs, row['field_error'], p.beam_at(s, True)
        )
        without = judged(
            row['condition'], m.knobs, row['field_error'], p.beam_at(s, False)
        )
        b = with_c.beam
        tunes_.append(
            {
                **base,
                'factor': f,
                'limit_turns': p.base * f,
                'wall_s': wall,
                'turn_max': int(s['turns'].max()),
                'transmission': float(b.transmission),
                'extracted_ua': float(b.extracted_ua),
                'w_mean': float(b.w_mean),
                'w_spread': float(b.w_spread),
                'lost_kw_with_circulating': float(b.lost_kw),
                'lost_kw_without_circulating': float(without.beam.lost_kw),
                'in_spec_with_circulating': bool(with_c.ready()),
                'in_spec_without_circulating': bool(without.ready()),
                'terms_with_circulating': json.dumps(with_c.spec_terms()),
                'terms_without_circulating': json.dumps(without.spec_terms()),
                **{
                    k: v
                    for fate in (EXTRACTED, Z_APERTURE, OUT_OF_MAP, TURN_LIMIT)
                    for k, v in classes(s, fate).items()
                },
            }
        )
        particles += [
            {
                **base,
                'factor': f,
                'particle': i,
                'phi_start_deg': float(np.degrees(p.phi[i])),
                'z0_m': float(p.z0[i]),
                'fate': str(s['fate'][i]),
                'w_end_mev': float(s['w_end'][i]),
                'turns': int(s['turns'][i]),
                'r_end_m': float(s['r'][i]),
            }
            for i in range(N)
        ]
    with open(path, 'wb') as out:
        pickle.dump((tunes_, particles), out)
    shutil.rmtree(work)
    print(
        row['tune'],
        row['method'],
        row['condition'],
        f'{wall:.0f} s',
        [
            (
                r['factor'],
                r['n_turn_limit'],
                r['in_spec_with_circulating'],
                r['in_spec_without_circulating'],
            )
            for r in tunes_
        ],
        flush=True,
    )
    return row['tune']


def limit_study(outdir, mp, workers, reference):
    """Run the turn-limit study over all tunes; write and print it."""
    rows = pq.read_table(outdir / 'case_v_tunes.parquet').to_pylist()
    parts = outdir / 'parts'
    parts.mkdir(parents=True, exist_ok=True)
    with multiprocessing.get_context('spawn').Pool(workers) as pool:
        list(
            pool.imap_unordered(limit_one, [(r, mp.dir, parts) for r in rows])
        )
    done = [
        pickle.loads((parts / f"{r['tune']:03d}.pkl").read_bytes())
        for r in rows
    ]
    v = [r for d in done for r in d[0]]
    units = (
        'lost_kw kW, w_mean w_* MeV, w_spread 1, extracted_ua uA, '
        'limit_turns turns, r_end_m m, phi_start_deg deg, z0_m m. '
        'turn_limit is a particle still circulating'
    )
    write(
        v,
        outdir / 'case_v_limit_verdict.parquet',
        'turn limit study per tune',
        mp,
        units,
    )
    write(
        [p for d in done for p in d[1]],
        outdir / 'case_v_limit_particles.parquet',
        'turn limit study per particle',
        mp,
        units,
    )
    shutil.rmtree(parts)
    was = {
        r['tune']: r
        for r in pq.read_table(
            reference / 'case_v_verdict.parquet'
        ).to_pylist()
    }
    one_x = [r for r in v if r['factor'] == 1]
    print(
        '1x against case_v: max |transmission diff|',
        max(
            abs(r['transmission'] - was[r['tune']]['lions_transmission'])
            for r in one_x
        ),
        'max |lost diff| kW',
        max(
            abs(
                r['lost_kw_with_circulating'] - was[r['tune']]['lions_lost_kw']
            )
            for r in one_x
        ),
        'verdict differs',
        sum(
            r['in_spec_with_circulating'] != was[r['tune']]['lions_in_spec']
            for r in one_x
        ),
    )
    for key in ('with_circulating', 'without_circulating'):
        by = {
            f: {r['tune']: r for r in v if r['factor'] == f}
            for f in LIMIT_FACTORS
        }
        print(
            key,
            'in spec',
            {
                f: sum(r[f'in_spec_{key}'] for r in by[f].values())
                for f in LIMIT_FACTORS
            },
            'of',
            len(rows),
        )
        for quantity in (f'lost_kw_{key}', 'transmission'):
            change = [
                abs(by[4][k][quantity] - by[2][k][quantity])
                / max(by[4][k][quantity], 1e-12)
                for k in by[4]
            ]
            print(
                ' ',
                quantity,
                'relative change 2x to 4x: median',
                float(np.median(change)),
                'max',
                float(np.max(change)),
            )
        print(
            '  tunes changing verdict 2x to 4x',
            sum(
                by[2][k][f'in_spec_{key}'] != by[4][k][f'in_spec_{key}']
                for k in by[4]
            ),
        )
    for f in LIMIT_FACTORS:
        s = [r for r in v if r['factor'] == f]
        print(
            f,
            'x:',
            {
                c: sum(r[f'n_{c}'] for r in s)
                for c in (EXTRACTED, Z_APERTURE, OUT_OF_MAP, TURN_LIMIT)
            },
            'of',
            N * len(s),
            'median lost kW with',
            float(np.median([r['lost_kw_with_circulating'] for r in s])),
            'without',
            float(np.median([r['lost_kw_without_circulating'] for r in s])),
        )


def verdict(outdir, mp, workers, start_mev, extraction, phase_rms_deg):
    """Run every tune of the harmonic through LIONS; write tables."""
    rows = [
        r
        for r in pq.read_table(
            outdir / 'case_v_tunes.parquet', use_threads=False
        ).to_pylist()
        if r['harmonic'] == case_p.HARMONIC
    ]
    parts = outdir / 'parts'
    parts.mkdir(parents=True, exist_ok=True)
    with multiprocessing.get_context('spawn').Pool(workers) as pool:
        list(
            pool.imap_unordered(
                one,
                [
                    (r, mp.dir, parts, start_mev, extraction, phase_rms_deg)
                    for r in rows
                ],
            )
        )
    done = [
        pickle.loads((parts / f"{r['tune']:03d}.pkl").read_bytes())
        for r in rows
    ]
    write(
        [d[0] for d in done],
        outdir / 'case_v_verdict.parquet',
        'verdict per tune',
        mp,
        'transmission 1, extracted_ua uA, w_mean MeV, w_spread 1, '
        'lost_kw kW, x y mm, rf_mhz MHz, b0_t T, dee_kv kV, '
        'phase_deg deg, wall_s s',
        phase_rms_deg,
    )
    write(
        [p for d in done for p in d[1]],
        outdir / 'case_v_particles.parquet',
        'one row per particle',
        mp,
        'rf_mhz MHz, b0_t T, dee_kv kV, phase_deg deg, '
        'phi_start_deg deg, phi_first_deg deg (phase at the first '
        'theta = 0 crossing), z0_m m, w_final_mev MeV, '
        'w_extraction_mev MeV, phase_by_bin_deg deg per radial bin of '
        'machine.run',
        phase_rms_deg,
    )
    write(
        [c for d in done for c in d[2]],
        outdir / 'case_v_curves.parquet',
        'probe and phase against radius',
        mp,
        'r_m m, probe particle-turns per macroparticle, phase deg',
        phase_rms_deg,
    )
    shutil.rmtree(parts)
    v = [d[0] for d in done]
    for method in sorted({r['method'] for r in v}):
        s = [r for r in v if r['method'] == method]
        print(
            method,
            len(s),
            'tunes; T0 in spec',
            sum(r['t0_in_spec'] for r in s),
            'LIONS',
            sum(r['lions_in_spec'] for r in s),
            'LIONS with the behind-RF rule',
            sum(r['rule_in_spec'] for r in s),
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'step', choices=('tunes', 'tunes_table', 'check', 'verdict', 'limit')
    )
    parser.add_argument(
        '--outroot', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_v'
    )
    parser.add_argument(
        '--mapdir', type=Path, default=ROOT / 'benchmarks' / 'runs' / 'case_p'
    )
    parser.add_argument(
        '--reference',
        type=Path,
        default=ROOT / 'benchmarks' / 'runs' / 'case_t_drift',
    )
    parser.add_argument('--start-mev', type=float)
    parser.add_argument('--phase-rms-deg', type=float)
    parser.add_argument(
        '--extraction', choices=('mean', 'energy'), default='mean'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=int(os.environ.get('SLURM_CPUS_PER_TASK', 1)),
    )
    args = parser.parse_args()
    if args.phase_rms_deg is not None and args.step != 'verdict':
        parser.error('--phase-rms-deg applies to verdict only')
    outdir = args.outroot.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    mp_ = Map(args.mapdir, args.start_mev, args.extraction)
    if args.step == 'tunes_table':
        sys.exit(
            1 if tunes_from_table(outdir, args.reference.resolve(), mp_) else 0
        )
    elif args.step == 'tunes':
        sys.exit(
            1
            if tunes(outdir, args.reference.resolve(), mp_, args.workers)
            else 0
        )
    elif args.step == 'check':
        check(outdir, mp_)
    elif args.step == 'limit':
        limit_study(outdir, mp_, args.workers, args.reference.resolve())
    else:
        verdict(
            outdir,
            mp_,
            args.workers,
            args.start_mev,
            args.extraction,
            args.phase_rms_deg,
        )
