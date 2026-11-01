"""One-turn maps at the theta = 0 section, closed orbits and tunes, for LIONS, OPAL-cycl and RF-Track.

Section state: (r, pr, z, pz) at y = 0, x > 0, momenta in beta gamma; the total momentum p is fixed.
Each code returns full-step trajectories (t, pos, mom); the section crossing is interpolated in time
with a 4-point Lagrange polynomial.
"""
import concurrent.futures
import contextlib
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import RF_Track as rft

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'benchmarks'))

import lions  # noqa: E402
import lions.cli  # noqa: E402
import lions.pusher  # noqa: E402
import synthetic as syn  # noqa: E402
import convert  # noqa: E402

OPAL_IMAGE = 'lions-bench/opal:2024.1.0'


def lagrange4(t, f, t0):
    """Value at t0 of the cubic through 4 samples (t, f), f of shape (4, ...)."""
    out = 0.0
    for i in range(4):
        w = np.prod([(t0 - t[j]) / (t[i] - t[j]) for j in range(4) if j != i])
        out = out + w * f[i]
    return out


def crossings(t, pos, mom):
    """Section states (r, pr, z, pz) and times at each upward crossing of y = 0 with x > 0."""
    y = pos[:, 1]
    idx = np.nonzero((y[:-1] < 0) & (y[1:] >= 0) & (pos[:-1, 0] > 0))[0]
    out = []
    for k in idx:
        if k < 1 or k + 2 >= len(t):
            continue
        s = slice(k - 1, k + 3)
        ts, ps, ms = t[s], pos[s], mom[s]
        tc = t[k] - y[k] * (t[k + 1] - t[k]) / (y[k + 1] - y[k])
        for _ in range(20):
            yc = lagrange4(ts, ps[:, 1], tc)
            h = 1e-6 * (t[k + 1] - t[k])
            dy = (lagrange4(ts, ps[:, 1], tc + h) - lagrange4(ts, ps[:, 1], tc - h)) / (2 * h)
            step = yc / dy
            tc -= step
            if abs(step) < 1e-15 * (t[k + 1] - t[k]) + 1e-300:
                break
        p = lagrange4(ts, ps, tc)
        m = lagrange4(ts, ms, tc)
        out.append((tc, p[0], m[0], p[2], m[2]))
    return np.array(out)


def initial_vectors(states, p):
    """Section states (n, 4) to positions and momenta (n, 3) at y = 0, moving towards +y."""
    s = np.atleast_2d(states)
    r, pr, z, pz = s.T
    pt = np.sqrt(p**2 - pr**2 - pz**2)
    return (np.column_stack([r, np.zeros_like(r), z]),
            np.column_stack([pr, pt, pz]))


class LionsTracker:
    """LIONS in a working directory holding bmap.dat, pot1.bin, pot2.dat. Many particles per run.

    w_mev: kinetic energy per nucleon; fhf_hz: RF frequency at harmonic 1 (one turn = steps_per_turn steps);
    species: (atomic mass in u, mass number, charge state), protons by default.
    interpolation: the LIONS deck option ('bilinear', the default, or 'hermite').
    """

    def __init__(self, workdir, w_mev, fhf_hz, steps_per_turn, species=(syn.HYDROGEN, 1, 1),
                 interpolation='bilinear'):
        self.workdir, self.w, self.fhf, self.n, self.species = Path(workdir), w_mev, fhf_hz, steps_per_turn, species
        self.tail = () if interpolation == 'bilinear' else ('no no', f'interpolation={interpolation}')

    def track(self, states, turns):
        """Full-step trajectories (t, pos, mom) of shape (steps, n, 3), momenta at full steps."""
        states = np.atleast_2d(states)
        n_part = len(states)
        with contextlib.chdir(self.workdir):
            syn.write_deck('deck.dat', fhf=self.fhf, r0=states[0, 0], w0=self.w, turns=int(np.ceil(turns)),
                           npaspt=self.n, nbrion=n_part, species=self.species, tail=self.tail)
            d = lions.State()
            d.run.swap_bytes = 'n'
            quiet = contextlib.redirect_stdout(io.StringIO())
            quiet.__enter__()
            lions.load(d, 'deck.dat')
            b = d.bunch
            pos, mom = initial_vectors(states, b.p[0])
            b.x[:], b.y[:], b.z[:] = pos.T
            # LIONS treats the deck momentum as p(-1/2): step it back half a step from p(0).
            b.px[:], b.py[:], b.pz[:] = mom.T
            b.gamma[:] = np.sqrt(1 + (mom**2).sum(axis=1))
            lions.cli.check_box(b, d.field, d.box)   # r and the in-map mask for compute_fields
            lions.cli.compute_fields(d)
            half_kick(d, -0.5)
            rows = []

            def record(state):
                # p(k) = half Boris kick of p(k - 1/2) with the field at x(k): the inverse of the
                # start, so the recorded map is conjugate to the Boris step map.
                bb = state.bunch
                kept = bb.px.copy(), bb.py.copy(), bb.pz.copy(), bb.gamma.copy()
                lions.cli.check_box(bb, state.field, state.box)
                lions.cli.compute_fields(state)
                half_kick(state, 0.5)
                rows.append(np.stack([np.column_stack([bb.x, bb.y, bb.z]),
                                      np.column_stack([bb.px, bb.py, bb.pz])]))
                bb.px[:], bb.py[:], bb.pz[:], bb.gamma[:] = kept

            try:
                lions.cli.track(d, on_step=record)
            finally:
                quiet.__exit__(None, None, None)
        a = np.array(rows)                      # (steps, 2, n, 3)
        t = np.arange(len(a)) * d.run.dt
        return t, a[:, 0], a[:, 1]


def half_kick(d, fraction):
    """Boris momentum update over fraction x dt with the current fields; positions unchanged."""
    b = d.bunch
    saved = b.x.copy(), b.y.copy(), b.z.copy()
    coefs = b.boris_b_coef.copy(), b.boris_e_coef.copy()
    b.boris_b_coef *= fraction
    b.boris_e_coef *= fraction
    lions.pusher.boris_push(b, d.run, d.traps)
    b.boris_b_coef, b.boris_e_coef = coefs
    b.x[:], b.y[:], b.z[:] = saved


class OpalTracker:
    """OPAL-cycl single-particle runs in Docker, CARBONCYCL map written from a LIONS FieldMap."""

    def __init__(self, workdir, field, mass_kg, p, f_rev_hz, steps_per_turn, integrator='LF2'):
        self.workdir, self.n, self.integrator, self.p = Path(workdir).resolve(), steps_per_turn, integrator, p
        self.mass_gev = mass_kg * syn.C**2 / syn.QE / 1e9
        self.f_rev_mhz = f_rev_hz / 1e6
        self.workdir.mkdir(parents=True, exist_ok=True)
        convert.write_opal_carbon(self.workdir / 'bmap_carbon.dat', field)
        (self.workdir / 'part.dat').write_text('1\n0 0 0 0 0 0\n')
        self.workers = int(subprocess.run(['docker', 'info', '--format', '{{.NCPU}}'], capture_output=True,
                                          text=True, check=True).stdout)

    def track_one(self, state, turns, tag='run'):
        r, pr, z, pz = state
        p = self.p
        stem = f'{tag}_{self.integrator}_{self.n}'
        (self.workdir / f'{stem}.in').write_text(f'''OPTION, ECHO=FALSE;
OPTION, PSDUMPFREQ=100000000;
OPTION, SPTDUMPFREQ=1;
OPTION, VERSION=20240;
Title, string="LIONS benchmark";
cyc: CYCLOTRON, TYPE=CARBONCYCL, CYHARMON=1, PHIINIT=0.0, PRINIT={pr:.17g}, RINIT={r * 1e3:.17g},
     ZINIT={z * 1e3:.17g}, PZINIT={pz:.17g}, SYMMETRY=1.0, RFFREQ={self.f_rev_mhz:.17g},
     FMAPFN="bmap_carbon.dat", MINZ=-1000, MAXZ=1000, MINR=0, MAXR=10000;
l1: LINE=(cyc);
dist: DISTRIBUTION, TYPE=FROMFILE, FNAME="part.dat";
fs: FIELDSOLVER, FSTYPE=NONE, MX=8, MY=8, MT=8, PARFFTX=true, PARFFTY=true, PARFFTT=false,
    BCFFTX=open, BCFFTY=open, BCFFTT=open;
beam1: BEAM, MASS={self.mass_gev:.17g}, CHARGE=1.0, PC={p * self.mass_gev:.17g}, NPART=1,
       BCURRENT=0, BFREQ={self.f_rev_mhz:.17g};
TRACK, LINE=l1, BEAM=beam1, MAXSTEPS={int(turns * self.n)}, STEPSPERTURN={self.n},
       TIMEINTEGRATOR="{self.integrator}";
RUN, METHOD="CYCLOTRON-T", BEAM=beam1, FIELDSOLVER=fs, DISTRIBUTION=dist;
ENDTRACK;
STOP;
''')
        subprocess.run(['docker', 'run', '--rm', '-v', f'{self.workdir}:/data', OPAL_IMAGE,
                        'bash', '-c', f'opal {stem}.in > {stem}.log 2>&1'], check=True)
        pos, mom = convert.read_opal_orbit(self.workdir / f'{stem}-trackOrbit.dat')
        t = np.arange(len(pos)) / (self.f_rev_mhz * 1e6 * self.n)
        return t, pos, mom

    def track(self, states, turns):
        with concurrent.futures.ThreadPoolExecutor(self.workers) as pool:
            out = list(pool.map(lambda a: self.track_one(a[1], turns, tag=f'p{a[0]}'),
                                enumerate(np.atleast_2d(states))))
        return stack_runs(out)


class _UserField(rft.UserField):
    """RF-Track field from b(x, y, z) -> (Bx, By, Bz) in T, LIONS frame and metres.

    RF-Track calls get_field in the element frame, z from 0 at the entrance; the element is centred on Z = 0.
    """

    def __init__(self, length, b):
        super().__init__(length)
        self.b, self.half_mm = b, length * 1e3 / 2

    def get_field(self, x, y, z, t):
        pos, _ = convert.from_rftrack([x, 0.0, y, 0.0, z - self.half_mm, 0.0], 1.0)
        return [0.0, 0.0, 0.0], list(convert.to_rftrack_field(self.b(*pos[0])))


class RFTrackTracker:
    """RF-Track `Volume` around a static magnetic field b(x, y, z), one particle per run.

    The volume spans |y| < half_size (RF-Track Z); positions are sampled samples_per_turn times per
    nominal period.
    """

    def __init__(self, b, half_size, mass_kg, p, period_s, samples_per_turn, tol, algorithm='rkf45'):
        self.mass_mev = mass_kg * syn.C**2 / syn.QE / 1e6
        self.p, self.tol, self.algorithm = p, tol, algorithm
        self.period_mm = period_s * syn.C * 1e3
        self.dt_mm = self.period_mm / samples_per_turn
        self.volume = rft.Volume()
        self.volume.add(_UserField(2 * half_size, b), 0.0, 0.0, -half_size)

    def track_one(self, state, turns):
        pos, mom = initial_vectors(state, self.p)
        opt = rft.TrackingOptions()
        opt.dt_mm = opt.tt_dt_mm = self.dt_mm
        opt.t_max_mm = turns * self.period_mm
        opt.odeint_algorithm = self.algorithm
        opt.odeint_epsabs = opt.odeint_epsrel = self.tol
        self.volume.track(rft.Bunch6dT(self.mass_mev, 1.0, 1.0, convert.to_rftrack(pos, mom, self.mass_mev)), opt)
        tt = self.volume.get_transport_table('%t %mean_X %mean_Px %mean_Y %mean_Py %mean_Z %mean_Pz')
        pos, mom = convert.from_rftrack(tt[:, 1:], self.mass_mev)
        return tt[:, 0] * 1e-3 / syn.C, pos, mom

    def track(self, states, turns):
        return stack_runs([self.track_one(s, turns) for s in np.atleast_2d(states)])


class BorisTracker:
    """The LIONS Boris step (magnetic field only) in numpy, with a pluggable field(x, y, z) -> (Bx, By, Bz) in T.

    Same conventions as LionsTracker: momenta in beta gamma at half steps, started by a backward half kick at x0,
    recorded at full steps after a forward half kick, so one-turn maps are comparable with LIONS step for step.
    """

    def __init__(self, field, mass_kg, charge_c, p, f_rev_hz, steps_per_turn):
        self.field, self.p, self.n = field, p, steps_per_turn
        self.dt = 1 / (f_rev_hz * steps_per_turn)
        self.bcoef = charge_c * self.dt / 2 / mass_kg      # LIONS boris_b_coef
        self.c_dt = syn.C * self.dt

    def kick(self, pos, mom, fraction):
        """Boris rotation of mom over fraction x dt with the field at pos (LIONS boris_push with E = 0)."""
        bx, by, bz = self.field(*pos.T)
        gamma = np.sqrt(1 + (mom**2).sum(axis=1))
        vec = fraction * self.bcoef / gamma
        a1, a2, a3 = vec * bx, vec * by, vec * bz
        v1, v2, v3 = a1 * a1, a2 * a2, a3 * a3
        vec = 1.0 / (1.0 + v1 + v2 + v3)
        s12 = a1 * a2
        s21, s12 = s12 - a3, s12 + a3
        s13 = a1 * a3
        s31, s13 = s13 + a2, s13 - a2
        s23 = a2 * a3
        s32, s23 = s23 - a1, s23 + a1
        px, py, pz = mom.T
        v4, v5, v6 = vec * (0.0 + 2.0 * px), vec * (0.0 + 2.0 * py), vec * (0.0 + 2.0 * pz)
        return np.column_stack([px - (v2 + v3) * v4 + s12 * v5 + s13 * v6,
                                py + s21 * v4 - (v1 + v3) * v5 + s23 * v6,
                                pz + s31 * v4 + s32 * v5 - (v1 + v2) * v6])

    def track(self, states, turns):
        pos, mom = initial_vectors(states, self.p)
        half = self.kick(pos, mom, -0.5)
        steps = int(np.ceil(turns)) * self.n
        out_pos, out_mom = [pos], [mom]
        for _ in range(steps):
            half = self.kick(pos, half, 1.0)
            pos = pos + self.c_dt * half / np.sqrt(1 + (half**2).sum(axis=1))[:, None]
            out_pos.append(pos)
            out_mom.append(self.kick(pos, half, 0.5))
        return np.arange(steps + 1) * self.dt, np.array(out_pos), np.array(out_mom)


class RK4Tracker:
    """Classical RK4, fixed step, on dx/dt = c p / gamma and dp/dt = (q / m) (p / gamma) x B, p in beta gamma.

    Same field call, time step and full-step recording as BorisTracker.
    """

    def __init__(self, field, mass_kg, charge_c, p, f_rev_hz, steps_per_turn):
        self.field, self.p, self.n = field, p, steps_per_turn
        self.dt = 1 / (f_rev_hz * steps_per_turn)
        self.qm = charge_c / mass_kg

    def rate(self, pos, mom):
        gamma = np.sqrt(1 + (mom**2).sum(axis=1))[:, None]
        b = np.column_stack(self.field(*pos.T))
        return syn.C * mom / gamma, self.qm * np.cross(mom, b) / gamma

    def track(self, states, turns):
        pos, mom = initial_vectors(states, self.p)
        steps = int(np.ceil(turns)) * self.n
        h = self.dt
        out_pos, out_mom = [pos], [mom]
        for _ in range(steps):
            k1x, k1p = self.rate(pos, mom)
            k2x, k2p = self.rate(pos + h / 2 * k1x, mom + h / 2 * k1p)
            k3x, k3p = self.rate(pos + h / 2 * k2x, mom + h / 2 * k2p)
            k4x, k4p = self.rate(pos + h * k3x, mom + h * k3p)
            pos = pos + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
            mom = mom + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
            out_pos.append(pos)
            out_mom.append(mom)
        return np.arange(steps + 1) * self.dt, np.array(out_pos), np.array(out_mom)


def stack_runs(out):
    """Single-particle runs [(t, pos, mom)] to (t, pos, mom) of shape (steps, n, 3), cut to the shortest."""
    steps = min(len(o[0]) for o in out)
    return (out[0][0][:steps], np.stack([o[1][:steps] for o in out], axis=1),
            np.stack([o[2][:steps] for o in out], axis=1))


def section_series(t, pos, mom):
    """Turn-by-turn section states for each particle: list of arrays (turns, 5)."""
    return [crossings(t, pos[:, i], mom[:, i]) for i in range(pos.shape[1])]


def one_turn(track, states):
    """Map section states (n, 4) through one turn. track(states, turns) -> (t, pos, mom)."""
    t, pos, mom = track(states, 1.2)
    series = section_series(t, pos, mom)
    lost = [i for i, s in enumerate(series) if len(s) == 0]
    if lost:
        raise RuntimeError(f'particles {lost} never crossed the section (lost from the map or invalid start)')
    return np.array([s[0, 1:] for s in series])


def closed_orbit(track, guess, iters=8, h=1e-7):
    """Newton on (r, pr) at z = pz = 0, finite-difference Jacobian. Returns state and residual."""
    s = np.array([guess[0], guess[1], 0.0, 0.0])
    for _ in range(iters):
        probes = [s] + [s + h * e for e in np.eye(4)[:2]] + [s - h * e for e in np.eye(4)[:2]]
        out = one_turn(track, np.array(probes))
        f0 = out[0, :2] - s[:2]
        jac = np.column_stack([(out[1 + i, :2] - out[3 + i, :2]) / (2 * h) for i in range(2)]) - np.eye(2)
        delta = np.linalg.solve(jac, -f0)
        s[:2] += delta
        if np.max(np.abs(delta)) < 1e-13:
            break
    return s, np.max(np.abs(f0))


def tune(signal):
    """Fractional tune of a turn-by-turn signal: Hann-windowed DFT peak, refined by golden section."""
    x = np.asarray(signal) - np.mean(signal)
    n = len(x)
    w = 0.5 * (1 - np.cos(2 * np.pi * np.arange(n) / n))
    k = np.arange(n)

    def amp(q):
        return -abs(np.sum(w * x * np.exp(-2j * np.pi * q * k)))

    spec = np.abs(np.fft.rfft(w * x))
    q0 = np.argmax(spec[1:]) + 1
    a, b = (q0 - 1) / n, (q0 + 1) / n
    g = (np.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    for _ in range(100):
        if amp(c) < amp(d):
            b = d
        else:
            a = c
        c, d = b - g * (b - a), a + g * (b - a)
    return 0.5 * (a + b)


S4 = np.kron(np.eye(2), np.array([[0.0, 1.0], [-1.0, 0.0]]))


def jacobian(track, s, h_r, h_p):
    """4x4 one-turn Jacobian at the section state s, central differences (steps h_r in m, h_p in p)."""
    h = np.array([h_r, h_p, h_r, h_p])
    steps = np.diag(h)
    out = one_turn(track, np.vstack([s + steps, s - steps]))
    return ((out[:4] - out[4:]) / (2 * h)[:, None]).T


def linear_tunes(m):
    """Fractional tunes from the r and z diagonal blocks: cos(2 pi Q) = trace / 2."""
    return tuple(np.arccos(np.clip(np.trace(m[i:i + 2, i:i + 2]) / 2, -1, 1)) / (2 * np.pi)
                 for i in (0, 2))


def symplectic_error(m):
    """max |M^T S M - S| and |det M - 1|."""
    return np.max(np.abs(m.T @ S4 @ m - S4)), abs(np.linalg.det(m) - 1)
