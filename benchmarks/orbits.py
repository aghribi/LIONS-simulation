"""One-turn maps at the theta = 0 section, closed orbits and tunes, for LIONS and OPAL-cycl.

Section state: (r, pr, z, pz) at y = 0, x > 0, momenta in beta gamma; the total momentum p is fixed.
Each code returns full-step trajectories (t, pos, mom); the section crossing is interpolated in time
with a 4-point Lagrange polynomial.
"""
import contextlib
import io
import subprocess
import sys
from pathlib import Path

import numpy as np

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
    """LIONS in a working directory holding bmap.dat, pot1.bin, pot2.dat. Many particles per run."""

    def __init__(self, workdir, w_mev, b_ref, steps_per_turn):
        self.workdir, self.w, self.b_ref, self.n = Path(workdir), w_mev, b_ref, steps_per_turn

    def track(self, states, turns):
        """Full-step trajectories (t, pos, mom) of shape (steps, n, 3), momenta at full steps."""
        states = np.atleast_2d(states)
        n_part = len(states)
        with contextlib.chdir(self.workdir):
            syn.write_deck('deck.dat', fhf=syn.revolution_frequency(self.w, self.b_ref),
                           r0=states[0, 0], w0=self.w, turns=int(np.ceil(turns)), npaspt=self.n,
                           nbrion=n_part)
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
        out = [self.track_one(s, turns, tag=f'p{i}') for i, s in enumerate(np.atleast_2d(states))]
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
