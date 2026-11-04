"""Figures, tables and numbers of the lecture notes, from the benchmark
tables in benchmarks/runs.

usage: uv run --group bench --with pandas python docs/lecture/figures.py
Writes docs/lecture/gen/: one PDF per figure, numbers.tex (one macro per
quoted value), tab_*.tex, sources.tex.
"""

import hashlib
import sys
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('pdf')
import matplotlib.pyplot as plt  # noqa: E402 (after sys.path)

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'benchmarks'))
import fieldmodels  # noqa: E402 (after sys.path)
from provenance import CODE  # noqa: E402 (after sys.path)

RUNS = ROOT / 'benchmarks' / 'runs'
GEN = Path(__file__).resolve().parent / 'gen'
NUMBERS, SOURCES = {}, {}
COLORS = {
    'lions': '#1f5f99',
    'opal_LF2': '#c8553d',
    'opal_RK4': '#8a3ffc',
    'rf_track': '#2e8b57',
    'spline': '#e6a100',
    'opal': '#c8553d',
    'rftrack': '#2e8b57',
}


def log_ticks(ax, ticks):
    """Label the x axis of ax at ticks, no minor tick labels."""
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


plt.rcParams.update(
    {
        'font.size': 9,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'figure.dpi': 150,
        'legend.fontsize': 7.5,
        'savefig.bbox': 'tight',
        'lines.markersize': 4,
    }
)


def psi_band():
    """Energy range [MeV] flagged as the Qr half-integer band in the G2
    step 3 table.
    """
    d = table('case_g2/step3/case_g2_step3.parquet')
    e = d[d.qr_half_integer_band].energy_mev
    return float(e.min()), float(e.max())


def table(path):
    """Read a Parquet table under RUNS; record its sha256 in SOURCES."""
    path = RUNS / path
    SOURCES[str(path.relative_to(ROOT))] = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()
    return pd.read_parquet(path)


def num(name, value, digits=2):
    """\\name expands to value in siunitx scientific form with `digits`
    significant digits.
    """
    if not name.isalpha():
        raise ValueError(f'macro name {name} must be letters only')
    value = float(value)
    if value != 0 and (abs(value) < 1e-2 or abs(value) >= 1e4):
        text = f'{value:.{digits - 1}e}'
    elif abs(value) >= 10**digits:  # e.g. 770, not 7.7e+02
        text = f'{round(value, digits - 1 - int(np.floor(np.log10(abs(value))))):.0f}'  # noqa: E501 (f-string, Python 3.10)
    else:
        text = f'{value:.{digits}g}'
    NUMBERS[name] = rf'\num{{{text}}}'


def raw(name, text):
    """Set the LaTeX macro name to text, as given."""
    NUMBERS[name] = text


def sci(value, digits=2):
    """Format value as a siunitx \\num in scientific form."""
    return rf'\num{{{value:.{digits - 1}e}}}'


def save(fig, name):
    """Save fig as GEN/name.pdf and close it."""
    fig.savefig(GEN / f'{name}.pdf')
    plt.close(fig)


def write_table(name, header, rows):
    """Write GEN/tab_name.tex: a booktabs table of header and rows."""
    lines = [
        r'\begin{tabular}{' + 'l' + 'r' * (len(header) - 1) + '}',
        r'\toprule',
        ' & '.join(header) + r' \\',
        r'\midrule',
        *(' & '.join(row) + r' \\' for row in rows),
        r'\bottomrule',
        r'\end{tabular}',
    ]
    (GEN / f'tab_{name}.tex').write_text('\n'.join(lines) + '\n')


def order(e1, e2):
    """Observed convergence order, log2 of |e1| / |e2|."""
    return np.log2(abs(e1) / abs(e2))


# Illustrations computed here (no benchmark table).


def boris_polygon():
    """Boris in a uniform field, n steps per turn: positions on a circle
    of radius rho sqrt(1 + (w dt/2)^2).
    """
    n, rho = 8, 1.0
    wdt = 2 * np.pi / n
    angle = (
        2 * np.arctan(wdt / 2)
        # Boris rotation of p per step; |p| and the step length rho w dt
        # are exact
    )
    phase = -angle / 2 + angle * np.arange(1, n + 1)
    steps = rho * wdt * np.column_stack([-np.sin(phase), np.cos(phase)])
    xs = np.vstack([[0.0, 0.0], np.cumsum(steps, axis=0)])
    (ax_, ay), (bx, by), (cx, cy) = xs[:3]
    den = 2 * (ax_ * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    s = (ax_**2 + ay**2, bx**2 + by**2, cx**2 + cy**2)
    xs -= (
        np.array(
            [
                s[0] * (by - cy) + s[1] * (cy - ay) + s[2] * (ay - by),
                s[0] * (cx - bx) + s[1] * (ax_ - cx) + s[2] * (bx - ax_),
            ]
        )
        / den
    )  # circumcentre
    phi = np.linspace(0, 2 * np.pi, 400)
    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    ax.plot(
        rho * np.cos(phi),
        rho * np.sin(phi),
        'k-',
        lw=0.8,
        label=r'exact orbit, radius $\rho$',
    )
    rb = np.hypot(*xs[0])
    ax.plot(
        rb * np.cos(phi),
        rb * np.sin(phi),
        'k:',
        lw=0.8,
        label=r'$\rho\sqrt{1+(\omega\Delta t/2)^2}$',
    )
    ax.plot(
        xs[:, 0],
        xs[:, 1],
        'o-',
        color=COLORS['lions'],
        lw=1,
        label='Boris positions, 8 steps',
    )
    ax.set_aspect('equal')
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.legend(loc='upper right', fontsize=6.5)
    ax.set_xlabel(r'$x/\rho$')
    ax.set_ylabel(r'$y/\rho$')
    save(fig, 'boris_polygon')
    num('BorisPolyRadius', rb / rho, 4)


def azimuth_error():
    """Azimuth recovered from (x, y) = r (cos t, sin t) in float64:
    arccos(x / r) against atan2(y, x).
    """
    theta = np.logspace(-13, -0.2, 2000)
    r = 1.0
    x, y = r * np.cos(theta), r * np.sin(theta)
    rr = np.hypot(x, y)
    e_acos = np.abs(np.arccos(np.clip(x / rr, -1, 1)) - theta) / theta
    e_atan = np.abs(np.arctan2(y, x) - theta) / theta
    eps = np.finfo(float).eps
    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    ax.loglog(
        theta,
        np.maximum(e_acos, 1e-18),
        '.',
        ms=1.5,
        color=COLORS['opal'],
        label=r'$\arccos(x/r)$',
    )
    ax.loglog(
        theta,
        np.maximum(e_atan, 1e-18),
        '.',
        ms=1.5,
        color=COLORS['lions'],
        label=r'atan2$(y, x)$',
    )
    ax.loglog(
        theta,
        np.minimum(eps / theta**2, 1),
        'k--',
        lw=0.7,
        label=r'$\min(1, \varepsilon/\theta^2)$',
    )
    ax.set_xlabel(r'$\theta$ [rad]')
    ax.set_ylabel('relative error of the azimuth')
    ax.set_ylim(1e-18, 3)
    ax.legend(loc='lower left')
    save(fig, 'azimuth_error')
    num('AcosErrAtOneMicro', np.interp(1e-6, theta, e_acos))


def field_maps():
    """Median-plane Bz of the two FROM_ORBIT masters."""
    fig, axes = plt.subplots(
        1, 2, figsize=(6.4, 3.0), subplot_kw={'projection': 'polar'}
    )
    for ax, name, label in zip(
        axes,
        ('3sector_master', '10sector_master'),
        ('3-sector, H$_2^+$', '10-sector, p'),
    ):
        path = RUNS / 'case_g1' / name / 'master.h5'
        SOURCES[str(path.relative_to(ROOT))] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        with h5py.File(path) as f:
            bz, a = f['bz'][()], dict(f.attrs)
        n_sec = int(a['sector_count'])
        r = a['ri_m'] + a['rstep_m'] * np.arange(bz.shape[0])
        th = a['theta_step_rad'] * np.arange(bz.shape[1] * n_sec + 1)
        full = np.hstack([bz] * n_sec + [bz[:, :1]])
        m = ax.pcolormesh(
            th,
            r,
            full,
            shading='auto',
            cmap='RdBu_r',
            vmin=-np.abs(full).max(),
            vmax=np.abs(full).max(),
            rasterized=True,
        )
        ax.set_ylim(r[0], r[-1])
        ax.set_rorigin(
            r[0] - 0.6 * (r[-1] - r[0])
        )  # annulus of the map, inner radius not to scale
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_title(
            f'{label}, $r$ from {r[0]:.1f} to {r[-1]:.1f} m', fontsize=8
        )
        fig.colorbar(m, ax=ax, shrink=0.7, pad=0.08, label=r'$B_z$ [T]')
    save(fig, 'field_maps')


def interpolation_kinks():
    """Bz and dBz/dtheta along theta on the 3-sector master at one
    radius: 1 deg nodes, bilinear against Hermite.
    """
    path = RUNS / 'case_g1' / '3sector_master' / 'master.h5'
    with h5py.File(path) as f:
        grids = {k: f[k][()] for k in ('bz', 'dbz_dr', 'dbz_dtheta')}
        attrs = dict(f.attrs)
    attrs['sector_count'] = int(attrs['sector_count'])
    every = 3
    r, th, bz, dbr, dbt = fieldmodels.node_arrays(grids, attrs, every)
    surface = fieldmodels.hermite_surface(r, th, bz, dbr, dbt)
    i = np.searchsorted(r, 1.2)
    rr = r[i] + 0.37 * (r[1] - r[0])
    t = np.linspace(np.deg2rad(24), np.deg2rad(29), 1200)
    h_bz, _, h_dt = surface(np.full_like(t, rr), t)
    # bilinear on the same nodes (LIONS weights), and its theta
    # derivative per cell
    dth = th[1] - th[0]
    j = (t / dth).astype(int)
    v = (t - th[j]) / dth
    u = (rr - r[i]) / (r[1] - r[0])
    bil = (1 - u) * ((1 - v) * bz[i, j] + v * bz[i, j + 1]) + u * (
        (1 - v) * bz[i + 1, j] + v * bz[i + 1, j + 1]
    )
    bil_dt = (
        (1 - u) * (bz[i, j + 1] - bz[i, j])
        + u * (bz[i + 1, j + 1] - bz[i + 1, j])
    ) / dth
    # reference: Hermite surface on the 1/3 deg master nodes
    rm, thm, bzm, _, dbtm = fieldmodels.node_arrays(grids, attrs, 1)
    exact = fieldmodels.hermite_surface(
        rm, thm, bzm, *fieldmodels.node_arrays(grids, attrs, 1)[3:]
    )
    e_bz, _, e_dt = exact(np.full_like(t, rr), t)
    deg = np.rad2deg(t)
    fig, axes = plt.subplots(2, 1, figsize=(4.6, 4.0), sharex=True)
    axes[0].plot(
        deg,
        -(bil - e_bz),
        '-',
        color=COLORS['opal'],
        lw=1,
        label='bilinear, 1 deg nodes',
    )
    axes[0].plot(
        deg,
        -(h_bz - e_bz),
        '--',
        color=COLORS['lions'],
        lw=1,
        label='Hermite, 1 deg nodes',
    )
    axes[0].set_ylabel(r'$B_z - B_z^{\rm ref}$ [T]')
    axes[0].legend(fontsize=6.5)
    axes[1].plot(
        deg, -e_dt, 'k-', lw=0.8, label='reference (1/3 deg master, Hermite)'
    )
    axes[1].plot(deg, -bil_dt, '-', color=COLORS['opal'], lw=1)
    axes[1].plot(deg, -h_dt, '--', color=COLORS['lions'], lw=1)
    axes[1].set_ylabel(r'$\partial B_z/\partial\theta$ [T/rad]')
    axes[1].set_xlabel(r'$\theta$ [deg]')
    axes[1].legend(fontsize=6.5)
    for x in np.rad2deg(th[(th > t[0]) & (th < t[-1])]):
        for ax in axes:
            ax.axvline(x, color='0.6', lw=0.4, zorder=0)
    save(fig, 'interpolation_kinks')
    raw('KinkRadius', f'{rr:.3f}')


def interpolation_orders():
    """Interpolation error bounds for one harmonic sin(2 pi theta /
    lambda) sampled with N nodes per period.
    """
    for n, key in ((18, 'Eighteen'), (54, 'FiftyFour')):
        kh = 2 * np.pi / n
        num(f'LinErr{key}', kh**2 / 8)
        num(f'HermErr{key}', kh**4 / 384)


# Benchmarks.


def case_a():
    """Figure case_a and its numbers, from case_a.parquet."""
    d = table('case_a/case_a.parquet')
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.7))
    fig.subplots_adjust(wspace=0.4)
    for (code, integ), g in d.groupby(['code', 'integrator']):
        if code == 'rf_track':
            continue
        key = 'lions' if code == 'lions' else f'opal_{integ}'
        for ax, metric in zip(axes, ('radius_rel_error', 'phase_rel_error')):
            s = g[g.metric == metric].sort_values('steps_per_turn')
            ax.loglog(
                s.steps_per_turn,
                np.maximum(s.value.abs(), 1e-17),
                'o-',
                color=COLORS[key],
                label=f'{code.upper()} {integ}',
            )
    n = np.array([90, 1440])
    axes[0].loglog(
        n,
        (2 * np.pi / n) ** 2 / 8,
        'k--',
        lw=0.7,
        label=r'$(\omega\Delta t)^2/8$',
    )
    axes[0].loglog(
        n,
        (2 * np.pi) ** 6 / (144 * n**5.0),
        'k:',
        lw=0.7,
        label=r'$(2\pi)^6/(144 n^5)$',
    )
    nn = np.geomspace(90, 1440, 50)
    x = 2 * np.pi / nn
    axes[1].loglog(
        nn,
        1 - 2 * np.arctan(x / 2) / x,
        'k--',
        lw=0.7,
        label=r'$1-2\arctan(\omega\Delta t/2)/\omega\Delta t$',
    )
    axes[0].set_ylabel('relative radius error')
    axes[1].set_ylabel('relative phase error')
    for ax in axes:
        ax.set_xlabel('steps per turn')
        ax.legend(fontsize=6)
    save(fig, 'case_a')
    v = d.set_index(
        ['code', 'integrator', 'steps_per_turn', 'metric']
    ).value.sort_index()
    num(
        'AlionsRadOrder',
        order(
            v['lions', 'boris', 720, 'radius_rel_error'],
            v['lions', 'boris', 1440, 'radius_rel_error'],
        ),
        3,
    )
    num('AlionsRadFourteen', v['lions', 'boris', 1440, 'radius_rel_error'])
    num(
        'AlionsPolygon',
        v['lions', 'boris']
        .xs('radius_vs_boris_polygon', level='metric')
        .abs()
        .max(),
        1,
    )
    num('ARkFourPThreeSixty', abs(v['opal', 'RK4', 360, 'p_rel_change']), 3)
    num(
        'ARkFourRadThreeSixty',
        abs(v['opal', 'RK4', 360, 'radius_rel_error']),
        3,
    )
    num('ARkFourPred', (2 * np.pi) ** 6 / (144 * 360**5), 3)
    num('LeapfrogPred', (1 + (2 * np.pi / 360) ** 2) ** 180 - 1, 4)
    num('ARkFourPhaseFloor', abs(v['opal', 'RK4', 1440, 'phase_rel_error']), 1)
    num(
        'ALFTwoP',
        v['opal', 'LF2'].xs('p_rel_change', level='metric').abs().max(),
        1,
    )


def case_b():
    """Figures case_b0, case_b_scatter, numbers: case_b.parquet."""
    d = table('case_b/case_b.parquet')
    b0 = d[d.case == 'B0']
    exact = b0[b0.code == 'exact'].set_index(['energy_mev', 'metric']).value
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    for (code, integ), g in b0[
        (b0.code != 'exact') & (b0.metric == 'r_eo')
    ].groupby(['code', 'integrator']):
        if code == 'rf_track':
            g = g[g.tolerance == 1e-12]
        err = g.apply(
            lambda s: abs(s.value / exact[s.energy_mev, 'r_eo'] - 1), axis=1
        )
        g = g.assign(err=err)
        if (
            code == 'rf_track'
        ):  # adaptive, sampled 1440 times per turn: one value per energy
            ax.axhline(
                g.err.max(),
                color=COLORS['rf_track'],
                ls='--',
                lw=1,
                label=r'RF-Track rkf45, tol. $10^{-12}$',
            )
            continue
        key = 'lions' if code == 'lions' else f'opal_{integ}'
        agg = g.groupby('steps_per_turn').err.agg(['min', 'max'])
        ax.fill_between(
            agg.index, agg['min'], agg['max'], color=COLORS[key], alpha=0.25
        )
        ax.loglog(
            agg.index,
            agg['max'],
            'o-',
            color=COLORS[key],
            label=f'{code.upper()} {integ}',
        )
    log_ticks(ax, [360, 720, 1440])
    ax.set_xlabel('steps per turn')
    ax.set_ylabel(r'$|r_{eo}/r_{exact}-1|$, 1 to 5 MeV')
    ax.legend(fontsize=6.5)
    save(fig, 'case_b0')
    lions = b0[(b0.code == 'lions') & (b0.metric == 'r_eo')]

    def err(s):
        return abs(s.value / exact[s.energy_mev, 'r_eo'] - 1)

    e = {
        n: lions[lions.steps_per_turn == n].apply(err, axis=1)
        for n in (720, 1440)
    }
    num(
        'BzeroLionsOrderMin',
        min(order(a, b) for a, b in zip(e[720], e[1440])),
        4,
    )
    num('BzeroLionsFourteenMax', e[1440].max())
    tun = b0[
        (b0.code == 'lions')
        & (b0.steps_per_turn == 1440)
        & b0.metric.isin(['qz'])
    ]
    num(
        'BzeroLionsDet',
        b0[
            (b0.code == 'lions')
            & (b0.steps_per_turn == 1440)
            & (b0.metric == 'det_error')
        ]
        .value.abs()
        .max(),
        1,
    )
    opal = b0[
        (b0.code == 'opal')
        & (b0.integrator == 'RK4')
        & (b0.steps_per_turn == 1440)
        & (b0.metric == 'r_eo')
    ]
    num('BzeroOpalOffset', opal.apply(err, axis=1).max(), 1)
    num(
        'BzeroLionsQzErrMax',
        max(
            abs(s.value - exact[s.energy_mev, 'qz']) for _, s in tun.iterrows()
        ),
    )
    lf = b0[(b0.code == 'lions') & (b0.steps_per_turn == 1440)].pivot_table(
        index='energy_mev', columns='metric', values='value'
    )
    num(
        'BzeroFFTAgree',
        max((lf.qr - lf.qr_fft).abs().max(), (lf.qz - lf.qz_fft).abs().max()),
        1,
    )
    num('BzeroParityBz', d[d.metric == 'parity_bz_max_abs'].value.max())
    par = d[d.metric == 'parity_bxy_max_abs'].value
    num('BzeroParityBxy', par.max())

    # case B, 1 MeV, Qr against steps per code, and the RF-Track run on
    # the closed-form field (case B2)
    b = d[(d.case == 'B') & (d.energy_mev == 1.0) & (d.map_dtheta_deg == 1.0)]
    b2 = table('case_b2/case_b2.parquet')
    smooth = b2[(b2.case == 'B-smooth') & (b2.metric == 'qr')]
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    for (code, integ), g in b[b.metric == 'qr'].groupby(
        ['code', 'integrator']
    ):
        if code == 'rf_track':
            continue
        key = 'lions' if code == 'lions' else f'opal_{integ}'
        g = g.sort_values('steps_per_turn')
        ax.plot(
            g.steps_per_turn,
            g.value,
            'o-',
            color=COLORS[key],
            label=f'{code.upper()} {integ}, bilinear map',
        )
    s = smooth.sort_values('steps_per_turn')
    ax.plot(
        s.steps_per_turn,
        s.value,
        's-',
        color=COLORS['rf_track'],
        label='RF-Track rkf45, closed form',
    )
    ax.set_xscale('log', base=2)
    log_ticks(ax, [360, 720, 1440])
    ax.set_xlabel('steps per turn')
    ax.set_ylabel(r'fractional $Q_r$, 1 MeV')
    ax.legend(fontsize=6.5)
    save(fig, 'case_b_scatter')
    q = b[(b.metric == 'qr') & (b.code != 'rf_track')]
    num(
        'BQrSpreadMax',
        q.groupby(['code', 'integrator']).value.agg(np.ptp).max(),
    )
    num('BsmoothQrSpread', np.ptp(s.value))
    num('BLionsQrSpread', np.ptp(q[q.code == 'lions'].value))
    det = b[(b.metric == 'det_error') & (b.code != 'rf_track')].value.abs()
    num('BDetMin', det.min(), 1)
    num('BDetMax', det.max(), 1)
    sp = d[(d.case == 'B') & (d.code == 'lions') & (d.metric == 'qr')]
    ratio = [
        np.ptp(g[g.map_dtheta_deg == 0.5].value)
        / np.ptp(g[g.map_dtheta_deg == 1.0].value)
        for _, g in sp.groupby('energy_mev')
    ]
    num('BHalfDegQrRatioMin', min(ratio))
    num('BHalfDegQrRatioMax', max(ratio), 3)
    # case B2: conditioning det(M_r - I) = 4 sin^2(pi Qr) (LIONS, 1440
    # steps) and r_eo spread across steps, per energy
    full_b = d[(d.case == 'B') & (d.map_dtheta_deg == 1.0)]
    grad = b2[b2.case == 'B-gradient']

    def cond(x):
        return (
            x[
                (x.code == 'lions')
                & (x.steps_per_turn == 1440)
                & (x.metric == 'qr')
            ]
            .set_index('energy_mev')
            .value.pipe(lambda q: 4 * np.sin(np.pi * q) ** 2)
        )

    c_ratio = cond(grad) / cond(full_b)
    num('BTwoCondMin', c_ratio.min())
    num('BTwoCondMax', c_ratio.max())
    # closed-orbit sensitivity |(M - I)^-1| ~ 1 / (2 sin(pi Q)) = 1 /
    # sqrt(det(M - I))
    num('BTwoSensMin', np.sqrt(c_ratio.min()))
    num('BTwoSensMax', np.sqrt(c_ratio.max()))
    sm = (
        b2[(b2.case == 'B-smooth') & (b2.steps_per_turn == 1440)]
        .set_index('metric')
        .value
    )
    num('BRfQz', sm['qz'], 3)
    num('BRfQr', sm['qr'], 3)

    def spread(x):
        return (
            x[(x.metric == 'r_eo') & x.code.isin(['lions', 'opal'])]
            .groupby(['code', 'integrator', 'energy_mev'])
            .value.agg(np.ptp)
        )

    s_ratio = (
        (spread(full_b) / spread(grad))
        .groupby(level=[0, 1])
        .agg(['min', 'max'])
    )
    for (code, integ), key in (
        (('lions', 'boris'), 'Lions'),
        (('opal', 'LF2'), 'LF'),
        (('opal', 'RK4'), 'RK'),
    ):
        num(f'BTwoSpread{key}Min', s_ratio.loc[(code, integ), 'min'])
        num(f'BTwoSpread{key}Max', s_ratio.loc[(code, integ), 'max'])
    lq = (
        b[(b.code == 'lions') & (b.steps_per_turn == 1440)]
        .set_index('metric')
        .value
    )
    num('BLionsQr', lq['qr'])
    num('BLionsQz', lq['qz'], 3)
    num('BCondDet', 4 * np.sin(np.pi * lq['qr']) ** 2)
    num('BSens', 1 / (2 * np.sin(np.pi * lq['qr'])))
    # smooth approximation (Zaremba and Kleeven, eqs. 14, 15), N = 4, k
    # = 0, xi = 0, flutter eps^2 / 2
    n_sec, flutter = 4, 0.3**2 / 2
    num('BFlutter', flutter, 2)
    num('BSmoothQz', np.sqrt(n_sec**2 / (n_sec**2 - 1) * flutter), 3)
    num(
        'BSmoothQrFrac',
        np.sqrt(1 + 3 * n_sec**2 / ((n_sec**2 - 1) * (n_sec**2 - 4)) * flutter)
        - 1,
        2,
    )


def case_d():
    """Table d2 and numbers, from case_d.parquet."""
    d = table('case_d/case_d.parquet')
    v = d.set_index(['code', 'integrator', 'metric']).value
    rows = []
    for code, integ, label in (
        ('lions', 'boris', 'LIONS Boris'),
        ('opal', 'LF2', 'OPAL LF2'),
        ('opal', 'RK4', 'OPAL RK4'),
        ('rf_track', 'rkf45', r'RF-Track rkf45, $10^{-12}$'),
    ):
        rows.append(
            [
                label,
                *(
                    sci(v[code, integ, m])
                    for m in ('l_secular_rel', 'l_max_dev_rel', 'p_rel_change')
                ),
            ]
        )
    write_table(
        'd2',
        [
            r'Code, integrator',
            r'secular $\Delta L/|L_0|$',
            r'$\max|L-L_0|/|L_0|$',
            r'$\Delta|p|/|p|$',
        ],
        rows,
    )
    num(
        'DOscRatio',
        v['lions', 'boris', 'l_max_dev_rel']
        / v['opal', 'LF2', 'l_max_dev_rel'],
    )
    num('DTurns', v['lions', 'boris', 'turns_tracked'], 1)
    num('DRkFourPred', 1e4 * (2 * np.pi) ** 6 / (144 * 360**5), 2)


def case_e():
    """Figure case_e_fd, tables e_grad, e_cost: case_e.parquet."""
    d = table('case_e/case_e.parquet')
    g = d[d.param != 'validation']
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    for case, marker in (('A', 'o'), ('B0', 's'), ('B', '^')):
        s = (
            g[(g.case == case) & (g.steps == 3600)]
            .groupby('fd_step')
            .rel_diff.max()
        )
        ax.loglog(
            s.index, s.values, marker + '-', label=f'case {case}, 3600 steps'
        )
    ax.axhline(1e-6, color='k', ls='--', lw=0.7, label='registered criterion')
    ax.set_xlabel(
        'finite-difference step $h$ (relative for $s$; m for $r_0$; '
        '$\\beta\\gamma$ for $p_{r0}$)'
    )
    ax.set_ylabel('max |autograd/FD - 1|')
    ax.legend(
        fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2
    )
    save(fig, 'case_e_fd')
    reg = g[g.registered]
    rows = [
        [
            case,
            sci(reg[(reg.case == case) & (reg.steps == 360)].rel_diff.max()),
            sci(reg[(reg.case == case) & (reg.steps == 3600)].rel_diff.max()),
        ]
        for case in ('A', 'B0', 'B')
    ]
    write_table('e_grad', ['Case', '360 steps', '3600 steps'], rows)
    ex = g[
        (g.case == 'B')
        & (g.steps == 3600)
        & (g.output == 'r')
        & (g.param == 's')
        & g.registered
    ].iloc[0]
    num('EExampleAutograd', ex.autograd, 4)
    num('EExampleFD', ex.fd, 4)
    sweep = (
        g[(g.case == 'B') & (g.steps == 3600)]
        .groupby('fd_step')
        .rel_diff.max()
    )
    num('ESweepLargeH', sweep[1e-4])
    num('ESweepSmallH', sweep[1e-8])
    num('ESweepTinyH', sweep[1e-9])
    fail = reg[(reg.case == 'A') & (reg.rel_diff > 1e-6)]
    num('EAFailCount', len(fail), 1)
    cost1 = (
        table('case_e/case_e_cost.parquet')
        .query('particles == 1 and not checkpoint')
        .set_index('steps')
        .peak_rss_mb
    )
    num(
        'EKBPerStep', (cost1[10000] - cost1[1000]) * 2**20 / 1e3 / 9000, 2
    )  # peak RSS in MiB
    large = (fail.output == 'phase') & (
        (fail.param == 'pr0') | ((fail.param == 'r0') & (fail.steps == 3600))
    )
    num('EAFailOtherMax', fail[~large].rel_diff.max(), 1)
    val = d[d.param == 'validation']
    num('EPortMax', val.rel_diff.max(), 1)
    c = table('case_e/case_e_cost.parquet')
    rows = []
    for (steps, particles), s in c.groupby(['steps', 'particles']):
        cells = []
        for ck in (False, True):
            x = s[s.checkpoint == ck].iloc[0]
            cells.append(
                'skipped'
                if isinstance(x.error, str)
                else f'{x.wall_s:.2g} s, {x.peak_rss_mb / 1000:.2g} GB'
            )
        rows.append([rf'\num{{{steps:.0e}}}', str(particles), *cells])
    both = (
        c[c.error.isna()]
        .pivot_table(
            index=['steps', 'particles'], columns='checkpoint', values='wall_s'
        )
        .dropna()
    )
    num('ECkptTimeMin', (both[True] / both[False]).min())
    num('ECkptTimeMax', (both[True] / both[False]).max())
    write_table(
        'e_cost',
        [
            'Steps',
            'Particles',
            'no checkpoint',
            r'checkpoint every $\sqrt{N}$',
        ],
        rows,
    )


E2_SITES = {
    'M3 CPU': 'case_e2_cpu/case_e2.parquet',
    'Apple GPU (MPS)': 'case_e2_mps/case_e2.parquet',
    'A100 (Jean-Zay)': 'sites/jeanzay/case_e2/'
    '1791641_jean-zay-iam02/case_e2.parquet',
    'H100 NVL (CC-IN2P3)': 'sites/cc/case_e2/'
    '62923301_ccwgslurm0200/case_e2.parquet',
}


def case_e2():
    """Figure case_e2 and table e2, from the E2_SITES tables."""
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    rows = []
    for (label, path), color in zip(
        E2_SITES.items(), ('#1f5f99', '#7f7f7f', '#c8553d', '#2e8b57')
    ):
        d = table(path)
        fwd = d[d['mode'] == 'forward']
        base = (
            fwd[fwd.code == 'lions']
            .set_index('particles')
            .particle_steps_per_s
        )
        torch = fwd[fwd.code == 'torch']
        best64 = (
            torch[torch.dtype == 'float64']
            .groupby('particles')
            .particle_steps_per_s.max()
        )
        best32 = torch[
            (torch.dtype == 'float32') & (torch.compiled == True)  # noqa: E712 (column)
        ].set_index('particles')
        if len(best64):
            ax.loglog(
                best64.index,
                best64.values,
                'o-',
                color=color,
                label=f'{label}, best float64',
            )
        else:
            ax.loglog(
                best32.index,
                best32.particle_steps_per_s,
                'o:',
                color=color,
                label=f'{label}, float32',
            )
        ax.loglog(base.index, base.values, 'x--', color=color, lw=0.7)
        fb = d[
            (d['mode'] == 'forward+backward')
            & (d.compiled == True)  # noqa: E712 (column)
            & (d.dtype == 'float64')
        ]

        def cell(n):
            return (
                f'{best64[n] / base[n]:.3g}'
                if n in best64.index
                else 'no float64'
            )

        rows.append(
            [
                label,
                cell(1000),
                cell(100000),
                f'{best32.particle_steps_per_s[100000] / base[100000]:.3g}'
                if 100000 in best32.index
                else '',
                f'{fb.wall_s.min():.3g} s' if len(fb) else 'fails',
            ]
        )
        key = ''.join(ch for ch in label.split()[0] if ch.isalpha())
        if len(best64):
            num(f'EtwoBest{key}', best64[100000], 3)
            num(f'EtwoRatio{key}', best64[100000] / base[100000], 3)
    a100 = table(E2_SITES['A100 (Jean-Zay)'])
    w = (
        a100[
            (a100['mode'] == 'forward')
            & (a100.code == 'torch')
            & (a100.dtype == 'float64')
            & (a100.compiled == True)  # noqa: E712 (column)
        ]
        .set_index('particles')
        .wall_s
    )
    num('EtwoAWallOne', w[1], 2)
    num('EtwoAWallMax', w[100000], 2)
    gpu = pd.concat(
        [
            table(E2_SITES[k])
            for k in ('A100 (Jean-Zay)', 'H100 NVL (CC-IN2P3)')
        ]
    )
    num(
        'EtwoGpuFloatThirtyTwoR',
        gpu[gpu.metric == 'cuda_float32_r_max_rel'].value.max(),
    )
    num(
        'EtwoGpuFloatSixtyFourR',
        gpu[
            gpu.metric.isin(
                ['cuda_float64_r_max_rel', 'cuda_float64_compiled_r_max_rel']
            )
        ].value.max(),
        1,
    )
    ax.set_xlabel('particles')
    ax.set_ylabel('particle-steps per second')
    ax.legend(fontsize=6, loc='upper left')
    ax.text(
        0.98,
        0.03,
        'dashed: plain LIONS on the same node',
        transform=ax.transAxes,
        ha='right',
        fontsize=6.5,
    )
    save(fig, 'case_e2')
    write_table(
        'e2',
        [
            'Hardware',
            r'float64, $10^3$',
            r'float64, $10^5$',
            r'float32 compiled, $10^5$',
            r'fwd+bwd, $10^3$',
        ],
        rows,
    )


def g1():
    """Figure g1 and its numbers, from the case G1 tables."""
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 4.8), sharex=True)
    fig.subplots_adjust(wspace=0.35)
    for col, (name, s2, s3) in enumerate(
        (
            (
                '3-sector',
                'case_g1/step2/case_g1_step2.parquet',
                'case_g1/step3/case_g1_step3.parquet',
            ),
            (
                '10-sector',
                'case_g1/10sector/step2/case_g1_step2.parquet',
                'case_g1/10sector/step3/case_g1_step3.parquet',
            ),
        )
    ):
        d2, d3 = table(s2), table(s3)
        d2 = d2[(d2['map'] == '1/3 deg')]
        d3 = d3[(d3['map'] == '1/3 deg')]
        for row, metric in enumerate(('qr_error', 'det_error')):
            ax = axes[row, col]
            s = (
                d2[d2.metric == metric]
                .assign(v=lambda x: x.value.abs())
                .groupby('steps_per_turn')
                .v.agg(['min', 'max'])
            )
            ax.fill_between(
                s.index, s['min'], s['max'], color=COLORS['lions'], alpha=0.25
            )
            ax.loglog(
                s.index,
                s['max'],
                'o-',
                color=COLORS['lions'],
                label='LIONS, bilinear',
            )
            for interp, color in (
                ('bicubic', COLORS['spline']),
                ('quintic', COLORS['rf_track']),
            ):
                t = d3[
                    (d3.interpolation == interp) & (d3.metric == metric)
                ].assign(v=lambda x: x.value.abs())
                t = t.groupby('steps_per_turn').v.agg(['min', 'max'])
                ax.fill_between(
                    t.index, t['min'], t['max'], color=color, alpha=0.2
                )
                ax.loglog(
                    t.index,
                    t['max'],
                    's-',
                    color=color,
                    label=f'Boris, {interp} spline',
                )
            ax.set_title(f'{name}, 1/3 deg', fontsize=8)
            ax.set_ylabel(
                r'$|\Delta Q_r|$'
                if metric == 'qr_error'
                else r'$|\det M - 1|$'
            )
        for ax in axes[1]:
            ax.set_xlabel('steps per turn')
    axes[0, 0].legend(fontsize=6)
    for ax in axes.flat:
        log_ticks(ax, [360, 720, 1440])
    save(fig, 'g1')
    for name, s2, s3 in (
        (
            'Three',
            'case_g1/step2/case_g1_step2.parquet',
            'case_g1/step3/case_g1_step3.parquet',
        ),
        (
            'Ten',
            'case_g1/10sector/step2/case_g1_step2.parquet',
            'case_g1/10sector/step3/case_g1_step3.parquet',
        ),
    ):
        d2, d3 = table(s2), table(s3)
        for tag, d in (
            (
                'Bil',
                d2[(d2['map'] == '1/3 deg') & (d2.steps_per_turn == 1440)],
            ),
            (
                'Spl',
                d3[
                    (d3['map'] == '1/3 deg')
                    & (d3.steps_per_turn == 1440)
                    & d3.interpolation.isin(['bicubic', 'quintic'])
                ],
            ),
        ):
            for metric, key in (
                ('qr_error', 'Qr'),
                ('det_error', 'Det'),
                ('r_eo_rel_error', 'Reo'),
            ):
                num(
                    f'G{name}{tag}{key}Max',
                    d[d.metric == metric].value.abs().max(),
                )
        spl = d3[
            (d3['map'] == '1/3 deg')
            & d3.interpolation.isin(['bicubic', 'quintic'])
        ]
        r = (
            spl[spl.metric == 'r_eo_rel_error']
            .assign(v=lambda x: x.value.abs())
            .pivot_table(
                index=['interpolation', 'a'],
                columns='steps_per_turn',
                values='v',
            )
        )
        ratios = pd.concat([r[360] / r[720], r[720] / r[1440]])
        num(f'G{name}SplOrbitRatioMin', ratios.min(), 3)
        num(f'G{name}SplOrbitRatioMax', ratios.max(), 3)
        tt = spl[
            (spl.metric == 'turn_time_rel_error')
            & (spl.steps_per_turn == 1440)
        ].value.abs()
        num(f'G{name}SplTurnTimeMax', tt.max())
        q360 = (
            d2[
                (d2['map'] == '1/3 deg')
                & (d2.steps_per_turn == 360)
                & (d2.metric == 'qr_error')
            ]
            .value.abs()
            .max()
        )
        num(f'G{name}BilQrMaxThreeSixty', q360)
        f4 = d2[d2['map'] == '4 deg']
        raw(
            f'G{name}FourDegFailed',
            f'{len(f4[f4.metric == "failed"])} of '
            f'{len(f4.groupby(["a", "steps_per_turn"]))}',
        )
    val = pd.concat(
        [
            table(p)[lambda x: x.metric == 'validation_vs_lions']
            for p in (
                'case_g1/step3/case_g1_step3.parquet',
                'case_g1/10sector/step3/case_g1_step3.parquet',
            )
        ]
    )
    num('GValidationMax', val.value.abs().max())
    t = table('case_g1/step1/case_g1_step1.parquet')
    t10 = table('case_g1/10sector/step1/case_g1_step1.parquet')
    num(
        'GTenRefQr',
        t10[t10.metric.isin(['qr_error', 'qz_error'])].value.abs().max(),
    )
    num('GTenRefOrbit', t10[t10.metric == 'r_eo_rel_error'].value.abs().max())
    num('GThreeRefQr', t[t.metric == 'qr_error'].value.abs().max())
    num('GThreeRefQz', t[t.metric == 'qz_error'].value.abs().max())
    num('GThreeRefOrbit', t[t.metric == 'r_eo_rel_error'].value.abs().max())


def g2():
    """Figure g2_tunes, table g2, numbers: the case G2 tables."""
    s2 = table('case_g2/step2/case_g2_step2.parquet')
    s3 = table('case_g2/step3/case_g2_step3.parquet')
    s3 = s3[
        ~((s3.code == 'opal') & (s3.energy_mev.round() == 512))
    ]  # OPAL returned a non-solution there
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 4.4), sharex=True)
    for code, d, style in (
        ('lions', s2, '-'),
        ('spline', s2, '-'),
        ('rftrack', s3, 'o'),
        ('opal', s3, 's'),
    ):
        for ax, metric in zip(axes, ('qr_error', 'qz_error')):
            s = d[(d.code == code) & (d.metric == metric)].sort_values(
                'energy_mev'
            )
            label = {
                'lions': 'LIONS (bilinear, Boris)',
                'spline': 'Boris on quintic spline',
                'rftrack': 'RF-Track rkf45 on spline',
                'opal': 'OPAL closed-orbit finder',
            }[code]
            ax.semilogy(
                s.energy_mev,
                s.value.abs(),
                style,
                color=COLORS[code],
                lw=0.8,
                ms=3,
                label=label,
            )
    for ax in axes:
        ax.axvspan(*psi_band(), color='0.85', zorder=0)
    axes[0].set_ylabel(r'$|Q_r - Q_r^{\rm FIXPO}|$')
    axes[1].set_ylabel(r'$|Q_z - Q_z^{\rm FIXPO}|$')
    axes[1].set_xlabel('kinetic energy [MeV]')
    axes[0].legend(fontsize=6, ncol=2)
    save(fig, 'g2_tunes')

    def interior(d):
        return d[
            (d.energy_mev > 80)
            & (d.energy_mev < 520)
            & ~d.qr_half_integer_band
        ]

    rows = []
    e3 = sorted(s3[s3.code == 'opal'].energy_mev.unique())
    for code, d, label in (
        ('lions', s2, 'LIONS (bilinear, Boris)'),
        ('spline', s2, 'Boris, quintic spline'),
        ('rftrack', s3, 'RF-Track, quintic spline'),
        ('opal', s3, 'OPAL'),
    ):
        d = interior(d[d.code == code])
        d = d[
            d.energy_mev.round(3).isin(np.round(e3, 3))
        ]  # the step 3 energies, without 512 MeV, for every code
        cells = []
        for metric in ('r_eo_rel_error', 'qr_error', 'qz_error'):
            v = d[d.metric == metric].value.abs()
            cells.append(f'{sci(v.median())} / {sci(v.max())}')
        rows.append([label, *cells])
        key = ''.join(ch for ch in code if ch.isalpha()).capitalize()
        num(f'GPsi{key}QrMed', d[d.metric == 'qr_error'].value.abs().median())
    s1 = table('case_g2/run2/step1/case_g2_step1.parquet')
    med = s1[s1.metric == 'dbz_dtheta_diff_median'].set_index('case').value
    num('GTwoDerivRatio', med['per_deg'] / med['per_rad'], 2)
    pr = s1[s1.run == '1b']
    good = pr[
        (pr.metric == 'pr_rel_error_pr_1e-3_m0c') & (pr.case == 'direction+1')
    ].value.abs()
    bad = pr[
        (pr.metric == 'pr_rel_error_pr_mrad') & (pr.case == 'direction+1')
    ].value.abs()
    num('GTwoPrGoodMax', good.max())
    num('GTwoPrMradMin', 100 * bad.min())
    num('GTwoPrMradMax', 100 * bad.max(), 3)
    write_table('g2', ['Code', r'$r_{eo}$ rel.', r'$Q_r$', r'$Q_z$'], rows)

    def q(d, code):
        return interior(d[d.code == code]).pipe(
            lambda x: (
                x[
                    (x.metric == 'qr_error')
                    & x.energy_mev.round(3).isin(np.round(e3, 3))
                ]
                .value.abs()
                .median()
            )
        )

    num('GPsiLionsOverOpal', q(s2, 'lions') / q(s3, 'opal'))
    lions = interior(s2[s2.code == 'lions'])
    num(
        'GPsiLionsDetMed',
        lions[lions.metric == 'det_error'].value.abs().median(),
    )
    num(
        'GPsiLionsDetMax',
        s2[(s2.code == 'lions') & (s2.metric == 'det_error')]
        .value.abs()
        .max(),
    )
    b = table('case_g2/step2b/case_g2_step2b.parquet')
    c = table('case_g2/step2c/case_g2_step2c.parquet')
    rows = []
    for e in sorted(b.energy_mev.unique()):

        def cell(d, code):
            return sci(
                d[
                    (d.code == code)
                    & (d.energy_mev == e)
                    & (d.metric == 'det_r_block')
                    & (d.jacobian_step_m == 1e-6)
                ].value.iloc[0]
            )  # stored as det - 1

        rows.append(
            [
                f'{e:.0f}',
                cell(b, 'lions'),
                cell(b, 'bilinear_numpy'),
                cell(c, 'rftrack_bilinear'),
                cell(b, 'spline'),
            ]
        )

    def det(d, code):
        return d[(d.code == code) & (d.metric == 'det_r_block')]

    lb = det(b, 'lions')
    var = lb.groupby('energy_mev').value.agg(
        lambda v: np.ptp(v) / v.abs().max()
    )
    num('GTwoBJacVar', 100 * var.max(), 1)
    ratio = (
        lb[lb.jacobian_step_m == 1e-6].set_index('energy_mev').value.abs()
        / det(c, 'rftrack_bilinear').set_index('energy_mev').value.abs()
    )
    num('GTwoCRfMax', det(c, 'rftrack_bilinear').value.abs().max(), 1)
    num('GTwoCRatioMin', ratio.min())
    num('GTwoCRatioMax', ratio.max(), 1)
    q2 = (
        pd.concat([b[b.jacobian_step_m == 1e-6], c])
        .query("metric == 'qr'")
        .pivot_table(index='energy_mev', columns='code', values='value')
    )
    num('GTwoCQrRfMax', (q2.rftrack_bilinear - q2.spline).abs().max())
    num('GTwoCQrLionsMin', (q2.lions - q2.spline).abs().min())
    num('GTwoCQrLionsMax', (q2.lions - q2.spline).abs().max())
    write_table(
        'g2_det',
        [
            '$E$ [MeV]',
            'LIONS',
            'numpy Boris, bilinear',
            'RF-Track, bilinear',
            'Boris, quintic spline',
        ],
        rows,
    )


H_MAPS = {
    '3-sector': 'case_h/case_h_matrix_3-sector.parquet',
    '10-sector': 'case_h/case_h_matrix_10-sector.parquet',
    'psi': 'case_h/case_h_matrix_psi.parquet',
}
H_FIELDS = {
    'F0_bilinear': 'F0 bilinear',
    'F1_hermite': 'F1 Hermite',
    'F2_cubic': 'F2 cubic spline',
    'F3_quintic': 'F3 quintic spline',
}


def h_worst(d):
    """max over cases of max(|Qr|, |Qz|, |det - 1|) per (spacing, field,
    integrator, steps); PSI band excluded.
    """
    energy = d.case.str.extract(r'E=([\d.]+)')[0].astype(float)
    d = d[~energy.between(*psi_band())]
    d = d[d.metric.isin(['qr_error', 'qz_error', 'det_error'])].assign(
        v=lambda x: x.value.abs()
    )
    return d.groupby(
        ['spacing', 'field', 'integrator', 'steps_per_turn']
    ).v.max()


def case_h():
    """Figure h_matrix, h_* tables, numbers: the H_MAPS tables."""
    worst = {m: h_worst(table(p)) for m, p in H_MAPS.items()}
    panels = (
        ('3-sector', '1/3 deg'),
        ('10-sector', '1/3 deg'),
        ('10-sector', '1 deg'),
        ('psi', '1/3 deg'),
    )
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.6), sharey=True)
    colors = dict(zip(H_FIELDS, ('#c8553d', '#1f5f99', '#e6a100', '#2e8b57')))
    for ax, (m, sp) in zip(axes, panels):
        w = worst[m]
        for field in H_FIELDS:
            for integ, style in (('boris', 'o-'), ('rk4', 's--')):
                s = w[sp, field, integ]
                ax.loglog(
                    s.index,
                    s.values,
                    style,
                    color=colors[field],
                    lw=0.9,
                    ms=3,
                    label=f'{H_FIELDS[field]}, '
                    f'{"Boris" if integ == "boris" else "RK4"}',
                )
        ax.axhline(1e-4, color='k', lw=0.7, ls=':')
        ax.set_title(f'{"PSI Ring" if m == "psi" else m}, {sp}', fontsize=8)
        log_ticks(ax, [720, 1440, 2880])
        ax.set_xlabel('steps per turn')
    axes[0].set_ylabel(r'worst of $|\Delta Q_r|, |\Delta Q_z|, |\det M-1|$')
    fig.legend(
        *axes[0].get_legend_handles_labels(),
        loc='lower center',
        ncol=4,
        fontsize=6.5,
        bbox_to_anchor=(0.5, -0.22),
    )
    save(fig, 'h_matrix')
    # PSI against FIXPO: matrix errors are against RF-Track (q - q_rft);
    # add RF-Track minus FIXPO from G2 step 3
    psi = table(H_MAPS['psi'])
    psi = psi.assign(
        energy=psi.case.str.extract(r'E=([\d.]+)')[0].astype(float)
    )
    s3 = table('case_g2/step3/case_g2_step3.parquet')
    rft = s3[s3.code == 'rftrack'].set_index(['energy_mev', 'metric']).value
    opal_e = s3[
        (s3.code == 'opal') & (s3.energy_mev.round() != 512)
    ].energy_mev.unique()
    keep = (
        psi.energy.isin(opal_e)
        & (psi.energy > 80)
        & (psi.energy < 520)
        & ~psi.energy.between(*psi_band())
    )
    for field, key in (
        ('F0_bilinear', 'Bil'),
        ('F1_hermite', 'Herm'),
        ('F2_cubic', 'Cub'),
        ('F3_quintic', 'Quin'),
    ):
        for metric, mk in (('qr_error', 'Qr'), ('qz_error', 'Qz')):
            x = psi[
                keep
                & (psi.field == field)
                & (psi.integrator == 'boris')
                & (psi.steps_per_turn == 1440)
                & (psi.metric == metric)
            ]
            err = (x.value + [rft[e, metric] for e in x.energy]).abs()
            num(f'HFixpo{key}{mk}Med', err.median())
    # det(M) - 1 alone, worst case per option, 1440 steps
    det = {
        m: table(p)
        .pipe(
            lambda d: d[
                ~d.case.str.extract(r'E=([\d.]+)')[0]
                .astype(float)
                .between(*psi_band())
            ]
        )
        .pipe(
            lambda d: d[(d.metric == 'det_error') & (d.steps_per_turn == 1440)]
        )
        .assign(v=lambda d: d.value.abs())
        .groupby(['spacing', 'field', 'integrator'])
        .v.max()
        for m, p in H_MAPS.items()
    }
    f0f1 = [
        det[m][sp, 'F0_bilinear', 'boris'] / det[m][sp, 'F1_hermite', 'boris']
        for m in det
        for sp in det[m].index.get_level_values(0).unique()
    ]
    num('HDetRatioMin', min(f0f1))
    num('HDetRatioMax', max(f0f1))
    br = [
        det[m][sp, 'F0_bilinear', 'boris'] / det[m][sp, 'F0_bilinear', 'rk4']
        for m in det
        for sp in det[m].index.get_level_values(0).unique()
    ]
    num('HBorisOverRkMin', min(br))
    num('HBorisOverRkMax', max(br))
    rk = [
        det[m][sp, 'F0_bilinear', 'rk4']
        for m in ('10-sector', 'psi')
        for sp in det[m].index.get_level_values(0).unique()
    ]
    num('HFzeroRkDetMin', min(rk))
    num('HFzeroRkDetMax', max(rk))
    num('HFoneDetPsi', det['psi']['1/3 deg', 'F1_hermite', 'boris'])
    psi_det = table(H_MAPS['psi']).pipe(
        lambda d: d[
            (d.field == 'F1_hermite')
            & (d.integrator == 'boris')
            & (d.steps_per_turn == 2880)
            & (d.metric == 'det_error')
        ]
    )
    psi_det = psi_det[
        ~psi_det.case.str.extract(r'E=([\d.]+)')[0]
        .astype(float)
        .between(*psi_band())
    ]
    num('HFoneDetPsiTwoEightEighty', psi_det.value.abs().max())
    cost = table('case_h/case_h_cost.parquet')
    c = (
        cost[(cost.metric == 's_per_particle_step') & (cost.particles == 1e5)]
        .set_index('code')
        .value
    )
    base = c['F0_bilinear_boris']
    rows = []
    for field, label in H_FIELDS.items():
        rows.append(
            [
                label,
                f'{sci(c[field + "_boris"], 3)} '
                f'({c[field + "_boris"] / base:.2g})',
                f'{sci(c[field + "_rk4"], 3)} '
                f'({c[field + "_rk4"] / base:.2g})',
            ]
        )
        num(
            'HCost'
            + {'F0': 'Bil', 'F1': 'Herm', 'F2': 'Cub', 'F3': 'Quin'}[
                field[:2]
            ],
            c[field + '_boris'] / base,
        )
    write_table(
        'h_cost', ['Field model', 'Boris [s] (ratio)', 'RK4 [s] (ratio)'], rows
    )
    num('HLionsCost', c['lions'] / base)
    rows = []
    for m, sp in (
        ('3-sector', '1/3 deg'),
        ('3-sector', '1 deg'),
        ('10-sector', '1/3 deg'),
        ('10-sector', '1 deg'),
        ('psi', '1/3 deg'),
    ):
        w = worst[m]
        cells = []
        for field in H_FIELDS:
            for integ in ('boris', 'rk4'):
                s = w[sp, field, integ]
                ok = s[s <= 1e-4]
                cells.append(str(int(ok.index.min())) if len(ok) else '--')
        rows.append([f'{"PSI" if m == "psi" else m}, {sp}', *cells])
    write_table(
        'h_decision',
        [
            'Map',
            *(f'{f.split("_")[0]} {i}' for f in H_FIELDS for i in ('B', 'RK')),
        ],
        rows,
    )
    v = table('case_h/case_h_validate.parquet')
    num('HFOneNodes', v[v.check == 'F1_nodes'].value.abs().max(), 2)
    p2 = table('case_h_rerun/case_h_phase2.parquet')
    ot = p2[p2.metric == 'one_turn_vs_harness'].value.abs()
    num('HPhaseTwoOneTurn', ot.max())
    dif = p2[p2.metric.isin(['qr_error', 'det_error'])]
    num('HPhaseTwoJacMax', (dif.value - dif.harness).abs().max())
    num('HTargetRatio', 1e-4 / (dif.value - dif.harness).abs().max(), 1)
    pc = p2[p2.metric == 's_per_particle_step'].set_index('case').value
    num(
        'HPhaseTwoCostRatio',
        pc['lions_hermite_100000'] / pc['lions_bilinear_100000'],
    )
    az = table('case_h/case_h_azimuth.parquet')
    rows = []
    for m, label in (('10-sector', '10-sector, 1 deg'), ('psi', 'PSI Ring')):
        a = az[(az['map'] == m) & (az.metric == 'max_abs_dM_rr_vs_numpy')]
        cells = []
        for azm in ('arccos', 'atan2'):
            s = a[a.azimuth == azm].set_index('jacobian_step_m').value
            cells.append(' / '.join(sci(s[h], 2) for h in (1e-4, 1e-6, 1e-7)))
        rows.append([label, *cells])
    write_table(
        'h_azimuth',
        ['Map', r'arccos, $h=10^{-4}/10^{-6}/10^{-7}$ m', r'atan2, same $h$'],
        rows,
    )


def main():
    """Write every figure, table and the .tex number files in GEN."""
    GEN.mkdir(exist_ok=True)
    for f in (
        boris_polygon,
        azimuth_error,
        interpolation_orders,
        field_maps,
        interpolation_kinks,
        case_a,
        case_b,
        case_d,
        case_e,
        case_e2,
        g1,
        g2,
        case_h,
    ):
        f()
    lines = [
        rf'\newcommand{{\{k}}}{{{v}}}' for k, v in sorted(NUMBERS.items())
    ]
    (GEN / 'numbers.tex').write_text('\n'.join(lines) + '\n')
    src = [
        r'\begin{tabular}{ll}',
        r'\toprule',
        r'Table & sha256 (first 16) \\',
        r'\midrule',
        *(
            rf'\texttt{{{k.replace("_", chr(92) + "_")}}} & '
            rf'\texttt{{{v[:16]}}} \\'
            for k, v in sorted(SOURCES.items())
        ),
        r'\bottomrule',
        r'\end{tabular}',
    ]
    (GEN / 'sources.tex').write_text('\n'.join(src) + '\n')
    (GEN / 'code.tex').write_text(
        rf'\newcommand{{\CodeVersion}}{{\texttt{{{CODE["code_version"]}}}}}'
        + '\n'
    )
    for k, v in sorted(NUMBERS.items()):
        print(f'{k:28s} {v}')


if __name__ == '__main__':
    main()
