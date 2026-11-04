"""Trainer: T0 closed forms against LIONS, then every stage played by a
script: it can be won, and it can be lost.
"""

import math

import numpy as np
import pytest

pytest.importorskip('yaml')

import synthetic as syn  # noqa: E402 (after importorskip)
from lions.game import physics  # noqa: E402 (after importorskip)
from lions.game.stages import Action, load_scenario  # noqa: E402 (after importorskip)

DT = 1 / 60
W0 = 10.0  # MeV
STAGES = {s.key: s for s in load_scenario()}


def track(turns, npaspt=360):
    """x, y, z of one proton, one row per step."""
    rows = []

    def record(d):
        rows.append((d.bunch.x[0], d.bunch.y[0], d.bunch.z[0]))

    syn.quiet_run('deck.dat', on_step=record)
    return np.array(rows)


def test_uniform_field_against_lions(workdir):
    b = -1.2
    rho = physics.radius(W0, b)
    syn.write_bmap('bmap.dat', b, rmin=0.05, rmax=1.0, nr=200)
    syn.write_deck(
        'deck.dat',
        fhf=physics.revolution_frequency(W0, b),
        r0=rho,
        w0=W0,
        turns=5,
    )
    xyz = track(5)
    assert np.hypot(xyz[-1, 0], xyz[-1, 1]) == pytest.approx(rho, rel=1e-4)
    assert xyz[-1, 0] == pytest.approx(
        rho, rel=1e-4
    )  # back at the start after whole turns: the frequency is right
    assert math.copysign(1, xyz[10, 1]) == physics.rotation(b)


def test_smooth_qz_against_lions(workdir):
    """4 sectors, no mean gradient, flutter eps^2 / 2: the tracked Qz
    within 3% of the smooth approximation.
    """
    b, eps = -1.0, 0.3
    syn.write_bmap(
        'bmap.dat',
        bz=lambda r, t: b * (1 + eps * np.cos(4 * t)),
        dbzx=lambda r, t: 4 * eps * b * np.sin(t) * np.sin(4 * t) / r,
        dbzy=lambda r, t: -4 * eps * b * np.cos(t) * np.sin(4 * t) / r,
        ntet=721,
        rmin=0.05,
        rmax=0.8,
        nr=151,
    )
    syn.write_deck(
        'deck.dat',
        fhf=physics.revolution_frequency(W0, b),
        r0=physics.radius(W0, b),
        w0=W0,
        z0=1e-4,
        turns=64,
        npaspt=720,
    )
    z = track(64)[719::720, 2]
    window = np.hanning(len(z))
    spectrum = np.abs(np.fft.rfft(z * window, 16 * len(z)))
    tracked = np.argmax(spectrum) / (16 * len(z))
    assert tracked == pytest.approx(
        math.sqrt(physics.qz_squared(0.0, 4, eps**2 / 2, 0.0)), rel=0.03
    )


def play(stage, policy, seed=1, limit=200.0):
    stage.reset(seed)
    t = 0.0
    while not stage.done and t < limit:
        stage.step(policy(stage), DT)
        t += DT
    assert stage.done, f'{stage.key} did not end within {limit} s'
    return stage


def larmor_bot(s):
    if s.flight is not None:
        return Action()
    need = (
        s.target_rotation
        * physics.momentum(s.w)
        * physics.PROTON_KG
        * physics.C
        / (physics.E_CHARGE * s.target)
    )
    b = -need  # rotation +1 needs Bz < 0
    if abs(b) > s.p['b_max']:
        return Action(x=-1)
    if abs(b - s.b) > s.p['b_rate'] * DT:
        return Action(y=1 if b > s.b else -1)
    return Action(fire=not s.last_action.fire)


def dee_bot(s):
    f_c = physics.revolution_frequency(s.w, s.p['b'])
    wanted = f_c * (
        1 - s.phi / math.pi / 4
    )  # slide the phase to the crest, then hold
    return Action(
        x=np.clip((wanted - s.f) / (s.p['f_rate'] * s.f_nominal * DT), -1, 1)
    )


def trim_bot(s):
    wanted = s.needed_delta() + s.phi / 40
    return Action(
        y=np.clip((wanted - s.delta) / (s.p['delta_rate'] * DT), -1, 1)
    )


def fallout_bot(s):
    s_k = math.sqrt(physics.isochronous_index(s.w))
    if s.last_action.y:
        return Action()
    unstable = s.z + s.dz / s_k  # the growing mode of z'' = k z
    if abs(unstable) > 1.0:
        return Action(y=-math.copysign(1, unstable))
    return Action()


def thomas_bot(s):
    wanted = 0.75**2
    return Action(y=np.clip((wanted - s.qz2) * 20, -1, 1))


@pytest.mark.parametrize(
    'key, bot',
    [
        ('larmor', larmor_bot),
        ('dee_beat', dee_bot),
        ('trim_coils', trim_bot),
        ('fallout', fallout_bot),
        ('thomas', thomas_bot),
    ],
)
@pytest.mark.parametrize('seed', range(5))
def test_stage_can_be_won(key, bot, seed):
    stage = play(STAGES[key], bot, seed)
    assert stage.success, (
        f'{key} seed {seed}: {stage.message} {stage.metrics()}'
    )
    assert stage.points > 0


@pytest.mark.parametrize(
    'key', ['dee_beat', 'trim_coils', 'fallout', 'thomas']
)
def test_stage_can_be_lost(key):
    stage = play(STAGES[key], lambda s: Action(), 1)
    assert not stage.success, stage.message


def test_score_records_the_budget_and_no_player():
    stage = play(STAGES['larmor'], larmor_bot, 3)
    record = stage.score()
    assert record['seed'] == 3 and record['budget']['steps'] == stage.steps
    assert {'scenario', 'stage', 'switches', 'budget', 'metrics'} <= set(
        record
    )
    assert not {'player', 'name', 'initials'} & set(record)


def test_same_seed_same_run():
    assert (
        play(STAGES['dee_beat'], dee_bot, 7).score()
        == play(STAGES['dee_beat'], dee_bot, 7).score()
    )


def test_unit_mismatch_is_an_error(monkeypatch):
    from lions.game import stages

    monkeypatch.setitem(stages.STAGES['larmor'].UNITS, 'b_max', 'kG')
    with pytest.raises(ValueError, match='b_max'):
        load_scenario()


def test_logged_units():
    from lions.game.runlog import units_of

    assert (
        units_of('dz_mm_rad'),
        units_of('z_mm'),
        units_of('phi_rad'),
        units_of('w_MeV'),
        units_of('qz2'),
    ) == ('mm/rad', 'mm', 'rad', 'MeV', '1')


def test_arcade_full_flow(tmp_path):
    """A headless session through the real widget: consent, one stage
    logged, game over, high-score table.
    """
    import os

    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest
    from PyQt5.QtWidgets import QApplication
    import pyarrow.parquet as pq
    from lions.game.arcade import Arcade, LIVES
    from lions.game.runlog import Store

    app = QApplication.instance() or QApplication([])
    window = Arcade(Store(tmp_path))
    window.timer.stop()
    window.show()
    QTest.keyClick(window, Qt.Key_Return)
    assert window.state == 'consent'
    QTest.keyClick(window, Qt.Key_Y)
    assert window.state == 'menu' and window.store.consent is True
    for _ in range(LIVES):  # lose the first stage until the lives run out
        while window.state != 'play':
            QTest.keyClick(window, Qt.Key_Return)
        window.stage.points = 100
        window.stage.finish(False, 'TEST')
        window.frame()
        window.grab()
    assert window.state == 'over'
    QTest.keyClick(window, Qt.Key_Return)
    assert window.state == 'initials'
    for _ in range(3):
        QTest.keyClick(window, Qt.Key_Space)
    assert window.state == 'attract'
    assert window.store.scores()[0]['initials'] == 'AAA'
    session = window.store.session
    QTest.keyClick(window, Qt.Key_Return)
    assert window.state == 'consent'  # the next player is asked again
    QTest.keyClick(window, Qt.Key_N)
    assert window.store.consent is False and window.store.session != session
    assert window.store.log_run(window.stages[0], [{'step': 1}]) is None
    logs = list((tmp_path / 'runs').glob('*.parquet'))
    assert len(logs) == LIVES
    meta = pq.read_schema(logs[0]).metadata
    assert b'score' in meta and b'AAA' not in b''.join(meta.values())
    app.processEvents()


# ---------- missions

MISSIONS = {s.mission_id: s for s in load_scenario('missions')}


def tuned(mission, seed=1, **changes):
    """A mission reset, with the knobs at their best settings."""
    mission.reset(seed)
    mission.knobs = {**mission.ideal(), **changes}
    mission.refresh()
    return mission


def test_extraction_energy_follows_the_field():
    """The energy at the extraction radius is set by B rho = p / q,
    whatever the path.
    """
    from lions.game import machine as mc

    s = tuned(MISSIONS['NEON'])
    assert s.beam.w_mean / 20 == pytest.approx(10.0, rel=0.01)
    assert mc.Ion('Ne', 20, 4).rigidity(200.0) == pytest.approx(
        s.ideal()['b0']
        * (1 + s.errors.field)
        * mc.Ion('Ne', 20, 4).gamma(200.0)
        / mc.Ion('Ne', 20, 4).gamma(0.4)
        * s.machine.radius,
        rel=1e-9,
    )


def test_turns_match_the_energy_gain():
    from lions.game import machine as mc

    s = tuned(MISSIONS['F18'])
    per_turn = mc.energy_gain(
        s.machine, s.ion, s.knobs['dee_kv'], s.harmonic, 0.0
    )
    assert 18.0 / per_turn < s.beam.turns < 3 * 18.0 / per_turn + 51


@pytest.mark.parametrize(
    'knob, factor', [('rf_mhz', 1.01), ('b0', 1.01), ('trim', 0.8)]
)
def test_off_tune_loses_the_beam(knob, factor):
    s = MISSIONS['EYE']
    good = tuned(s).beam.transmission
    ideal = s.ideal()
    assert (
        tuned(s, **{knob: ideal[knob] * factor}).beam.transmission < 0.5 * good
    )


def test_the_rf_locks_on_two_harmonics():
    """Harmonic 2 and 3 both accelerate, with their own transit factors;
    between them nothing does.
    """
    s = MISSIONS['NEON']
    rf3 = tuned(s).ideal()['rf_mhz']
    h3, h2 = tuned(s).beam, tuned(s, rf_mhz=rf3 * 2 / 3).beam
    assert (h3.harmonic, h2.harmonic) == (3, 2)
    assert (
        h3.transmission > 0.3 and h2.transmission > 0.3 and h2.turns > h3.turns
    )
    assert tuned(s, rf_mhz=rf3 * 5 / 6).beam.transmission == 0


def test_buncher_increases_transmission():
    s = MISSIONS['F18']
    assert (
        tuned(s, buncher=1.0).beam.transmission
        > 1.3 * tuned(s, buncher=0.0).beam.transmission
    )


def test_higher_dee_voltage_separates_turns():
    s = MISSIONS['F18']
    assert (
        tuned(s, dee_kv=60.0).beam.lost_kw < tuned(s, dee_kv=20.0).beam.lost_kw
    )


@pytest.mark.parametrize('w, cm', [(70.0, 4.0), (230.0, 33.0)])
def test_proton_range(w, cm):
    assert physics.proton_range(w) == pytest.approx(cm, rel=0.05)


def test_bragg_peak_at_the_residual_range():
    s = tuned(MISSIONS['DEEP'])
    assert s.depth_axis[int(np.argmax(s.bragg))] == pytest.approx(
        s.p['depth_cm'], abs=0.8
    )


def test_f18_saturates_at_the_sourced_yield(monkeypatch):
    s = tuned(MISSIONS['F18'])
    s.beam.w_mean, s.beam.extracted_ua, s.beam.x, s.beam.y = (
        18.0,
        10.0,
        0.0,
        0.0,
    )
    monkeypatch.setitem(s.p, 'activity_gbq', 1e9)
    for _ in range(5000):
        s.deliver(1.0)
    assert s.activity == pytest.approx(80.0, rel=1e-3)  # 8 GBq/uA x 10 uA


def deliver_all(s, dt=0.05):
    """Beam on at the tuned setting, off at the prescription for a
    therapy mission.
    """
    t = 0.0
    while not s.done and t < s.p['time_limit']:
        ready = hasattr(s, 'healthy') and s.delivered >= s.p['dose_gy']
        s.step(Action(fire=not ready), dt)
        s.knobs = {**s.ideal(), **s.fixes}
        t += dt
    return s


FIXES = {
    'EYE': {},
    'DEEP': {},
    'QUAL': {'quad': 0.5},
    'F18': {'dee_kv': 60.0},
    'NEON': {'slit': 0.0015, 'source_ua': 1.0},
    'ISOL': {'quad': 0.25},
}


@pytest.mark.parametrize('name', sorted(MISSIONS))
@pytest.mark.parametrize('seed', range(2))
def test_mission_can_be_won_when_tuned(name, seed):
    s = tuned(MISSIONS[name], seed)
    s.fixes = FIXES[name]
    s.knobs = {**s.ideal(), **s.fixes}
    s.refresh()
    deliver_all(s)
    assert s.success, f'{name}: {s.message} {s.metrics()}'


@pytest.mark.parametrize('name', sorted(MISSIONS))
def test_untuned_start_delivers_little(name):
    s = MISSIONS[name]
    s.reset(1)
    t = 0.0
    while not s.done and t < 5.0:
        s.step(Action(fire=True), DT)
        t += DT
    assert not s.success, s.message


# ---------- coach and environment


def follow_coach(s, seed):
    from lions.game.coach import Coach

    s.reset(seed)
    s.tutorial = True
    coach = Coach(s)
    keys = [x['key'] for x in s.specs]
    while not s.done and s.time < 500:
        a = coach.advise()
        x, y = 0, 0
        if a.knob:
            i = keys.index(a.knob)
            if i != s.selected:
                y = (
                    0
                    if s.last_action.y
                    else (1 if (s.selected - i) % len(keys) == 1 else -1)
                )
            else:
                x = a.direction * (0.1 if a.fine else 1)
        s.step(Action(x, y, a.beam), 0.1)
    s.tutorial = False
    return s


@pytest.mark.parametrize('name', sorted(MISSIONS))
def test_following_the_coach_finishes_the_mission(name):
    s = follow_coach(MISSIONS[name], 11)
    assert s.success, f'{name}: {s.message} {s.metrics()}'


def test_coach_reads_no_hidden_state():
    import inspect
    from lions.game import coach

    source = inspect.getsource(coach)
    assert '.errors' not in source and 'ideal(' not in source


def test_no_beam_reads_no_position_and_scores_on_the_probe_only():
    """Gaming check: losing the beam cannot read as a clean, centred
    beam."""
    s = MISSIONS['EYE']
    s.reset(1)
    s.knobs['rf_mhz'] = 104.0
    s.refresh()
    assert (
        not s.has_beam() and math.isnan(s.readings['x_mm']) and not s.ready()
    )
    assert s.tune_score() <= 0.3
    assert tuned(s, 1).tune_score() > 0.3


@pytest.mark.parametrize('name', sorted(MISSIONS))
def test_parking_on_the_knob_edges_scores_low(name):
    """Gaming check: every knob at a limit of its range gains
    nothing."""
    s = MISSIONS[name]
    s.reset(1)
    for side in ('min', 'max'):
        for spec in s.specs:
            s.knobs[spec['key']] = spec[side]
        s.refresh()
        assert not s.ready() and s.tune_score() <= 0.5


def test_saturating_a_limit_does_not_read_as_success():
    """Gaming check: more current than the target takes is out of spec,
    not a high score.
    """
    s = MISSIONS['F18']
    assert tuned(s, 1, source_ua=400.0).ready() is False


def test_tuned_machines_meet_the_spec():
    for name, m in MISSIONS.items():
        tuned(m, 1, **FIXES[name])
        assert m.ready(), f'{name}: {m.spec_terms()}'


def test_environment_api():
    pytest.importorskip('stable_baselines3')
    from stable_baselines3.common.env_checker import check_env
    from lions.game.env import MissionEnv

    env = MissionEnv('F18', seeds=[1, 2, 3])
    check_env(env, warn=True)
    obs, info = env.reset(seed=0)
    assert env.observation_space.contains(obs) and info['condition'] in (
        1,
        2,
        3,
    )
    assert not any('error' in str(k) for k in vars(env))


@pytest.mark.parametrize('name', sorted(MISSIONS))
@pytest.mark.parametrize('tutorial', [False, True])
def test_every_mission_draws(name, tutorial):
    import os

    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    import tempfile
    from PyQt5.QtWidgets import QApplication
    from lions.game.arcade import Arcade
    from lions.game.runlog import Store

    app = QApplication.instance() or QApplication([])
    w = Arcade(Store(tempfile.mkdtemp()))
    w.timer.stop()
    w.new_game()
    w.stages, w.index = [m for m in w.missions if m.mission_id == name], 0
    w.state = 'brief'
    w.grab()
    w.start_stage(tutorial=tutorial)
    w.grab()  # no beam at the start
    w.stage.knobs = w.stage.ideal()
    w.stage.refresh()
    w.grab()  # beam
    assert (w.coach is not None) == tutorial
    app.processEvents()


@pytest.mark.parametrize('key', sorted(STAGES))
def test_following_the_helper_finishes_each_stage(key):
    from lions.game.coach import coach_for

    stage = STAGES[key]
    stage.tutorial = True
    stage.reset(5)
    coach = coach_for(stage)
    while not stage.done and stage.time < 300:
        a = coach.advise()
        stage.step(Action(a.x, a.y, a.fire), DT)
    stage.tutorial = False
    assert stage.success, f'{key}: {stage.message} {stage.metrics()}'


def test_every_knob_has_help():
    from lions.game.coach import KNOB_HELP

    assert {s['key'] for m in MISSIONS.values() for s in m.specs} <= set(
        KNOB_HELP
    )
    assert (
        max(len(text) for text in KNOB_HELP.values()) <= 80
    )  # one line of the helper box


def test_tutorial_for_each_challenge_with_real_keys(tmp_path):
    """Down to each line, T: that challenge's own walkthrough, then its
    helper. A failed tutorial costs no life.
    """
    import os

    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest
    from PyQt5.QtWidgets import QApplication
    from lions.game.arcade import Arcade
    from lions.game.runlog import Store

    app = QApplication.instance() or QApplication([])
    w = Arcade(Store(tmp_path))
    w.timer.stop()
    w.show()
    w.new_game()
    items = w.menu_items()
    for line, (key, title, stages) in enumerate(items):
        w.state, w.choice = 'menu', 0
        for _ in range(line):
            QTest.keyClick(w, Qt.Key_Down)
        assert w.choice == line
        w.grab()
        QTest.keyClick(w, Qt.Key_T)
        if line:
            assert w.state == 'walkthrough' and w.stages[0].mission_id == key
            w.grab()
            QTest.keyClick(w, Qt.Key_Return)
        assert w.state == 'play' and w.tutorial and w.stage is stages[0]
        w.frame()
        w.grab()
        lives = w.lives
        w.stage.finish(False, 'TEST')
        w.frame()
        assert w.state == 'clear' and w.lives == lives and w.total == 0
    w.state, w.choice = 'menu', 1
    QTest.keyClick(w, Qt.Key_A)  # A and ENTER accept the highlighted line
    assert w.state == 'brief' and w.stages[0].mission_id == items[1][0]
    app.processEvents()


@pytest.mark.parametrize('name', sorted(MISSIONS))
def test_walkthrough_covers_every_knob_and_the_helper(name):
    """Each challenge has its own steps, naming every knob it offers;
    the helper's advice always maps to one.
    """
    from lions.game.coach import Coach

    s = MISSIONS[name]
    named = [k for step in s.walkthrough for k in step['knobs']]
    missing = {spec['key'] for spec in s.specs} - set(named)
    assert not missing, f'{name}: {missing}'
    s.reset(3)
    advice = Coach(s).advise()
    assert advice.knob in named or not advice.knob


TUNING_MISSION = 'F18'


def fresh_mission(seed=3):
    from lions.game.env import mission_by_id

    m = mission_by_id(TUNING_MISSION)
    m.reset(seed)
    return m


def test_a_pass_drifts_the_field():
    from lions.game import tuners

    m = fresh_mission()
    before = m.errors.field
    run = tuners.Run(m)
    run.m.knobs.update(m.knobs)
    try:
        run.passed()
    except tuners.Run.Done:
        pass
    assert (
        m.errors.field != before
        and abs(m.errors.field - before) < 10 * m.p['drift']
    )


def test_the_environment_drifts_too():
    pytest.importorskip('gymnasium')
    from lions.game.env import MissionEnv

    env = MissionEnv(TUNING_MISSION, seeds=[3])
    env.reset(seed=0)
    before = env.m.errors.field
    env.step(np.zeros(env.action_space.shape))
    assert env.m.errors.field != before


@pytest.mark.parametrize(
    'name',
    [
        'COACH',
        'NELDER-MEAD',
        'NSGA-II (MOGA)',
        'CMA-ES',
        'RANDOM',
        'BAYESIAN OPT.',
    ],
)
def test_every_method_keeps_the_budget_and_repeats(name):
    from lions.game import tuners

    method = next(x for x in tuners.METHODS if x.name == name)
    if method.missing():
        pytest.skip(f'needs {method.missing()}')
    budget = 25
    first, second = (
        tuners.tune(fresh_mission(), method.fn, 7, budget) for _ in range(2)
    )
    assert (
        first.evals <= budget
        and first.evals == second.evals
        and first.best == second.best
    )
    assert first.reached is None or first.reached == first.evals


def test_nsga_objectives_rank_no_beam_behind_any_beam(monkeypatch):
    from lions.game import tuners

    m = fresh_mission()
    run = tuners.Run(m)
    m.knobs['rf_mhz'] = next(s for s in m.specs if s['key'] == 'rf_mhz')['min']
    monkeypatch.setattr(run, 'set', lambda u: m.refresh())
    terms = run.terms(run.start)
    assert (
        not m.has_beam()
        and min(terms) > 1e9
        and len(terms) == len(m.spec_terms())
    )


def test_a_policy_runs_on_the_same_budget(tmp_path):
    pytest.importorskip('stable_baselines3')
    from lions.game import tuners

    path = tmp_path / 'F18_ppo_0.zip'
    tuners.train_policy(TUNING_MISSION, 'ppo', 0, path, steps=64)
    run = tuners.tune(
        fresh_mission(),
        tuners.policy,
        1,
        budget=5,
        model=tuners.load_policy(path),
    )
    assert run.evals <= 5


def arcade_window(tmp_path, monkeypatch):
    import os

    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtWidgets import QApplication
    from lions.game import tuners
    from lions.game.arcade import Arcade
    from lions.game.runlog import Store

    monkeypatch.setattr(tuners, 'POLICIES', tmp_path / 'policies')
    app = QApplication.instance() or QApplication([])
    w = Arcade(Store(tmp_path))
    w.timer.stop()
    w.show()
    w.new_game()
    return app, w


@pytest.mark.parametrize(
    'name', ['COACH', 'NELDER-MEAD', 'NSGA-II (MOGA)', 'CMA-ES', 'RANDOM']
)
def test_each_challenge_runs_a_method_with_real_keys(
    tmp_path, monkeypatch, name
):
    """Down to the line, O, choose, ENTER: passes on the screen, in the
    end the real machine agrees. No points, no life, no log.
    """
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest
    from lions.game.tuners import METHODS

    app, w = arcade_window(tmp_path, monkeypatch)
    QTest.keyClick(w, Qt.Key_Down)
    QTest.keyClick(w, Qt.Key_O)
    assert w.state == 'pick'
    w.grab()
    for _ in range([m.name for m in METHODS].index(name)):
        QTest.keyClick(w, Qt.Key_Down)
    assert METHODS[w.method_line].name == name
    seed = w.machine_seed
    QTest.keyClick(w, Qt.Key_Return)
    assert w.state == 'tuning'
    w.grab()
    w.tuning.thread.join(timeout=120)
    while not w.tuning_over:
        w.frame()
    w.grab()
    run = w.tuning.run
    assert w.tuning.error is None and w.replayed == run.evals
    assert w.tuning.m.measured_ready() == (run.reached is not None)
    assert (
        w.stage.errors.steer_x == w.tuning.m.errors.steer_x
    )  # the screen's machine is the method's machine
    lives, total = w.lives, w.total
    QTest.keyClick(w, Qt.Key_O)
    assert w.state == 'pick' and w.machine_seed == seed and name in w.tuned
    w.grab()
    assert (w.lives, w.total) == (lives, total) and not list(
        tmp_path.rglob('*.parquet')
    )
    QTest.keyClick(w, Qt.Key_Escape)
    assert w.state == 'menu'
    w.close()
    app.processEvents()


def test_take_the_beam_after_a_method(tmp_path, monkeypatch):
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest

    app, w = arcade_window(tmp_path, monkeypatch)
    QTest.keyClick(w, Qt.Key_Down)
    QTest.keyClick(w, Qt.Key_O)
    QTest.keyClick(w, Qt.Key_Return)
    w.tuning.thread.join(timeout=60)
    while not w.tuning_over:
        w.frame()
    QTest.keyClick(w, Qt.Key_Return)
    assert w.state == 'play' and w.tutorial and w.coach and w.tuning is None
    w.frame()
    w.grab()
    w.close()


def test_level_0_has_no_tuning_method(tmp_path, monkeypatch):
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest

    _, w = arcade_window(tmp_path, monkeypatch)
    QTest.keyClick(w, Qt.Key_O)
    assert w.state == 'menu'
    w.close()


def test_closing_during_a_method_ends_its_thread(tmp_path, monkeypatch):
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest

    _, w = arcade_window(tmp_path, monkeypatch)
    QTest.keyClick(w, Qt.Key_Down)
    QTest.keyClick(w, Qt.Key_O)
    w.method_line = 2  # Bayesian optimisation: long enough to be running
    QTest.keyClick(w, Qt.Key_Return)
    thread = w.tuning.thread
    w.close()
    thread.join(timeout=10)
    assert not thread.is_alive()


def test_behind_window_opens_on_b_and_says_what_runs(tmp_path, monkeypatch):
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest
    from lions.game.behind import CODES

    _, w = arcade_window(tmp_path, monkeypatch)
    QTest.keyClick(w, Qt.Key_B)
    assert w.behind.isVisible()
    assert 'does not run LIONS or OPAL' in w.behind.text.toPlainText()
    QTest.keyClick(w, Qt.Key_Down)
    QTest.keyClick(w, Qt.Key_Return)
    w.behind.refresh()
    assert w.stages[0].title in w.behind.text.toPlainText()
    QTest.keyClick(w, Qt.Key_Return)
    w.behind.refresh()
    assert 'lions.game.machine' in w.behind.live.text()
    assert {'LIONS (tracker)', 'OPAL, OPALX'} <= {
        c[0] for c in CODES if c[1] == 'NOT RUNNING'
    }
    w.grab()
    w.close()
    assert not w.behind.isVisible()


def test_the_game_imports_no_tracking_code():
    """The Behind window says LIONS and OPAL do not run. If a level
    starts to call them, this fails: update the window.
    """
    import ast
    from pathlib import Path
    import lions.game

    allowed = ('lions.game', 'lions.constants')
    for path in Path(lions.game.__file__).parent.glob('*.py'):
        for node in ast.walk(ast.parse(path.read_text())):
            names = (
                [n.name for n in node.names]
                if isinstance(node, ast.Import)
                else [node.module]
                if isinstance(node, ast.ImportFrom)
                else []
            )
            for name in names:
                assert not name.startswith('lions') or name.startswith(
                    allowed
                ), f'{path.name} imports {name}'
                assert not any(
                    code in name.lower()
                    for code in ('opal', 'rf_track', 'rftrack')
                ), f'{path.name} imports {name}'


def test_train_saves_a_policy_with_held_out_result_with_real_keys(
    tmp_path, monkeypatch
):
    pytest.importorskip('stable_baselines3')
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest
    from lions.game import train, tuners

    app, w = arcade_window(tmp_path, monkeypatch)
    train.main(
        ['F18'], steps=2048
    )  # past the first progress write, on a machine with no policy folder yet
    report = tuners.policy_report('F18')
    assert (
        report['conditions'] == len(tuners.HELD_OUT)
        and report['steps'] == 2048
        and report['seed'] == 0
    )
    assert not tuners.policy_path('F18').with_suffix('.progress').exists()
    w.choice = [key for key, _, _ in w.menu_items()].index('F18')
    QTest.keyClick(w, Qt.Key_O)
    w.method_line = [m.fn for m in tuners.METHODS].index(tuners.policy)
    w.grab()
    QTest.keyClick(w, Qt.Key_Return)
    assert w.state == 'tuning'
    w.tuning.thread.join(timeout=60)
    assert w.tuning.error is None and w.tuning.run.evals >= 1
    w.close()


def test_training_from_the_game_is_a_subprocess_that_esc_stops(
    tmp_path, monkeypatch
):
    pytest.importorskip('stable_baselines3')
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest
    from lions.game import tuners

    monkeypatch.setenv('HOME', str(tmp_path))
    _, w = arcade_window(tmp_path, monkeypatch)
    QTest.keyClick(w, Qt.Key_Down)
    QTest.keyClick(w, Qt.Key_O)
    w.method_line = [m.fn for m in tuners.METHODS].index(tuners.policy)
    QTest.keyClick(w, Qt.Key_Return)
    assert w.state == 'training'
    process = w.training
    w.frame()
    w.grab()
    assert process.poll() is None
    QTest.keyClick(w, Qt.Key_Escape)
    assert (
        w.state == 'pick' and w.training is None and process.poll() is not None
    )
    w.close()


def test_reduced_model_returns_a_beam_from_the_tables_alone():
    import dataclasses
    from pathlib import Path

    import numpy as np
    import pytest

    from lions.game import machine as mc
    from lions.game.env import mission_by_id
    from lions.game.reduced import Reduced

    table = (
        Path(mc.__file__).parent / 'scenarios' / 'f18_synthetic_orbits.parquet'
    )
    model = Reduced(table)
    m = mission_by_id('F18')
    m.reset(10000)
    m.knobs = m.ideal()
    t0 = mc.run(m.machine, m.ion, m.knobs, m.errors, m.draws)
    beam = model.run(m.machine, m.ion, m.knobs, m.errors, m.draws)
    assert isinstance(beam, mc.Beam) and {
        f.name for f in dataclasses.fields(beam)
    } == {f.name for f in dataclasses.fields(t0)}
    assert 0 < len(beam.energies) <= mc.MACROPARTICLES
    assert abs(beam.energies.mean() / m.p['w_target'] - 1) < 0.02
    # Same knobs, same result, and T0 is untouched by it.
    again = model.run(m.machine, m.ion, m.knobs, m.errors, m.draws)
    assert np.array_equal(beam.energies, again.energies)
    assert np.array_equal(
        t0.energies,
        mc.run(m.machine, m.ion, m.knobs, m.errors, m.draws).energies,
    )
    with pytest.raises(ValueError, match='harmonic'):
        model.run(
            m.machine,
            m.ion,
            {**m.knobs, 'rf_mhz': m.knobs['rf_mhz'] * 1.5},
            m.errors,
            m.draws,
        )


@pytest.mark.parametrize(
    'mission_id', ['EYE', 'DEEP', 'F18', 'NEON', 'ISOL', 'QUAL']
)
def test_measured_terms_are_the_true_terms_without_noise(mission_id):
    from lions.game.env import mission_by_id

    m = mission_by_id(mission_id)
    m.reset(5)
    m.knobs = m.ideal()
    m.p['noise'] = 0.0
    m.refresh()
    assert m.has_beam()
    assert m.measured_terms() == pytest.approx(
        m.spec_terms(), rel=1e-12, abs=1e-12
    )
    assert m.measured_ready() == m.ready()


@pytest.mark.parametrize(
    'mission_id', ['EYE', 'DEEP', 'F18', 'NEON', 'ISOL', 'QUAL']
)
def test_the_measured_terms_read_the_noise_and_the_true_terms_do_not(
    mission_id,
):
    from lions.game.env import mission_by_id

    m = mission_by_id(mission_id)
    m.reset(5)
    m.knobs = m.ideal()
    m.p['noise'] = 0.05
    true = set()
    measured = {k: set() for k in m.spec_terms()}
    for _ in range(12):
        m.refresh()
        true.add(tuple(m.spec_terms().values()))
        for k, v in m.measured_terms().items():
            measured[k].add(v)
    assert len(true) == 1 and all(
        len(v) > 1 for k, v in measured.items() if k != 'spread'
    )


def test_run_stops_on_the_measured_test_unless_told_otherwise():
    from lions.game import tuners

    calls = []

    def stop():
        calls.append(1)
        return True

    run = tuners.tune(fresh_mission(), tuners.random_search, 1, stop=stop)
    assert run.reached == 1 and len(calls) == 1
    assert tuners.Run(fresh_mission()).stop.__name__ == 'measured_ready'


def test_a_wider_box_moves_only_its_knob():
    from lions.game import tuners

    m = fresh_mission()
    plain = tuners.Run(m)
    key, spec = plain.keys[0], m.specs[0]
    below, above = (
        2 * tuners.BOX_ERRORS * spec['start_error'],
        3 * tuners.BOX_ERRORS * spec['start_error'],
    )
    wide = tuners.Run(m, box={key: (below, above)})
    assert wide.lo[0] == max(m.knobs[key] - below, spec['min']) and wide.hi[
        0
    ] == min(m.knobs[key] + above, spec['max'])
    assert list(wide.lo[1:]) == list(plain.lo[1:]) and list(
        wide.hi[1:]
    ) == list(plain.hi[1:])


def test_log_terms_rank_no_beam_behind_any_beam(monkeypatch):
    from lions.game import tuners

    m = fresh_mission()
    run = tuners.Run(m)
    m.knobs['rf_mhz'] = next(s for s in m.specs if s['key'] == 'rf_mhz')['min']
    monkeypatch.setattr(run, 'set', lambda u: m.refresh())
    terms = run.log_terms(run.start)
    assert (
        not m.has_beam()
        and min(terms) > math.log(tuners.TERM_CLAMP)
        and len(terms) == len(m.spec_terms())
    )


@pytest.mark.parametrize('name', ['turbo', 'constrained_bo', 'mobo'])
def test_botorch_methods_spend_the_budget(name):
    pytest.importorskip('botorch')
    from lions.game import tuners

    m = fresh_mission()
    budget = 2 * len(m.specs) + 22
    run = tuners.tune(m, getattr(tuners, name), 1, budget=budget)
    assert (run.evals == budget or run.reached) and run.model_failures == 0


def test_the_screen_says_when_only_the_readings_are_in_spec(
    tmp_path, monkeypatch
):
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest
    from lions.game.tuners import METHODS

    _, w = arcade_window(tmp_path, monkeypatch)
    QTest.keyClick(w, Qt.Key_Down)
    QTest.keyClick(w, Qt.Key_O)
    assert {'TURBO', 'CONSTRAINED BO', 'MOBO'} <= {m.name for m in METHODS}
    s, method = w.stages[w.index], METHODS[0]
    w.tuned[method.name] = (7, 0.9, False)
    assert (
        w.method_state(method, s) == 'READS IN SPEC, 7 PASSES. TRUE BEAM OUT'
    )
    w.tuned[method.name] = (7, 0.9, True)
    assert w.method_state(method, s) == 'IN SPEC, 7 PASSES'
    w.close()


def test_model_slot_holds_t0_in_game_and_table_in_benchmark(
    monkeypatch,
):
    from lions.game import machine as mc
    from lions.game.missions import Mission
    from lions.game.reduced import install_from_env

    m = MISSIONS['NEON']
    m.reset(1)
    assert Mission.run_model is mc.run
    monkeypatch.setattr(Mission, 'run_model', Mission.run_model)
    monkeypatch.setenv('CASE_REDUCED', 'neon_synthetic_orbits.parquet:10')
    assert install_from_env()
    m.knobs |= m.ideal()
    m.refresh()
    assert m.beam.transmission > 0.5
    m.knobs['rf_mhz'] *= 2 / 3
    m.refresh()
    assert not m.has_beam()
