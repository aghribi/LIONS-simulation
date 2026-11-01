"""Trainer: T0 closed forms against LIONS, then every stage played by a script: it can be won, and it can be lost."""
import math

import numpy as np
import pytest

pytest.importorskip('yaml')

import lions  # noqa: E402
import synthetic as syn  # noqa: E402
from lions.game import physics  # noqa: E402
from lions.game.stages import Action, load_scenario  # noqa: E402

DT = 1 / 60
W0 = 10.0    # MeV
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
    syn.write_deck('deck.dat', fhf=physics.revolution_frequency(W0, b), r0=rho, w0=W0, turns=5)
    xyz = track(5)
    assert np.hypot(xyz[-1, 0], xyz[-1, 1]) == pytest.approx(rho, rel=1e-4)
    assert xyz[-1, 0] == pytest.approx(rho, rel=1e-4)      # back at the start after whole turns: the frequency is right
    assert math.copysign(1, xyz[10, 1]) == physics.rotation(b)


def test_smooth_qz_against_lions(workdir):
    """4 sectors, no mean gradient, flutter eps^2 / 2: the tracked Qz within 3% of the smooth approximation."""
    b, eps = -1.0, 0.3
    syn.write_bmap('bmap.dat', bz=lambda r, t: b * (1 + eps * np.cos(4 * t)),
                   dbzx=lambda r, t: 4 * eps * b * np.sin(t) * np.sin(4 * t) / r,
                   dbzy=lambda r, t: -4 * eps * b * np.cos(t) * np.sin(4 * t) / r, ntet=721, rmin=0.05, rmax=0.8, nr=151)
    syn.write_deck('deck.dat', fhf=physics.revolution_frequency(W0, b), r0=physics.radius(W0, b), w0=W0, z0=1e-4,
                   turns=64, npaspt=720)
    z = track(64)[719::720, 2]
    window = np.hanning(len(z))
    spectrum = np.abs(np.fft.rfft(z * window, 16 * len(z)))
    tracked = np.argmax(spectrum) / (16 * len(z))
    assert tracked == pytest.approx(math.sqrt(physics.qz_squared(0.0, 4, eps**2 / 2, 0.0)), rel=0.03)


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
    need = s.target_rotation * physics.momentum(s.w) * physics.PROTON_KG * physics.C / (physics.E_CHARGE * s.target)
    b = -need     # rotation +1 needs Bz < 0
    if abs(b) > s.p['b_max']:
        return Action(x=-1)
    if abs(b - s.b) > s.p['b_rate'] * DT:
        return Action(y=1 if b > s.b else -1)
    return Action(fire=not s.last_action.fire)


def dee_bot(s):
    f_c = physics.revolution_frequency(s.w, s.p['b'])
    wanted = f_c * (1 - s.phi / math.pi / 4)    # slide the phase to the crest, then hold
    return Action(x=np.clip((wanted - s.f) / (s.p['f_rate'] * s.f_nominal * DT), -1, 1))


def trim_bot(s):
    wanted = s.needed_delta() + s.phi / 40
    return Action(y=np.clip((wanted - s.delta) / (s.p['delta_rate'] * DT), -1, 1))


def fallout_bot(s):
    s_k = math.sqrt(physics.isochronous_index(s.w))
    if s.last_action.y:
        return Action()
    unstable = s.z + s.dz / s_k       # the growing mode of z'' = k z
    if abs(unstable) > 1.0:
        return Action(y=-math.copysign(1, unstable))
    return Action()


def thomas_bot(s):
    wanted = 0.75**2
    return Action(y=np.clip((wanted - s.qz2) * 20, -1, 1))


@pytest.mark.parametrize('key, bot', [('larmor', larmor_bot), ('dee_beat', dee_bot), ('trim_coils', trim_bot),
                                      ('fallout', fallout_bot), ('thomas', thomas_bot)])
@pytest.mark.parametrize('seed', range(5))
def test_stage_can_be_won(key, bot, seed):
    stage = play(STAGES[key], bot, seed)
    assert stage.success, f'{key} seed {seed}: {stage.message} {stage.metrics()}'
    assert stage.points > 0


@pytest.mark.parametrize('key', ['dee_beat', 'trim_coils', 'fallout', 'thomas'])
def test_stage_can_be_lost(key):
    stage = play(STAGES[key], lambda s: Action(), 1)
    assert not stage.success, stage.message


def test_score_records_the_budget_and_no_player():
    stage = play(STAGES['larmor'], larmor_bot, 3)
    record = stage.score()
    assert record['seed'] == 3 and record['budget']['steps'] == stage.steps
    assert {'scenario', 'stage', 'switches', 'budget', 'metrics'} <= set(record)
    assert not {'player', 'name', 'initials'} & set(record)


def test_same_seed_same_run():
    assert play(STAGES['dee_beat'], dee_bot, 7).score() == play(STAGES['dee_beat'], dee_bot, 7).score()


def test_unit_mismatch_is_an_error(monkeypatch):
    from lions.game import stages
    monkeypatch.setitem(stages.STAGES['larmor'].UNITS, 'b_max', 'kG')
    with pytest.raises(ValueError, match='b_max'):
        load_scenario()


def test_logged_units():
    from lions.game.runlog import units_of
    assert (units_of('dz_mm_rad'), units_of('z_mm'), units_of('phi_rad'), units_of('w_MeV'), units_of('qz2')) == \
        ('mm/rad', 'mm', 'rad', 'MeV', '1')


def test_arcade_full_flow(tmp_path):
    """A headless session through the real widget: consent, one stage logged, game over, high-score table."""
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
    assert window.state == 'brief' and window.store.consent is True
    for _ in range(LIVES):    # lose the first stage until the lives run out
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
    assert window.state == 'consent'      # the next player is asked again
    QTest.keyClick(window, Qt.Key_N)
    assert window.store.consent is False and window.store.session != session
    assert window.store.log_run(window.stages[0], [{'step': 1}]) is None
    logs = list((tmp_path / 'runs').glob('*.parquet'))
    assert len(logs) == LIVES
    meta = pq.read_schema(logs[0]).metadata
    assert b'score' in meta and b'AAA' not in b''.join(meta.values())
    app.processEvents()

