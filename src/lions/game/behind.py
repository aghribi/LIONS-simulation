"""What runs behind the game: a separate window. It states what the game
executes and what it does not. tests/test_game.py checks the imports.
"""

import time
from importlib import metadata

from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QLabel, QTextBrowser, QVBoxLayout, QWidget

from lions.game import machine as mc
from lions.game.missions import BUDGET, Mission
from lions.game.tuners import METHODS, POLICIES, policy_report

# name, in the game?, what it does here, how far it is checked
CODES = [
    (
        'lions.game.physics',
        'RUNNING',
        'Closed forms: rho = p / (q B), revolution frequency, smooth Qz and Qr, proton range, straggling, scattering. All of level 0 and parts of the missions.',
        'Radius, frequency, rotation and smooth Qz tracked by LIONS in tests/test_game.py. Range formula from Bortfeld 1997.',
    ),
    (
        'lions.game.machine',
        'RUNNING',
        'Mission machine, tier T0: 120 macroparticles, turn by turn, no field map.',
        'Not validated against LIONS. Tests check its own consistency only: extraction energy, turn count, losses off tune, harmonics.',
    ),
    (
        'LIONS (tracker)',
        'NOT RUNNING',
        'Used by two tests of the closed forms. No level of the game calls it. SPEC levels 1 to 5 would (tiers T1, T3): not built.',
        'Fortran bit-identity and benchmarks, outside the game.',
    ),
    (
        'OPAL, OPALX',
        'NOT RUNNING',
        'Offline reference codes for benchmarks of LIONS on field maps. Nothing from them enters the game.',
        'benchmarks/, benchmarks/DESIGN.md.',
    ),
    (
        'RF-Track, FROM_ORBIT',
        'NOT RUNNING',
        'Offline references, as above.',
        'benchmarks/.',
    ),
    (
        'PyTorch fast model (T2)',
        'NOT BUILT',
        'SPEC tier for gradients and RL on a map. PPO here trains on lions.game.machine instead.',
        'Needs the Hermite field first.',
    ),
]
OMITTED = 'No field map, no isochronism error from a real map, no misalignment of the dees, no hysteresis, no finite actuator speed beyond the knob rate, no space charge beyond one term in Qz, no extraction tracking.'
PASS = [
    'Knobs, set by you or by a method',
    'machine.run: for each turn, energy gain from the dee voltage and phase; phase slip from RF, field and gamma; Qz from flutter, spiral, gradient and current; losses on the dee aperture and when behind the RF',
    'Septum, energy slit, deflector and steerers: position and spot at the target, by closed forms',
    'Mission physics: Bragg peak, activity, flux, power density',
    'Readings with noise: what the control room shows',
    'measured_terms: one normalised error per client requirement, from the readings. tune_score: the same terms as one number. spec_terms and ready: the same on the true beam, not seen by any method',
]
PACKAGES = (
    'numpy',
    'scipy',
    'scikit-learn',
    'pymoo',
    'cma',
    'stable-baselines3',
    'torch',
    'PyQt5',
)


def version(package):
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return 'not installed'


def page(stage):
    """The explanation, for the challenge in front of the player."""
    rows = ''.join(
        f'<tr><td><b>{name}</b></td><td style="color:{"#3dff7a" if state == "RUNNING" else "#ff4d4d"}">{state}</td><td>{what}</td><td>{checked}</td></tr>'
        for name, state, what, checked in CODES
    )
    steps = ''.join(f'<li>{step}</li>' for step in PASS)
    methods = ''.join(
        f'<tr><td><b>{m.name}</b></td><td>{m.what}</td></tr>' for m in METHODS
    )
    packages = ', '.join(f'{p} {version(p)}' for p in PACKAGES)
    here = ''
    if isinstance(stage, Mission):
        report = policy_report(stage.mission_id)
        trained = (
            f"{report['in_spec']}/{report['conditions']} held-out conditions in spec, seed {report['seed']}, {report['steps']} training passes"
            if report
            else 'none trained yet'
        )
        switches = ', '.join(
            f'{key}: {value}' for key, value in stage.switches.items()
        )
        here = (
            f'<h3>This challenge: {stage.title}</h3><p>{stage.ion.name}, harmonic {stage.harmonic}, {len(stage.specs)} knobs. A training machine built to agree with the formulas, not a copy of a real one. '
            f'Parameters with no source: {switches}. Sources are listed at the top of scenarios/missions.yaml.</p><p>PPO policy: {trained}.</p>'
        )
    return f"""<h2>What runs behind ION RUSH</h2>
<p><b>The game runs closed-form models. It does not run LIONS or OPAL.</b> A score here says nothing about a real machine.</p>
<table border="1" cellpadding="5" cellspacing="0"><tr><th>Code</th><th>In the game</th><th>What it does</th><th>Checked by</th></tr>{rows}</table>
{here}
<h3>One pass of the machine</h3><ol>{steps}</ol>
<h3>What the model leaves out</h3><p>{OMITTED}</p>
<h3>Tuning methods (key O in a challenge)</h3><p>Each gets {BUDGET} passes, sees the readings or the score, never the hidden errors. They work on a copy of the machine and the screen replays each pass. No points, no lives, nothing logged. Comparison on held-out conditions: docs/training/RL_RUNS.md.</p>
<table border="1" cellpadding="5" cellspacing="0">{methods}</table>
<p>Trained policies: {POLICIES}</p>
<p>Installed: {packages}.</p>"""


class Behind(QWidget):
    """Reads the game, never writes to it. Hidden errors are not shown: it
    would give the answer away.
    """

    def __init__(self, arcade):
        super().__init__()
        self.arcade, self.shown = arcade, object()
        self.setWindowTitle('ION RUSH: what runs behind')
        self.resize(900, 720)
        self.live = QLabel()
        self.live.setWordWrap(True)
        self.text = QTextBrowser()
        layout = QVBoxLayout(self)
        layout.addWidget(self.live)
        layout.addWidget(self.text)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.refresh)
        self.timer.start(500)
        self.refresh()

    def refresh(self):
        stage = self.arcade.current_stage()
        if stage is not self.shown:
            self.shown = stage
            self.text.setHtml(page(stage))
        self.live.setText(f'<b>Now:</b> {self.status(stage)}')

    def status(self, stage):
        tuning = self.arcade.tuning
        if not isinstance(stage, Mission):
            return (
                f'{stage.title}: lions.game.physics, closed forms.'
                if stage
                else 'No stage running.'
            )
        if not hasattr(stage, 'beam'):
            return f'{stage.title}: not started. It will run on lions.game.machine. LIONS and OPAL: not running.'
        started = time.perf_counter()
        beam = mc.run(
            stage.machine, stage.ion, stage.knobs, stage.errors, stage.draws
        )
        driver = (
            f'{tuning.method.name}, {len(tuning.passes)} passes computed, {BUDGET} at most'
            if tuning
            else 'no tuning method running'
        )
        return f'{stage.title}: lions.game.machine, {beam.turns} turns, one pass takes {1e3 * (time.perf_counter() - started):.0f} ms. {driver}. LIONS and OPAL: not running.'
