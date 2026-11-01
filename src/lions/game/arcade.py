"""ION RUSH: the arcade front end of the trainer. Draws the stages of lions.game.stages; the physics is all there."""
import math
import random
import sys

from PyQt5.QtCore import QPointF, QRectF, Qt, QTimer
from PyQt5.QtGui import QColor, QFont, QPainter, QPainterPath, QPen, QPixmap
from PyQt5.QtWidgets import QApplication, QWidget

from lions.game import physics
from lions.game.runlog import Store
from lions.game.missions import Mission
from lions.game.stages import Action, load_scenario

W, H = 960, 720    # virtual screen, scaled to the window
DT = 1 / 60
LIVES = 3
CYAN, MAGENTA, YELLOW, GREEN, RED, WHITE = (QColor(c) for c in ('#27e8ff', '#ff3df2', '#ffe14d', '#3dff7a', '#ff4d4d', '#e8f0ff'))
DIM = QColor('#3a4a78')
BACK = QColor('#05060f')
FONT = 'Menlo'
CONTROLS = {'larmor': 'UP DOWN: FIELD B   LEFT RIGHT: ENERGY   SPACE: FIRE', 'dee_beat': 'LEFT RIGHT: RF FREQUENCY',
            'trim_coils': 'UP DOWN: FIELD TRIM', 'fallout': 'TAP UP DOWN: CORRECTOR KICK',
            'thomas': 'UP DOWN: FLUTTER F   LEFT RIGHT: SPIRAL XI'}
MISSION_CONTROLS = 'UP DOWN: PICK A KNOB   LEFT RIGHT: TURN IT   SHIFT: FINE   HOLD SPACE: BEAM ON'


def wrapped(text, width):
    lines, line = [], ''
    for word in text.split():
        if len(line) + len(word) + 1 > width:
            lines.append(line)
            line = word
        else:
            line = f'{line} {word}'.strip()
    return [*lines, line]


class Arcade(QWidget):
    def __init__(self, store=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle('ION RUSH')
        self.resize(W, H)
        self.setFocusPolicy(Qt.StrongFocus)
        self.store = store or Store()
        self.level0, self.missions = load_scenario(), load_scenario('missions')
        self.stages = self.level0
        self.keys, self.rows = set(), []
        self.clock, self.state = 0.0, 'attract'
        self.seed_source = random.Random()
        self.scan = self.scanlines()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.frame)
        self.timer.start(16)

    def scanlines(self):
        pixmap = QPixmap(W, H)
        pixmap.fill(Qt.transparent)
        p = QPainter(pixmap)
        p.setPen(QPen(QColor(0, 0, 0, 70), 1))
        for y in range(0, H, 3):
            p.drawLine(0, y, W, y)
        p.end()
        return pixmap

    # ---------- flow

    def frame(self):
        self.clock += DT
        if self.state == 'play':
            self.play_step()
        self.update()

    def new_game(self):
        self.lives, self.total, self.choice = LIVES, 0, 0
        self.completed = set()
        self.state = 'menu'

    def menu_items(self):
        return [('LEVEL 0', 'HOW A CYCLOTRON WORKS', self.level0), *((m.mission_id, m.title, [m]) for m in self.missions)]

    def choose(self):
        self.stages, self.index = self.menu_items()[self.choice][2], 0
        self.state = 'brief'

    def start_stage(self):
        self.stage = self.stages[self.index]
        self.stage.reset(self.seed_source.randrange(2**31))
        self.rows = []
        self.state = 'play'

    def play_step(self):
        ax = (Qt.Key_Right in self.keys) - (Qt.Key_Left in self.keys)
        ax = ax * 0.1 if Qt.Key_Shift in self.keys else ax
        ay = (Qt.Key_Up in self.keys) - (Qt.Key_Down in self.keys)
        action = Action(ax, ay, Qt.Key_Space in self.keys)
        observation, _, done = self.stage.step(action, DT)
        self.rows.append({'step': self.stage.steps, 'action_x': ax, 'action_y': ay, 'action_fire': action.fire, **observation})
        if done:
            self.store.log_run(self.stage, self.rows)
            self.total += self.stage.points
            if self.stage.success or not self.stage.mandatory:
                self.state = 'clear'
            else:
                self.lives -= 1
                self.state = 'lost' if self.lives else 'over'

    def press_start(self):
        state = self.state
        if state == 'attract':
            self.state = 'consent'
        elif state == 'brief':
            self.start_stage()
        elif state == 'clear':
            self.index += 1
            if self.index < len(self.stages):
                self.state = 'brief'
            else:
                self.completed.add(self.menu_items()[self.choice][0])
                self.state = 'won' if len(self.completed) == len(self.menu_items()) else 'menu'
        elif state == 'menu':
            self.choose()
        elif state == 'lost':
            self.state = 'brief'
        elif state in ('won', 'over'):
            if self.store.qualifies(self.total):
                self.initials, self.letter = ['A', 'A', 'A'], 0
                self.state = 'initials'
            else:
                self.state = 'attract'

    def keyPressEvent(self, event):
        key = event.key()
        if event.isAutoRepeat():
            return
        self.keys.add(key)
        if self.state == 'consent' and key in (Qt.Key_Y, Qt.Key_N):
            self.store.new_player(key == Qt.Key_Y)
            self.new_game()
        elif self.state == 'initials':
            self.edit_initials(key)
        elif self.state == 'menu' and key in (Qt.Key_Up, Qt.Key_Down):
            self.choice = (self.choice + (1 if key == Qt.Key_Down else -1)) % len(self.menu_items())
        elif key in (Qt.Key_Return, Qt.Key_Enter) or (key == Qt.Key_Space and self.state != 'play'):
            self.press_start()
        elif key == Qt.Key_Escape:
            self.state = 'over' if self.state == 'menu' else 'menu' if self.state in ('brief', 'play', 'lost', 'clear') else 'attract'

    def keyReleaseEvent(self, event):
        if not event.isAutoRepeat():
            self.keys.discard(event.key())

    def edit_initials(self, key):
        letter = ord(self.initials[self.letter]) - 65
        if key in (Qt.Key_Up, Qt.Key_Down):
            self.initials[self.letter] = chr(65 + (letter + (1 if key == Qt.Key_Up else -1)) % 26)
        elif key in (Qt.Key_Return, Qt.Key_Enter, Qt.Key_Space, Qt.Key_Right):
            self.letter += 1
            if self.letter == 3:
                self.store.add_score(''.join(self.initials), self.total)
                self.state = 'attract'

    # ---------- drawing tools

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.fillRect(self.rect(), BACK)
        scale = min(self.width() / W, self.height() / H)
        p.translate((self.width() - W * scale) / 2, (self.height() - H * scale) / 2)
        p.scale(scale, scale)
        p.setClipRect(QRectF(0, 0, W, H))
        self.p = p
        getattr(self, f'draw_{self.state}')()
        p.drawPixmap(0, 0, self.scan)
        p.end()

    def text(self, x, y, string, size=20, color=WHITE, align='center'):
        p = self.p
        font = QFont(FONT, size, QFont.Bold)
        font.setStyleHint(QFont.Monospace)
        p.setFont(font)
        flag = {'center': Qt.AlignHCenter, 'left': Qt.AlignLeft, 'right': Qt.AlignRight}[align]
        left = {'center': x - W, 'left': x, 'right': x - 2 * W}[align]
        rect = QRectF(left, y, 2 * W, size * 2)
        glow = QColor(color)
        glow.setAlpha(60)
        p.setPen(glow)
        for dx, dy in ((-1.5, 0), (1.5, 0), (0, -1.5), (0, 1.5)):
            p.drawText(rect.translated(dx, dy), flag, string)
        p.setPen(color)
        p.drawText(rect, flag, string)

    def block(self, y, string, width, size=16, color=WHITE):
        for line in wrapped(string, width):
            self.text(W / 2, y, line, size, color)
            y += size * 1.7
        return y

    def neon(self, draw, color, width=2):
        glow = QColor(color)
        glow.setAlpha(55)
        self.p.setBrush(Qt.NoBrush)
        self.p.setPen(QPen(glow, width * 4, Qt.SolidLine, Qt.RoundCap))
        draw()
        self.p.setPen(QPen(color, width, Qt.SolidLine, Qt.RoundCap))
        draw()

    def line(self, x1, y1, x2, y2, color, width=2):
        self.neon(lambda: self.p.drawLine(QPointF(x1, y1), QPointF(x2, y2)), color, width)

    def circle(self, x, y, r, color, width=2):
        self.neon(lambda: self.p.drawEllipse(QPointF(x, y), r, r), color, width)

    def ion(self, x, y, color=WHITE):
        glow = QColor(color)
        glow.setAlpha(70)
        self.p.setPen(Qt.NoPen)
        self.p.setBrush(glow)
        self.p.drawEllipse(QPointF(x, y), 11, 11)
        self.p.setBrush(color)
        self.p.drawEllipse(QPointF(x, y), 4.5, 4.5)
        self.p.setBrush(Qt.NoBrush)

    def bar(self, x, y, width, height, fraction, color, label=''):
        self.p.setPen(QPen(DIM, 2))
        self.p.setBrush(Qt.NoBrush)
        self.p.drawRect(QRectF(x, y, width, height))
        self.p.fillRect(QRectF(x + 2, y + 2, (width - 4) * min(max(fraction, 0), 1), height - 4), color)
        if label:
            self.text(x, y - 24, label, 13, color, 'left')

    def hud(self):
        s = self.stage
        self.text(16, 8, f'SCORE {self.total + s.points:07d}', 18, CYAN, 'left')
        self.text(W / 2, 8, f'{self.index + 1}/{len(self.stages)} {s.title}', 18, YELLOW)
        self.text(W - 16, 8, 'ION ' * self.lives, 18, MAGENTA, 'right')
        if s.message:
            self.text(W / 2, H - 42, s.message, 16, YELLOW)

    # ---------- screens

    def draw_attract(self):
        self.spiral(W / 2, 300, 135)
        self.text(W / 2, 40, 'ION RUSH', 58, CYAN)
        self.text(W / 2, 128, 'A CYCLOTRON ARCADE, ON THE LIONS PHYSICS', 14, MAGENTA)
        y = 450
        self.text(W / 2, y, 'HIGH SCORES', 18, YELLOW)
        for row in self.store.scores()[:5]:
            y += 28
            self.text(W / 2, y, f"{row['initials']}  {row['points']:07d}", 16)
        if int(self.clock * 2) % 2:
            self.text(W / 2, H - 56, 'PRESS START', 22, GREEN)

    def spiral(self, cx, cy, radius):
        turns, n = 5, 400
        phase = self.clock * 1.2 % 1
        path = QPainterPath(QPointF(cx, cy))
        for i in range(1, n + 1):
            a, r = 2 * math.pi * turns * i / n, radius * math.sqrt(i / n)
            path.lineTo(cx + r * math.cos(a), cy + r * math.sin(a))
        self.neon(lambda: self.p.drawPath(path), DIM, 1.5)
        at = (self.clock * 0.08) % 1
        a, r = 2 * math.pi * turns * at, radius * math.sqrt(at)
        self.ion(cx + r * math.cos(a), cy + r * math.sin(a), MAGENTA if phase < 0.5 else CYAN)

    def draw_consent(self):
        self.text(W / 2, 170, 'HELP THE RESEARCH?', 32, CYAN)
        y = self.block(270, 'Log your runs on this machine: knobs, readings, scores, with a random id. No name and no initials. '
                            'They can be used to compare human tuning with optimisers and policies.', 56)
        self.text(W / 2, y + 60, 'Y: YES     N: NO', 24, GREEN)

    def draw_menu(self):
        self.text(W / 2, 50, 'THE SHIFT', 40, CYAN)
        self.text(W / 2, 110, f'SCORE {self.total:07d}   IONS {self.lives}', 16, MAGENTA)
        y = 170
        for i, (key, title, stages) in enumerate(self.menu_items()):
            selected = i == self.choice
            color = GREEN if key in self.completed else YELLOW if selected else WHITE
            if selected:
                self.p.fillRect(QRectF(90, y - 6, W - 180, 44), QColor(255, 225, 77, 26))
            self.text(110, y, key, 18, color, 'left')
            self.text(260, y, title, 18, color, 'left')
            if key in self.completed:
                self.text(W - 110, y, 'DONE', 16, GREEN, 'right')
            y += 58
        request = self.menu_items()[self.choice][2][0]
        client = getattr(request, 'client', request.brief if self.choice else 'Five stages, from the Lorentz force to Thomas focusing.')
        self.block(H - 170, client, 66, 14, CYAN)
        self.text(W / 2, H - 52, 'UP DOWN: CHOOSE   START: ACCEPT   ESC: END THE SHIFT', 13, DIM)

    def draw_brief(self):
        s = self.stages[self.index]
        mission = isinstance(s, Mission)
        self.text(W / 2, 100, s.mission_id if mission else f'STAGE {self.index + 1}', 20, MAGENTA)
        self.text(W / 2, 140, s.title, 40, CYAN)
        y = self.block(230, s.client, 56, 17, YELLOW) + 18 if mission else 250
        y = self.block(y, s.brief, 56, 16)
        self.text(W / 2, y + 40, MISSION_CONTROLS if mission else CONTROLS[s.key], 13, GREEN)
        if int(self.clock * 2) % 2:
            self.text(W / 2, H - 80, 'PRESS START', 22, GREEN)

    def draw_clear(self):
        s = self.stage
        self.text(W / 2, 130, s.title, 28, YELLOW)
        self.text(W / 2, 180, s.message, 22, GREEN if s.success else RED)
        self.text(W / 2, 240, f'+{s.points} POINTS', 24, CYAN)
        self.block(330, s.lesson, 54, 18)
        self.text(W / 2, H - 80, 'PRESS START', 20, GREEN)

    def draw_lost(self):
        s = self.stage
        self.text(W / 2, 200, s.message, 28, RED)
        self.text(W / 2, 280, f'{self.lives} IONS LEFT', 22, MAGENTA)
        self.block(370, s.brief, 56, 15, DIM)
        self.text(W / 2, H - 80, 'PRESS START', 20, GREEN)

    def draw_over(self):
        self.text(W / 2, 230, 'GAME OVER', 56, RED)
        self.text(W / 2, 340, f'SCORE {self.total:07d}', 26, CYAN)
        self.text(W / 2, H - 80, 'PRESS START', 20, GREEN)

    def draw_won(self):
        self.text(W / 2, 170, 'LEVEL 0 CLEAR', 46, GREEN)
        self.text(W / 2, 280, f'SCORE {self.total:07d}', 26, CYAN)
        self.block(370, 'Next: the real machine. Central region, deflector, the CIME layout. Needs the CIME map.', 54, 16, YELLOW)
        self.text(W / 2, H - 80, 'PRESS START', 20, GREEN)

    def draw_initials(self):
        self.text(W / 2, 170, 'NEW HIGH SCORE', 38, YELLOW)
        self.text(W / 2, 250, f'{self.total:07d}', 28, CYAN)
        for i, letter in enumerate(self.initials):
            self.text(W / 2 + (i - 1) * 70, 340, letter, 48, GREEN if i == self.letter else WHITE)
        self.text(W / 2, 480, 'UP DOWN: LETTER   SPACE: NEXT', 16, DIM)

    def draw_play(self):
        if isinstance(self.stage, Mission):
            self.control_room()
        else:
            getattr(self, f'view_{self.stage.key}')()
        self.hud()

    # ---------- control room

    def control_room(self):
        s = self.stage
        self.knob_panel(s)
        self.plot(350, 62, 580, 110, s.beam.probe_r, s.beam.probe, CYAN, 'RADIAL PROBE, CURRENT AGAINST RADIUS', 0)
        self.plot(350, 192, 580, 90, s.beam.probe_r, s.beam.phase, MAGENTA, 'PHASE PROBE, -90 TO +90 DEG', -90, 90, s.beam.probe > 0)
        r = s.readings
        y = 300
        for label, value, color in (('CURRENT', self.current_text(r['extracted_uA']), CYAN), ('ENERGY', f"{r['w_MeV']:.2f} MEV", WHITE),
                                    ('TRANSMISSION', f"{100 * r['transmission']:.1f}%", GREEN), ('LOST', f"{1e3 * r['lost_kW']:.1f} W", RED)):
            self.text(30, y + 325, label, 11, color, 'left')
            self.text(326, y + 323, value, 13, WHITE, 'right')
            y += 22
        self.bar(350, 300, 180, 16, s.heat, RED, '')
        self.text(350, 318, 'INTERLOCK', 11, RED, 'left')
        self.text(W - 30, 298, 'BEAM ON' if s.beam_on else 'BEAM OFF', 18, RED if s.beam_on else DIM, 'right')
        self.text(W / 2 + 120, 298, f'TIME {max(s.p["time_limit"] - s.time, 0):.0f}', 16, YELLOW)
        getattr(self, f'panel_{s.key}')(s, 350, 345, 580, 320)

    @staticmethod
    def current_text(ua):
        return f'{ua:.1f} UA' if ua >= 1 else f'{1e3 * ua:.2f} NA'

    def knob_panel(self, s):
        y = 60
        row = min(50, 560 / len(s.specs))
        for i, spec in enumerate(s.specs):
            selected = i == s.selected
            color = YELLOW if selected else WHITE
            if selected:
                self.p.fillRect(QRectF(16, y - 4, 318, row - 2), QColor(255, 225, 77, 28))
            value = s.knobs[spec['key']]
            digits = 3 if spec['unit'] in ('MHz', 'T') else 2
            self.text(24, y - 2, spec['label'], 12, color, 'left')
            self.text(326, y - 2, f"{value:.{digits}f} {spec['unit'] if spec['unit'] != '1' else ''}".strip(), 13, color, 'right')
            self.bar(24, y + 22, 302, 8, (value - spec['min']) / (spec['max'] - spec['min']), CYAN if selected else DIM)
            y += row

    def plot(self, x, y, width, height, xs, ys, color, title, low=None, high=None, mask=None):
        self.p.setPen(QPen(DIM, 1))
        self.p.setBrush(Qt.NoBrush)
        self.p.drawRect(QRectF(x, y, width, height))
        self.text(x + 6, y + 2, title, 10, DIM, 'left')
        top = max(ys.max(), 1e-12) if high is None else high
        bottom = 0 if low is None else low
        pos = lambda u, v: QPointF(x + width * (u - xs[0]) / (xs[-1] - xs[0]), y + height * (1 - (min(max(v, bottom), top) - bottom) / (top - bottom)))  # noqa: E731
        path, started = QPainterPath(), False
        for i, (u, v) in enumerate(zip(xs, ys)):
            if mask is not None and not mask[i]:
                started = False
                continue
            path.lineTo(pos(u, v)) if started else path.moveTo(pos(u, v))
            started = True
        self.neon(lambda: self.p.drawPath(path), color, 1.5)
        return pos

    def crosshair(self, x, y, radius, s, limit_mm, title, spot=None):
        scale = radius / (2 * limit_mm)
        self.circle(x, y, limit_mm * scale, GREEN)
        self.line(x - radius, y, x + radius, y, DIM, 1)
        self.line(x, y - radius, x, y + radius, DIM, 1)
        px, py = s.readings['x_mm'] * scale, s.readings['y_mm'] * scale
        px, py = max(min(px, radius), -radius), max(min(py, radius), -radius)
        self.circle(x + px, y - py, max((spot or s.spot()) * scale, 3), MAGENTA, 1.5)
        self.ion(x + px, y - py, YELLOW)
        self.text(x, y - radius - 24, title, 11, CYAN)

    def lights(self, x, y, rows):
        for label, good in rows:
            self.p.setPen(Qt.NoPen)
            self.p.setBrush(GREEN if good else RED)
            self.p.drawEllipse(QPointF(x, y + 10), 6, 6)
            self.p.setBrush(Qt.NoBrush)
            self.text(x + 16, y, label, 12, WHITE, 'left')
            y += 26

    def panel_therapy(self, s, x, y, width, height):
        p = s.p
        d0, d1 = s.depth_axis[0], s.depth_axis[-1]
        w = width - 200
        low = x + w * (p['depth_cm'] - p['half_width_cm'] - d0) / (d1 - d0)
        high = x + w * (p['depth_cm'] + p['half_width_cm'] - d0) / (d1 - d0)
        self.p.fillRect(QRectF(low, y, max(high - low, 2), 170), QColor(61, 255, 122, 70))
        self.plot(x, y, w, 170, s.depth_axis, s.bragg, YELLOW, 'DOSE AGAINST DEPTH IN WATER; GREEN: TUMOUR')
        self.text(x, y + 176, f'PEAK {s.peak:.2f} CM   TUMOUR {p["depth_cm"]:.1f} CM   SPOT {s.spot():.1f} MM', 13, WHITE, 'left')
        self.crosshair(x + width - 90, y + 90, 80, s, p['lateral_mm'], 'LATERAL')
        self.bar(x, y + 230, width, 22, s.delivered / (p['dose_gy'] * (1 + p['dose_tolerance'])), GREEN,
                 f"TUMOUR {s.delivered:.1f} / {p['dose_gy']:.1f} GY   RATE {s.rate:.2f} GY/S")
        ok = x + width * (1 - p['dose_tolerance']) / (1 + p['dose_tolerance'])
        self.line(ok, y + 226, ok, y + 256, YELLOW)
        self.bar(x, y + 288, width, 14, s.healthy / p['healthy_limit_gy'], RED, f'OUTSIDE THE TUMOUR {s.healthy:.2f} GY')

    def panel_isotope(self, s, x, y, width, height):
        p = s.p
        self.bar(x, y + 28, width - 200, 24, s.activity / p['activity_gbq'], GREEN, f"F-18 {s.activity:.1f} / {p['activity_gbq']:.0f} GBQ")
        self.bar(x, y + 92, width - 200, 18, s.target_kw / p['target_kw'], RED if s.target_kw > 0.9 * p['target_kw'] else CYAN,
                 f"POWER ON THE WATER {s.target_kw:.2f} / {p['target_kw']:.1f} KW")
        limit = s.machine.max_lost_kw
        self.bar(x, y + 150, width - 200, 18, s.beam.lost_kw / limit, RED, f'LOST IN THE MACHINE {s.beam.lost_kw:.2f} / {limit:.1f} KW')
        self.text(x, y + 190, f"SATURATION {p['yield_gbq_ua']:.0f} GBQ/UA AT {p['w_ref']:.0f} MEV", 12, DIM, 'left')
        self.crosshair(x + width - 90, y + 90, 80, s, p['window_mm'], 'TARGET WINDOW', spot=None)

    def panel_experiment(self, s, x, y, width, height):
        p, b = s.p, s.beam
        energy = abs(b.w_mean / p['w_target'] - 1) <= p['w_tolerance'] if b.w_mean else False
        self.lights(x, y, ((f"ENERGY {b.w_mean / s.ion.mass_number:.3f} MEV/U, WANT {p['w_target'] / s.ion.mass_number:.2f}", energy),
                           (f"SPREAD {100 * b.w_spread:.2f}%, MAX {100 * p['spread_max']:.1f}%", b.w_mean > 0 and b.w_spread <= p['spread_max']),
                           (f"CURRENT {b.extracted_ua:.2f} UA, {p['current_min']:.1f} TO {p['current_max']:.1f}", p['current_min'] <= b.extracted_ua <= p['current_max']),
                           (f"ON AXIS {math.hypot(b.x, b.y):.1f} MM, MAX {p['position_mm']:.0f}", math.hypot(b.x, b.y) <= p['position_mm']),
                           (f"SPOT {b.spot:.1f} MM, MAX {p['spot_max']:.1f}", b.spot <= p['spot_max'])))
        self.crosshair(x + width - 90, y + 90, 80, s, p['position_mm'], 'ON TARGET')
        self.bar(x, y + 190, width, 22, s.delivered / p['hold_s'], GREEN, f"BEAM IN SPEC {s.delivered:.1f} / {p['hold_s']:.0f} S   {s.ion.name}   HARMONIC {b.harmonic}")

    def panel_isol(self, s, x, y, width, height):
        p = s.p
        self.bar(x, y + 28, width - 200, 22, s.delivered / p['energy_goal_kws'], GREEN, f"ON TARGET {s.delivered:.1f} / {p['energy_goal_kws']:.0f} KJ")
        self.bar(x, y + 90, width - 200, 18, s.power() / p['power_limit_kw'], CYAN, f"POWER {s.power():.2f} / {p['power_limit_kw']:.1f} KW")
        self.bar(x, y + 148, width - 200, 18, s.density() / p['density_limit'], RED, f"POWER DENSITY {1e3 * s.density():.1f} / {1e3 * p['density_limit']:.0f} W/MM2")
        self.crosshair(x + width - 90, y + 90, 80, s, p['target_radius_mm'], 'TARGET')

    # ---------- stage views

    def view_larmor(self):
        s = self.stage
        cx, cy, scale = W / 2 - 110, H / 2 + 20, 255 / s.p['target_max']
        r = s.target * scale
        self.circle(cx, cy, r, GREEN, 2)
        self.p.setPen(QPen(GREEN, 2))
        tip = s.target_rotation * 0.35
        self.arrow(cx + r * math.cos(tip), cy - r * math.sin(tip), tip + s.target_rotation * math.pi / 2, GREEN)
        if s.shot_radius is not None:
            rs = s.shot_radius * scale
            fraction = 1 if s.flight is None else s.flight
            span = int(s.shot_rotation * fraction * 360 * 16)
            rect = QRectF(cx - rs, cy - rs, 2 * rs, 2 * rs)
            self.neon(lambda: self.p.drawArc(rect, 0, span), MAGENTA if s.flight is not None else DIM)
            if s.flight is not None:
                a = s.shot_rotation * 2 * math.pi * s.flight
                self.ion(cx + rs * math.cos(a), cy - rs * math.sin(a), MAGENTA)
        x = W - 210
        self.text(x, 90, 'FIELD B', 15, CYAN, 'left')
        self.text(x, 112, f'{s.b:+.2f} T', 28, WHITE, 'left')
        self.text(x, 180, 'ENERGY', 15, CYAN, 'left')
        self.text(x, 202, f'{s.w:.1f} MEV', 28, WHITE, 'left')
        self.text(x, 280, 'TARGET', 15, GREEN, 'left')
        self.text(x, 302, f'{s.target:.2f} M', 28, WHITE, 'left')
        self.text(x, 370, 'RINGS', 15, YELLOW, 'left')
        self.text(x, 392, f"{s.hits}/{s.p['targets']}", 28, WHITE, 'left')
        self.text(x, 460, 'SHOTS LEFT', 15, YELLOW, 'left')
        self.text(x, 482, f"{s.p['shots'] - s.shots}", 28, WHITE, 'left')
        self.text(x, 560, 'RHO = P / (Q B)', 13, DIM, 'left')

    def arrow(self, x, y, angle, color):
        for side in (2.6, -2.6):
            self.line(x, y, x + 16 * math.cos(angle + side), y - 16 * math.sin(angle + side), color)

    def gap_view(self, reference_radius):
        s = self.stage
        cx, cy, radius = 330, H / 2 + 20, 280
        scale = radius / reference_radius
        for sign, color in ((-1, CYAN), (1, MAGENTA)):
            path = QPainterPath(QPointF(cx - radius - 12, cy + sign * 8))
            path.lineTo(cx + radius + 12, cy + sign * 8)
            rect = QRectF(cx - radius - 12, cy - radius - 12 + sign * 8, 2 * radius + 24, 2 * radius + 24)
            path.arcTo(rect, 0, sign * -180)
            path.closeSubpath()
            dim = QColor(color)
            dim.setAlpha(150)
            self.neon(lambda path=path: self.p.drawPath(path), dim, 1.5)
        at = len(s.arcs) - 1
        visible = s.arcs[-90:]
        for i, (x0, x1) in enumerate(visible):
            n = at - len(visible) + 1 + i
            rect = QRectF(cx + scale * min(x0, x1), cy - scale * abs(x1 - x0) / 2, scale * abs(x1 - x0), scale * abs(x1 - x0))
            shade = QColor(WHITE)
            shade.setAlpha(40 + int(180 * i / len(visible)))
            self.p.setPen(QPen(shade, 1.2))
            self.p.setBrush(Qt.NoBrush)
            self.p.drawArc(rect, 0, 180 * 16 if n % 2 == 0 else -180 * 16)
        x0, x1 = s.arcs[-1]
        c, r = (x0 + x1) / 2, abs(x1 - x0) / 2
        progress = s.frac * math.pi
        side = -1 if at % 2 == 0 else 1
        a = progress if x1 < x0 else math.pi - progress
        self.ion(cx + scale * (c + r * math.cos(a)), cy + side * scale * r * math.sin(a), YELLOW)

    def phase_dial(self, x, y, ghost=None):
        s = self.stage
        r = 72
        self.circle(x, y, r, DIM)
        crest = math.degrees(math.acos(s.p['crest_cos']))
        rect = QRectF(x - r, y - r, 2 * r, 2 * r)
        self.neon(lambda: self.p.drawArc(rect, int((90 - crest) * 16), int(2 * crest * 16)), GREEN, 4)
        if ghost is not None:
            self.line(x, y, x + r * 0.8 * math.sin(ghost), y - r * 0.8 * math.cos(ghost), DIM)
        color = GREEN if math.cos(s.phi) >= s.p['crest_cos'] else RED if math.cos(s.phi) < 0 else YELLOW
        self.line(x, y, x + r * 0.9 * math.sin(s.phi), y - r * 0.9 * math.cos(s.phi), color, 3)
        self.text(x, y - r - 30, 'RF PHASE', 14, CYAN)

    def view_dee_beat(self):
        s = self.stage
        self.gap_view(physics.radius(s.p['w_target'], s.p['b']))
        x = W - 150
        self.phase_dial(x, 170, None)
        self.text(x, 280, 'RF', 14, CYAN)
        self.text(x, 302, f'{s.f / 1e6:.3f} MHZ', 20)
        self.text(x, 370, 'COMBO', 14, YELLOW)
        self.text(x, 392, f'x{min(s.combo, 20)}', 26, YELLOW)
        self.bar(x - 110, 490, 220, 22, s.w / s.p['w_target'], GREEN, f'{s.w:.2f} MEV')
        self.text(x, 540, f"TIME {s.p['time_limit'] - s.time:.0f}", 15)

    def view_trim_coils(self):
        s = self.stage
        gamma_ratio = physics.gamma(s.p['w_target']) / physics.gamma(s.p['w_start'])
        self.gap_view(physics.radius(s.p['w_target'], s.p['b'] * gamma_ratio))
        x = W - 150
        self.phase_dial(x, 170, s.ghost_phi)
        self.text(x, 262, 'GREY: UNTRIMMED', 11, DIM)
        self.text(x, 300, 'FIELD TRIM', 14, CYAN)
        self.text(x, 322, f'{100 * s.delta:+.2f}%', 24)
        self.text(x, 370, 'GAMMA', 14, MAGENTA)
        self.text(x, 392, f'{physics.gamma(s.w):.3f}', 24)
        self.bar(x - 110, 490, 220, 22, s.w / s.p['w_target'], GREEN, f'{s.w:.1f} MEV')
        self.text(x, 540, f"TIME {s.p['time_limit'] - s.time:.0f}", 15)

    def tunnel(self, top=110, bottom=470):
        s = self.stage
        mid, scale = (top + bottom) / 2, (bottom - top) / 2 / s.p['gap']
        for y in (top, bottom):
            self.line(40, y, W - 40, y, CYAN, 3)
        self.line(40, mid, W - 40, mid, DIM, 1)
        if s.trail:
            path = QPainterPath(QPointF(W - 120 - len(s.trail), mid - s.trail[0] * scale))
            for i, z in enumerate(s.trail):
                path.lineTo(W - 120 - len(s.trail) + i, mid - max(min(z, s.p['gap']), -s.p['gap']) * scale)
            self.neon(lambda: self.p.drawPath(path), MAGENTA)
        self.ion(W - 120, mid - max(min(s.z, s.p['gap']), -s.p['gap']) * scale, YELLOW)

    def view_fallout(self):
        s = self.stage
        self.tunnel()
        k = physics.isochronous_index(s.w)
        self.text(60, 505, f'QZ^2 = -K = {-k:+.3f}', 18, RED, 'left')
        self.text(60, 540, f'{s.w:.1f} MEV', 18, WHITE, 'left')
        self.text(W / 2, 505, f"TURN {s.turns:.1f}/{s.p['turns_target']}", 18, YELLOW)
        self.text(W - 60, 505, 'CHARGES', 15, CYAN, 'right')
        self.bar(W - 280, 535, 220, 22, s.charges / s.p['charges'], CYAN)

    def view_thomas(self):
        s = self.stage
        self.tunnel()
        self.bar(60, 535, 200, 20, s.flutter / s.p['f_max'], CYAN, f'FLUTTER F {s.flutter:.2f}')
        self.bar(60, 600, 200, 20, s.xi_deg / s.p['xi_max'], MAGENTA, f'SPIRAL XI {s.xi_deg:.0f} DEG')
        self.bar(W - 260, 535, 200, 20, s.w / s.p['w_target'], GREEN, f'{s.w:.0f} MEV')
        x0, x1, y, top = 330, 680, 600, 1.3
        at = lambda q: x0 + (x1 - x0) * min(q, top) / top  # noqa: E731
        self.p.fillRect(QRectF(at(s.p['band_low']), y - 4, at(s.p['band_high']) - at(s.p['band_low']), 24), QColor(61, 255, 122, 70))
        self.line(x0, y + 8, x1, y + 8, DIM)
        for q, color in ((0.5, RED), (1.0, RED)):
            self.line(at(q), y - 10, at(q), y + 26, color)
            self.text(at(q), y + 28, f'{q:g}', 12, color)
        if s.qz2 > 0:
            self.ion(at(s.qz), y + 8, YELLOW)
            label = f'QZ = {s.qz:.2f}'
        else:
            self.line(x0, y - 14, x0, y + 30, RED, 3)
            label = 'QZ^2 < 0: DEFOCUSING'
        self.text((x0 + x1) / 2, y - 50, label, 16, YELLOW if s.qz2 > 0 else RED)


def main():
    app = QApplication(sys.argv)
    window = Arcade()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
