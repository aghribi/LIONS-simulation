"""Local files of the trainer: consent, high scores, anonymous run logs (SPEC section 8, decision D6).

A run log holds observations, actions and the score record, with a random run id. It never holds the initials of the
high-score table. Nothing is written without consent.
"""
import json
import subprocess
import uuid
from importlib import metadata
from pathlib import Path

HOME = Path.home() / '.lions_arcade'
TABLE_SIZE = 8
UNITS = {'_T': 'T', '_MeV': 'MeV', '_m': 'm', '_rad': 'rad', '_Hz': 'Hz', '_mm': 'mm', '_mm_rad': 'mm/rad', '_deg': 'deg'}


def code_version():
    try:
        return subprocess.run(['git', '-C', str(Path(__file__).parent), 'describe', '--always', '--dirty'], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return metadata.version('lions')


class Store:
    def __init__(self, home=HOME):
        self.home = Path(home)
        self.home.mkdir(parents=True, exist_ok=True)
        self.session = uuid.uuid4().hex

    def _read(self, name, default):
        path = self.home / name
        return json.loads(path.read_text()) if path.exists() else default

    def _write(self, name, value):
        (self.home / name).write_text(json.dumps(value, indent=1))

    @property
    def consent(self):
        """None until asked."""
        return self._read('settings.json', {}).get('log_runs')

    @consent.setter
    def consent(self, value):
        self._write('settings.json', {'log_runs': bool(value)})

    def scores(self):
        return self._read('scores.json', [])

    def qualifies(self, points):
        table = self.scores()
        return points > 0 and (len(table) < TABLE_SIZE or points > table[-1]['points'])

    def add_score(self, initials, points):
        table = sorted([*self.scores(), {'initials': initials, 'points': points}], key=lambda r: -r['points'])
        self._write('scores.json', table[:TABLE_SIZE])

    def log_run(self, stage, rows):
        """One Parquet file per stage run. rows: one dict per step. Does nothing without consent."""
        if not self.consent or not rows:
            return None
        import pyarrow as pa
        import pyarrow.parquet as pq
        record = stage.score()
        units = {k: next((u for suffix, u in UNITS.items() if k.endswith(suffix)), '1') for k in rows[0]}
        meta = {'score': json.dumps(record), 'units': json.dumps(units), 'code_version': code_version(), 'session': self.session,
                'source': 'lions.game arcade, human player, anonymous', 'time_step_s': json.dumps(stage.time / max(stage.steps, 1))}
        columns = {k: [r[k] for r in rows] for k in rows[0]}
        path = self.home / 'runs' / f'{self.session}_{stage.key}_{stage.seed}.parquet'
        path.parent.mkdir(exist_ok=True)
        pq.write_table(pa.table(columns).replace_schema_metadata(meta), path)
        return path
