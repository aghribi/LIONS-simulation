"""The job image installs environment/requirements.gpu.txt: it must be the uv export of uv.lock."""

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(shutil.which('uv') is None, reason='uv not on PATH')
def test_gpu_requirements_match_lock():
    export = subprocess.run(
        [
            'uv',
            'export',
            '--frozen',
            '--no-default-groups',
            '--group',
            'gpu',
            '--no-emit-project',
            '--no-header',
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert (
        export == (ROOT / 'environment' / 'requirements.gpu.txt').read_text()
    ), (
        'regenerate: uv export --frozen --no-default-groups --group gpu --no-emit-project --no-header '
        '-o environment/requirements.gpu.txt'
    )
