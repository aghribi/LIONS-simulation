import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import synthetic  # noqa: E402


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    """Empty directory with zero central and outer potentials, set as cwd."""
    monkeypatch.chdir(tmp_path)
    synthetic.write_pot_central('pot1.bin')
    synthetic.write_pot_outer('pot2.dat')
    return tmp_path
