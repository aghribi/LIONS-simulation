import pytest

import synthetic


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    """Empty directory with zero central and outer potentials, set as
    cwd."""
    monkeypatch.chdir(tmp_path)
    synthetic.write_pot_central('pot1.bin')
    synthetic.write_pot_outer('pot2.dat')
    return tmp_path
