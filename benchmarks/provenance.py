"""Code and input provenance recorded in every benchmark table."""
import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def git(*args):
    return subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True, text=True,
                          check=True).stdout.strip()


def provenance(*inputs):
    """code_branch, code_version (git describe), code_commit, and path:sha256 of each input file."""
    out = {'code_branch': git('rev-parse', '--abbrev-ref', 'HEAD'),
           'code_version': git('describe', '--always', '--dirty', '--tags'),
           'code_commit': git('rev-parse', 'HEAD')}
    for path in inputs:
        out[f'input_sha256:{path}'] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return out
