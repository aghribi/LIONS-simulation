"""Code and input provenance recorded in every benchmark table."""
import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def git(*args):
    return subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True, text=True,
                          check=True).stdout.strip()


# Read once at import: a run records the code it started with, even if a commit lands before it writes its table.
CODE = {'code_branch': git('rev-parse', '--abbrev-ref', 'HEAD'),
        'code_version': git('describe', '--always', '--dirty', '--tags'),
        'code_commit': git('rev-parse', 'HEAD')}


def provenance(*inputs):
    """code_branch, code_version (git describe), code_commit at process start, and path:sha256 of each input file."""
    out = dict(CODE)
    for path in inputs:
        out[f'input_sha256:{path}'] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return out
