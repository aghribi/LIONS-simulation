"""Run index of a site results store: one row per run directory
<case>/<run>/, written to <store>/run_index.parquet.

usage: python benchmarks/build_run_index.py STORE [--site SITE]
[--dry-run] A run is complete when it holds a case_*.parquet; its code
and input provenance come from that table's schema metadata
(benchmarks/provenance.py). Runs without one are listed as incomplete.
"""

import argparse
import datetime as dt
import os
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

BUILDER_VERSION = '1'
FIELDS = [
    ('site', pa.string()),
    ('case', pa.string()),
    ('run', pa.string()),
    ('path', pa.string()),
    ('status', pa.string()),
    ('table', pa.string()),
    ('code_branch', pa.string()),
    ('code_version', pa.string()),
    ('code_commit', pa.string()),
    ('device_name', pa.string()),
    ('machine', pa.string()),
    ('inputs', pa.string()),
    ('n_files', pa.int64()),
    ('size_bytes', pa.int64()),
    ('modified', pa.timestamp('s', tz='UTC')),
]


def run_row(site, store, run):
    """Return the index row of one run directory of a store."""
    files = [f for f in run.rglob('*') if f.is_file()]
    tables = sorted(run.glob('case_*.parquet'))
    meta = {}
    if tables:
        raw = pq.read_schema(tables[0]).metadata or {}
        meta = {k.decode(): v.decode() for k, v in raw.items()}
    mtime = max(
        (f.stat().st_mtime for f in files), default=run.stat().st_mtime
    )
    return {
        'site': site,
        'case': run.parent.name,
        'run': run.name,
        'path': str(run.relative_to(store)),
        'status': 'complete' if tables else 'incomplete',
        'table': tables[0].name if tables else None,
        'code_branch': meta.get('code_branch'),
        'code_version': meta.get('code_version'),
        'code_commit': meta.get('code_commit'),
        'device_name': meta.get('device_name'),
        'machine': meta.get('machine'),
        'inputs': ';'.join(
            f'{k.removeprefix("input_sha256:")}={v}'
            for k, v in sorted(meta.items())
            if k.startswith('input_sha256:')
        )
        or None,
        'n_files': len(files),
        'size_bytes': sum(f.stat().st_size for f in files),
        'modified': dt.datetime.fromtimestamp(mtime, dt.timezone.utc),
    }


def main():
    """Write the run index of a store (--dry-run: print the count)."""
    parser = argparse.ArgumentParser()
    parser.add_argument('store', type=Path)
    parser.add_argument('--site', default=os.environ.get('XSH_SITE', 'local'))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    store = args.store.resolve()
    runs = sorted(
        r
        for case in store.glob('case_*')
        if case.is_dir()
        for r in case.iterdir()
        if r.is_dir()
    )
    rows = [run_row(args.site, store, r) for r in runs]
    table = pa.Table.from_pylist(
        rows, schema=pa.schema(FIELDS)
    ).replace_schema_metadata(
        {
            'builder': 'benchmarks/build_run_index.py',
            'builder_version': BUILDER_VERSION,
            'store': str(store),
            'site': args.site,
            'built': dt.datetime.now(dt.timezone.utc).isoformat(
                timespec='seconds'
            ),
            'units': 'size_bytes: bytes; modified: UTC',
        }
    )
    print(
        f'{len(rows)} runs '
        f'({sum(r["status"] == "complete" for r in rows)} complete) '
        f'in {store}'
    )
    if args.dry_run:
        return
    fd, tmp = tempfile.mkstemp(dir=store, suffix='.parquet.tmp')
    os.close(fd)
    pq.write_table(table, tmp)
    os.replace(tmp, store / 'run_index.parquet')


if __name__ == '__main__':
    main()
