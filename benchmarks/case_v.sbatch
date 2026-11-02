#!/bin/bash
#SBATCH --job-name=lions_case_v
#SBATCH --output=benchmarks/runs/slurm/%x_%j.log
# From the repo root, environment synced on the login node. Site, account, size and time at submission: sbatch -p <part> -c <n> --mem=<m> -t <t> benchmarks/case_v.sbatch "tunes check"
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH" OMP_NUM_THREADS=1
git log --oneline -1
RUN="${PY:-uv run --frozen --extra game --extra rl --group bench --group dev python}"
for step in $1; do
    echo "== $step"
    $RUN benchmarks/case_v.py "$step"
done
