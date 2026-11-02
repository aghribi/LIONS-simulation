#!/bin/bash
#SBATCH --job-name=lions_case_p
#SBATCH --output=benchmarks/runs/slurm/%x_%j.log
# From the repo root, environment synced on the login node. Site, account, size and time at submission: sbatch -A <acct> -p <part> -c 2 --mem=4G -t 01:00:00 benchmarks/case_p.sbatch
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH" OMP_NUM_THREADS=1
git log --oneline -1
RUN="uv run --frozen --extra game --group bench --group dev"
for step in design dees slip offaxis; do
    echo "== $step"
    $RUN python benchmarks/case_p.py "$step"
done
