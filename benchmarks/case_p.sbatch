#!/bin/bash
#SBATCH --job-name=lions_case_p
#SBATCH --account=ganil
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/runs/slurm/%x_%j.log
# From the repo root on CC-IN2P3, which sets a dead proxy in the login profile: sbatch benchmarks/case_p.sbatch
set -euo pipefail
unset HTTPS_PROXY HTTP_PROXY https_proxy http_proxy
export PATH="$HOME/.local/bin:$PATH" OMP_NUM_THREADS=1
git log --oneline -1
RUN="uv run --frozen --extra game --group bench --group dev"
for step in design dees slip offaxis; do
    echo "== $step"
    $RUN python benchmarks/case_p.py "$step"
done
