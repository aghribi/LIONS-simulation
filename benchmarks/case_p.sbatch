#!/bin/bash
#SBATCH --job-name=lions_case_p
#SBATCH --output=benchmarks/runs/slurm/%x_%j.log
# From the repo root, environment synced on the login node. Site, account, size and time at submission: sbatch -A <acct> -p <part> -c 2 --mem=4G -t 01:00:00 benchmarks/case_p.sbatch. PY=<env>/bin/python picks another environment, CASE_P_MISSION the machine, OUTROOT the output.
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH" OMP_NUM_THREADS=1
git log --oneline -1
RUN="${PY:-uv run --frozen --extra game --extra rl --group bench --group dev python}"
for step in design dees slip offaxis published; do
    echo "== $step"
    $RUN benchmarks/case_p.py "$step" ${OUTROOT:+--outroot "$OUTROOT"}
done
