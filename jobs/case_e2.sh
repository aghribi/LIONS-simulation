#!/bin/bash
# XSH name=case_e2 time=01:00:00 ntasks=1 cpus=8 mem=32G
# Case E2 speed spike (benchmarks/DESIGN.md section 20). Submit with --class gpu. Parameter: DEVICE (default cuda).
source _job_env.sh || exit 3
python ../benchmarks/case_e2.py "$XSH_RESULTS/case_e2/${XSH_JOB_ID}_$(hostname -s)" --device "${DEVICE:-cuda}"
