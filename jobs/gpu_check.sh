#!/bin/bash
# XSH name=gpu_check time=00:15:00 ntasks=1 cpus=4 mem=16G
# Site check before GPU runs: image, GPU, and the case E port in float64 on the GPU against the CPU (360 steps).
source _job_env.sh || exit 3
cd "$TMPDIR" 2>/dev/null || cd "$(mktemp -d)"
python - <<'PY'
import sys
sys.path.insert(0, f'{__import__("os").environ["XSH_ROOT"]}/benchmarks')
import torch
import case_e2 as c
print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), torch.cuda.get_device_name() if torch.cuda.is_available() else '')
c.write_case('B')
d = c.load(360)
out = {}
for dev in ('cpu', 'cuda'):
    c.DEVICE = dev
    with torch.no_grad():
        out[dev] = c.make_port(d, torch.float64, False).run(c.start(d, 1000, torch.float64), 1.0, 360)[0].cpu()
print('cuda vs cpu float64, 360 steps, max rel r:', ((out['cuda'] - out['cpu']).abs() / out['cpu']).max().item())
PY
