# Job image: LIONS dependencies and the gpu group (torch, pyarrow), from environment/requirements.gpu.txt
# (uv export of uv.lock, hash-pinned; tests/test_environment.py keeps it in sync). Built by
# .github/workflows/env-image.yml, tagged with tools/env-hash. The LIONS code comes from the site clone (PYTHONPATH).
# Pinned: base image and uv by digest, Python packages by hash.
FROM python:3.12-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337
COPY --from=ghcr.io/astral-sh/uv:0.13.0@sha256:cdc6093146eb3ff6a40107b38f008b789e050e77ad87865e381d9917da55a168 /uv /usr/local/bin/uv

# torch.compile builds C++ (CPU) and Triton launchers (GPU) at run time; git gives the code version.
RUN apt-get update && apt-get install -y --no-install-recommends gcc g++ libc6-dev git && rm -rf /var/lib/apt/lists/*

COPY environment/Dockerfile environment/requirements.gpu.txt /opt/lock/
ENV UV_PYTHON_DOWNLOADS=never UV_LINK_MODE=copy
RUN uv venv --python /usr/local/bin/python3 /opt/venv \
 && uv pip install --python /opt/venv/bin/python --require-hashes --no-deps -r /opt/lock/requirements.gpu.txt \
 && uv cache clean

# Identity of the image, checked by jobs/_job_env.sh against the clone's recipe.
RUN cat /opt/lock/Dockerfile /opt/lock/requirements.gpu.txt | sha256sum | cut -c1-12 > /opt/lock/HASH \
 && uv pip freeze --python /opt/venv/bin/python > /opt/lock/freeze.txt

ENV PATH=/opt/venv/bin:$PATH \
    PYTHONNOUSERSITE=1
