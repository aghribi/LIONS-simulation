# Job image: LIONS dependencies and the gpu group (torch, pyarrow) from uv.lock. Built by
# .github/workflows/env-image.yml, tagged with tools/env-hash. The LIONS code comes from the site clone (PYTHONPATH).
# Pinned: base image and uv by digest, Python packages by uv.lock.
FROM python:3.12-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337
COPY --from=ghcr.io/astral-sh/uv:0.13.0@sha256:cdc6093146eb3ff6a40107b38f008b789e050e77ad87865e381d9917da55a168 /uv /usr/local/bin/uv

# torch.compile builds C++ (CPU) and Triton launchers (GPU) at run time; git gives the code version.
RUN apt-get update && apt-get install -y --no-install-recommends gcc g++ libc6-dev git && rm -rf /var/lib/apt/lists/*

COPY environment/Dockerfile pyproject.toml uv.lock /opt/lock/
ENV UV_PROJECT_ENVIRONMENT=/opt/venv UV_PYTHON_DOWNLOADS=never UV_LINK_MODE=copy
RUN cd /opt/lock && uv sync --frozen --no-install-project --no-default-groups --group gpu && uv cache clean

# Identity of the image, checked by jobs/_job_env.sh against the clone's recipe.
RUN cat /opt/lock/Dockerfile /opt/lock/pyproject.toml /opt/lock/uv.lock | sha256sum | cut -c1-12 > /opt/lock/HASH \
 && uv pip freeze --python /opt/venv/bin/python > /opt/lock/freeze.txt

ENV PATH=/opt/venv/bin:$PATH \
    PYTHONNOUSERSITE=1
