# Sourced by every job script, from jobs/: site settings, container image, job identity.
# tools/xsh-submit exports XSH_SITE; sites/<site>.env gives the paths and the container runtime.
: "${XSH_SITE:?XSH_SITE is not set: submit with tools/xsh-submit}"
XSH_ROOT="$(cd .. && pwd)"
source "$XSH_ROOT/sites/$XSH_SITE.env"
[ -f "$XSH_ROOT/sites/$XSH_SITE.local.env" ] && source "$XSH_ROOT/sites/$XSH_SITE.local.env"
eval "${CONTAINER_SETUP:-}"

# The image must be the one built from this commit's recipe (environment/, pyproject.toml, uv.lock).
export XSH_ENV_HASH="$(sh "$XSH_ROOT/tools/env-hash")"
export XSH_IMAGE="$IMAGE_DIR/lions-sim-$XSH_ENV_HASH.sif"
[ -f "$XSH_IMAGE" ] || { echo "no image $XSH_IMAGE: run tools/site.sh $XSH_SITE setup" >&2; exit 3; }
_img_hash="$("$CONTAINER" exec "$XSH_IMAGE" cat /opt/lock/HASH)"
[ "$_img_hash" = "$XSH_ENV_HASH" ] || { echo "image $XSH_IMAGE has recipe $_img_hash, code expects $XSH_ENV_HASH" >&2; exit 3; }

mkdir -p "$XSH_RESULTS"   # bind sources must exist
_binds="$XSH_ROOT,$CLONE,$XSH_RESULTS${TMPDIR:+,$TMPDIR}"   # CLONE: git data of worktrees
_nv=""; command -v nvidia-smi >/dev/null && nvidia-smi -L >/dev/null 2>&1 && _nv="--nv"
python() { "$CONTAINER" exec $_nv --bind "$_binds" "$XSH_IMAGE" python "$@"; }

export XSH_SITE XSH_RESULTS
export XSH_JOB_ID="${SLURM_JOB_ID:-local}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export PYTHONPATH="$XSH_ROOT/src" PYTHONUNBUFFERED=1
export TORCHINDUCTOR_CACHE_DIR="${TMPDIR:-/tmp}/inductor" TRITON_CACHE_DIR="${TMPDIR:-/tmp}/triton"
export MPLCONFIGDIR="${TMPDIR:-/tmp}"   # home may be read-only inside the image
mkdir -p logs
echo "site=$XSH_SITE job=$XSH_JOB_ID host=$(hostname) image=$XSH_IMAGE gpu=${_nv:-none} code=$(git -C "$XSH_ROOT" describe --tags --always --dirty)"
