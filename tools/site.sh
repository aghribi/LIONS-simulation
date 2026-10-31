#!/bin/bash
# site.sh: drive a compute site from the laptop clone. Sites: sites/<site>.env. Same layout as xsuite_heb v22.
#
#   tools/site.sh <site> deploy               clone or update the site copy of the current branch
#   tools/site.sh <site> setup                install the image of the current recipe (environment/), then envcheck
#   tools/site.sh <site> envcheck             image present, recipe hash as in the clone
#   tools/site.sh <site> submit [opts] jobs/<job>.sh [KEY=VALUE ...]   (options: tools/xsh-submit --help)
#   tools/site.sh <site> status               your jobs in the site queue
#   tools/site.sh <site> logs <jobid>         tail of the job's output and error logs
#   tools/site.sh <site> fetch <path>         copy <path> under the site results root to benchmarks/runs/sites/<site>/
#
# deploy and submit refuse a dirty or unpushed laptop branch; submit also refuses
# a site copy that is not at the laptop commit, so every job runs pushed code.
set -euo pipefail

die() { echo "site.sh: $*" >&2; exit 1; }
root="$(cd "$(dirname "$0")/.." && pwd)"
[ $# -ge 2 ] || { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
site="$1" cmd="$2"; shift 2
[ -f "$root/sites/$site.env" ] || die "unknown site '$site' (no sites/$site.env)"
source "$root/sites/$site.env"
[ -f "$root/sites/$site.local.env" ] && source "$root/sites/$site.local.env"

branch="$(git -C "$root" symbolic-ref --short HEAD)"
head="$(git -C "$root" rev-parse HEAD)"

# Anonymous https clone URL from the laptop's origin, and the image registry of the repository owner.
url="$(git -C "$root" remote get-url origin | sed -E 's#^git@([^:]+):#https://\1/#')"
owner="$(echo "$url" | sed -E 's#^https://github\.com/([^/]+)/.*#\1#' | tr '[:upper:]' '[:lower:]')"
registry="ghcr.io/$owner/lions-sim"
hash="$("$root/tools/env-hash")"
image="lions-sim-$hash.sif"

clean_and_pushed() {
  [ -z "$(git -C "$root" status --porcelain --untracked-files=no)" ] || die "uncommitted changes: commit first"
  git -C "$root" fetch -q origin "$branch"
  [ "$head" = "$(git -C "$root" rev-parse "origin/$branch")" ] || die "branch $branch is not pushed to origin: push first"
}

# Remote shell snippet: cd into the site copy of $branch, or fail.
cd_copy='dir=$(git -C "$CLONE" worktree list --porcelain 2>/dev/null | awk -v b="refs/heads/$BRANCH" '"'"'/^worktree /{w=$2} $0=="branch "b{print w}'"'"')
[ -n "$dir" ] || { echo "no site copy of $BRANCH: run tools/site.sh '"$site"' deploy" >&2; exit 3; }
cd "$dir"'

remote() {   # run a bash script on the site, in a login shell (site proxies and modules), with the site variables set
  ssh -o BatchMode=yes "$SSH" "CLONE='$CLONE' BRANCH='$branch' HEAD='$head' URL='$url' SITE='$site' bash -l -s" <<< "$1"
}

case "$cmd" in
deploy)
  clean_and_pushed
  remote 'set -e
if [ ! -d "$CLONE/.git" ]; then
  mkdir -p "$(dirname "$CLONE")"
  git clone -q --branch "$BRANCH" "$URL" "$CLONE"
fi
git -C "$CLONE" fetch -q origin
dir=$(git -C "$CLONE" worktree list --porcelain | awk -v b="refs/heads/$BRANCH" '"'"'/^worktree /{w=$2} $0=="branch "b{print w}'"'"')
if [ -z "$dir" ]; then
  dir="$CLONE.$BRANCH"
  git -C "$CLONE" worktree add -q "$dir" "$BRANCH"
fi
cd "$dir"
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "site copy $dir has local changes: not touching it" >&2; exit 4; }
git merge -q --ff-only "origin/$BRANCH"
[ "$(git rev-parse HEAD)" = "$HEAD" ] || { echo "site copy is at $(git rev-parse --short HEAD), laptop at ${HEAD:0:7}" >&2; exit 5; }
echo "$SITE" > sites/.site
echo "deployed $BRANCH $(git describe --tags --always) to $(hostname -s):$dir"'
  ;;
setup)
  if [ -n "${IMAGE_FROM:-}" ]; then   # no runtime on the login nodes: copy the image from another site
    src_ssh="$(source "$root/sites/$IMAGE_FROM.env"; echo "$SSH")"; src_dir="$(source "$root/sites/$IMAGE_FROM.env"; echo "$IMAGE_DIR")"
    remote "mkdir -p '$IMAGE_DIR' '$XSH_RESULTS'"
    ssh -o BatchMode=yes "$SSH" "test -f '$IMAGE_DIR/$image'" || scp -3 -q "$src_ssh:$src_dir/$image" "$SSH:$IMAGE_DIR/$image"
  else
    remote "set -e
${CONTAINER_SETUP:-}
mkdir -p '$IMAGE_DIR' '$XSH_RESULTS' ${IMAGE_STAGE:+'$IMAGE_STAGE'}
dest='${IMAGE_STAGE:-$IMAGE_DIR}/$image'
[ -f '$IMAGE_DIR/$image' ] || [ -f \"\$dest\" ] || $CONTAINER pull \"\$dest\" 'docker://$registry:$hash'
${IMAGE_INSTALL:+[ -f '$IMAGE_DIR/$image' ] || $IMAGE_INSTALL \"\$dest\"}"
  fi
  exec "$0" "$site" envcheck
  ;;
envcheck)
  # Sites whose login nodes cannot run containers: show the image; every job checks it before Python starts.
  [ "${COMPUTE_ONLY_RUNTIME:-}" ] && { echo "$site: containers run on compute nodes only; image $(ssh -o BatchMode=yes "$SSH" "ls -l '$IMAGE_DIR/$image'" 2>&1 | tail -1)"; exit 0; }
  remote "set -e
${CONTAINER_SETUP:-}
img='$IMAGE_DIR/$image'
[ -f \"\$img\" ] || { echo \"no image \$img: run tools/site.sh $site setup\" >&2; exit 3; }
echo \"image \$img, recipe \$($CONTAINER exec \"\$img\" cat /opt/lock/HASH) (code expects $hash), torch \$($CONTAINER exec \"\$img\" python -c 'import torch; print(torch.__version__)')\""
  ;;
submit)
  clean_and_pushed
  args=""; [ $# -gt 0 ] && args="$(printf '%q ' "$@")"
  remote "set -e
$cd_copy
[ \"\$(git rev-parse HEAD)\" = \"\$HEAD\" ] || { echo \"site copy is at \$(git rev-parse --short HEAD), laptop at \${HEAD:0:7}: run tools/site.sh $site deploy\" >&2; exit 5; }
XSH_SITE=$site tools/xsh-submit $args"
  ;;
status)
  remote 'squeue -u "$USER" -o "%.12i %.24j %.10P %.8T %.10M %.12l %R"'
  ;;
logs)
  [ $# -ge 1 ] || die "usage: tools/site.sh $site logs <jobid>"
  remote "$cd_copy
for f in jobs/logs/*_$1.out jobs/logs/*_$1.err; do [ -f \"\$f\" ] && { echo \"== \$f\"; tail -n 40 \"\$f\"; }; done; true"
  ;;
fetch)
  [ $# -ge 1 ] || die "usage: tools/site.sh $site fetch <path under $XSH_RESULTS>"
  case "$1" in /*|*..*) die "fetch takes a relative path without '..'" ;; esac
  mkdir -p "$root/benchmarks/runs/sites/$site/$(dirname "$1")"
  rsync -a "$SSH:$XSH_RESULTS/$1" "$root/benchmarks/runs/sites/$site/$(dirname "$1")/"
  echo "copied to $root/benchmarks/runs/sites/$site/$1"
  ;;
*) die "unknown command '$cmd'" ;;
esac
