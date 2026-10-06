#!/usr/bin/env bash
# Prepares the code that a pull request is compared with ("the base", normally the tip of
# main): a checkout of the base commit with the CI scripts of the current checkout, so
# that the same CI scripts and the same synthetic input are used for both.
#
#   .github/ci/prepare_base.sh BASE_COMMIT DIR
#
# Run from the checkout under test (the pull request). Afterwards:
#   DIR                the files of BASE_COMMIT (a detached git worktree of this checkout)
#   DIR/.github/ci     replaced by this checkout's .github/ci (the base may have none)
#   DIR/FCCAnalyses    if BASE_COMMIT pins the same FCCAnalyses commit as this checkout:
#                      a symbolic link to this checkout's FCCAnalyses, so one build serves
#                      both; otherwise the base's own FCCAnalyses, to be built separately.
# Prints the base commit and the FCCAnalyses decision; with GITHUB_OUTPUT set (in
# GitHub Actions) also writes sha=, short=, fccanalyses=same|different and fccanalyses_sha=.

set -euo pipefail
[[ $# -eq 2 ]] || { sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1; }
CI_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$CI_DIR/../.." && pwd)
BASE=$(git -C "$REPO" rev-parse --verify "$1^{commit}")
DIR=$2
[[ -e $DIR ]] && { echo "prepare_base: $DIR exists already"; exit 1; }

git -C "$REPO" worktree add --quiet --detach "$DIR" "$BASE"
DIR=$(cd "$DIR" && pwd)
# the worktree's link to the repository made relative, like the FCCAnalyses link below, so
# that git finds the repository wherever the checkout is mounted (the CI container)
echo "gitdir: $(realpath --relative-to="$DIR" "$(git -C "$DIR" rev-parse --absolute-git-dir)")" > "$DIR/.git"
rm -rf "$DIR/.github/ci"
mkdir -p "$DIR/.github"
cp -r "$CI_DIR" "$DIR/.github/ci"

head_sub=$(git -C "$REPO" rev-parse HEAD:FCCAnalyses)
base_sub=$(git -C "$REPO" rev-parse "$BASE:FCCAnalyses")
if [[ $base_sub == "$head_sub" ]]; then
  sub=same
  rmdir "$DIR/FCCAnalyses"
  ln -s ../FCCAnalyses "$DIR/FCCAnalyses"   # relative: valid wherever the checkout is mounted
else
  sub=different
  git -C "$DIR" submodule --quiet update --init --recursive --depth 1 FCCAnalyses \
    || git -C "$DIR" submodule --quiet update --init --recursive FCCAnalyses
fi

short=$(git -C "$REPO" rev-parse --short=7 "$BASE")
echo "base: $short ($(git -C "$REPO" log -1 --format=%s "$BASE")) in $DIR"
echo "FCCAnalyses of the base: ${base_sub:0:12} ($sub from this checkout's ${head_sub:0:12})"
if [[ -n ${GITHUB_OUTPUT:-} ]]; then
  {
    echo "sha=$BASE"
    echo "short=$short"
    echo "fccanalyses=$sub"
    echo "fccanalyses_sha=$base_sub"
  } >> "$GITHUB_OUTPUT"
fi
