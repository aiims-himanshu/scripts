#!/usr/bin/env bash
# -------------------------------------------------------------
# get_full_templateflow_no_symlinks.sh
# Usage: ./get_full_templateflow_no_symlinks.sh  /ABS/PATH/TO/TEMPLATEFLOW  [N_JOBS]
# -------------------------------------------------------------
set -euo pipefail

TF_DIR=${1:-"$HOME/.cache/templateflow"}    # destination (default: ~/.cache/…)
N_JOBS=${2:-4}                              # parallel jobs for datalad get

echo ">>> Target directory: $TF_DIR"
mkdir -p "$TF_DIR"

# ---------------------------------------------------- 0) make sure tools exist
if ! command -v datalad >/dev/null 2>&1; then
  echo ">>> Installing DataLad + git-annex (Debian/Ubuntu)…"
  sudo apt-get update
  sudo apt-get install -y git git-annex datalad
fi

# ---------------------------------------------------- 1) environment
export TEMPLATEFLOW_HOME="$TF_DIR"
export GIT_ANNEX_FULL_HASH=1                # safer hashing on some filesystems

# ---------------------------------------------------- 2) clone TemplateFlow super-dataset
echo ">>> Cloning TemplateFlow super-dataset (recursive)…"
cd "$(dirname "$TF_DIR")"
datalad install -r -g \
  https://github.com/templateflow/templateflow.git \
  "$(basename "$TF_DIR")"                   # -r recurse, -g get data

# ---------------------------------------------------- 3) fetch every annexed object
echo ">>> Downloading all files with $N_JOBS parallel jobs…"
cd "$TF_DIR"
datalad get -r -J "$N_JOBS" .

# ---------------------------------------------------- 4) replace symlinks with real files
echo ">>> Converting git-annex symlinks → actual files (unlock)…"
git annex unlock . --recursive

# ---------------------------------------------------- 5) harden for concurrent read-only use
git config --local annex.merge-annex-branches false
find "$TF_DIR" -maxdepth 2 -type d -name '.git' \
  -execdir git config --local annex.merge-annex-branches false \;

echo "✅  TemplateFlow ready at $TF_DIR with NO symlinks."
echo "   Add 'export TEMPLATEFLOW_HOME=$TF_DIR' to your ~/.bashrc"

# Hint: on filesystems that cannot hard-link into the worktree, the unlock step
# will copy data, roughly doubling disk usage (≈ 45 GB → ≈ 90 GB).
# Check 'du -sh $TF_DIR' afterwards and plan storage accordingly.
