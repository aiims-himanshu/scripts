#!/usr/bin/env bash
set -euo pipefail

# Rename fMRIPrep stroke tissue segment files so they match the healthy-control
# naming convention used by CONN:
#   *_label-GM_desc-preproc_probseg.nii.gz -> *_label-GM_probseg.nii.gz
#   *_label-WM_desc-preproc_probseg.nii.gz -> *_label-WM_probseg.nii.gz
#   *_label-CSF_desc-preproc_probseg.nii.gz -> *_label-CSF_probseg.nii.gz
#
# Usage:
#   bash /home/sskgroup/Git/scripts/renamne_segment_wgs.sh
# Optional: set DRY_RUN=1 inside the script before running.

die() {
  echo "ERROR: $*" >&2
  exit 1
}

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

DRY_RUN="${DRY_RUN:-0}"
BATCH_DIR="/media/sskgroup/Thesis/functional/fmriprep_24/batch_8"

[[ -d "$BATCH_DIR" ]] || die "Batch directory not found: $BATCH_DIR"

rename_one() {
  local src="$1"
  local dst="${src/_desc-preproc/}"

  if [[ "$src" == "$dst" ]]; then
    return 0
  fi

  if [[ -e "$dst" ]]; then
    log "SKIP exists: $dst"
    return 0
  fi

  if [[ $DRY_RUN -eq 1 ]]; then
    log "DRY-RUN mv: $src -> $dst"
  else
    mv -- "$src" "$dst"
    log "RENAMED: $src -> $dst"
  fi
}

total=0
for label in GM WM CSF; do
  while IFS= read -r -d '' file; do
    rename_one "$file"
    total=$((total + 1))
  done < <(
    find "$BATCH_DIR" -type f \
      \( -name "*_label-${label}_desc-preproc_probseg.nii.gz" -o -name "*_label-${label}_desc-preproc_probseg.nii" \) \
      -print0
  )
done

log "Processed ${total} tissue files under $BATCH_DIR"
