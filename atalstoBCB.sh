#!/usr/bin/env bash
set -euo pipefail

# ============================
# HARD-CODED PATHS (EDIT)
# ============================
REF="/home/sskgroup/Documents/BCBToolKit/Tools/extraFiles/MNI152.nii.gz"
BRAINMASK="/home/sskgroup/Documents/BCBToolKit/Tools/extraFiles/MNI152_for_antsBrainExtractionMask.nii.gz"

# Raw atlas sources to resample INTO 1mm
# Uncomment if you actually have meta maps to process:
# ATLAS_META_RAW="/media/sskgroup/fast/modelling/atlas_meta"
ATLAS_DISC85_RAW="/media/sskgroup/fast/modelling/atlas_disconnectome"
ATLAS_COMP_RAW="/media/sskgroup/fast/modelling/Stroke_disconnection_components"

OUTROOT="/media/sskgroup/fast/modelling/atlas"
mkdir -p "$OUTROOT/ATLAS_meta" "$OUTROOT/ATLAS_disconnectome" "$OUTROOT/ATLAS_components"

# ============================
# Helpers
# ============================
res_lin_mask () {
  in="$1"; out="$2"
  # Skip if output exists and is newer than input
  if [ -e "$out" ] && [ "$out" -nt "$in" ]; then
    echo "  [skip] $(basename "$out") (up-to-date)"
    return
  fi
  antsApplyTransforms -d 3 -i "$in" -r "$REF" -o "$out" -n Linear -t identity
  fslmaths "$out" -mas "$BRAINMASK" "$out"
}

process_dir () {
  in_dir="$1"; out_sub="$2"
  if [ ! -d "$in_dir" ]; then
    echo "  [warn] missing: $in_dir — skipping"
    return
  fi
  out_dir="$OUTROOT/$out_sub"
  mkdir -p "$out_dir"

  # Find .nii / .nii.gz that are NOT already *_res-1mm.*
  i=0
  total=$(find "$in_dir" -maxdepth 1 \( -name '*.nii' -o -name '*.nii.gz' \) \
          ! -name '*_res-1mm.nii' ! -name '*_res-1mm.nii.gz' | wc -l | awk '{print $1}')
  if [ "$total" -eq 0 ]; then
    echo "  [info] no inputs in $in_dir — skipping"
    return
  fi

  # Use -print0 to be safe with spaces
  find "$in_dir" -maxdepth 1 \( -name '*.nii' -o -name '*.nii.gz' \) \
       ! -name '*_res-1mm.nii' ! -name '*_res-1mm.nii.gz' -print0 |
  while IFS= read -r -d '' f; do
    i=$((i+1))
    base="$(basename "$f")"
    base="${base%.nii}"; base="${base%.gz}"
    out="$out_dir/${base}_res-1mm.nii.gz"
    echo "  [$i/$total] $base → $(basename "$out")"
    res_lin_mask "$f" "$out"
  done
}

# ============================
# Run
# ============================

echo "[1/3] ATLAS_meta → 1mm (if configured)"
if [ "${ATLAS_META_RAW-}" ]; then
  process_dir "$ATLAS_META_RAW" "ATLAS_meta"
else
  echo "  [info] ATLAS_META_RAW not set — skipping meta"
fi

echo "[2/3] ATLAS_disconnectome_85 → 1mm"
process_dir "$ATLAS_DISC85_RAW" "ATLAS_disconnectome"

echo "[3/3] Stroke+disconnectome+components → 1mm"
process_dir "$ATLAS_COMP_RAW" "ATLAS_components"

echo "[OK] Resampled atlases in: $OUTROOT"
echo "Counts:"
echo "  meta:        $(ls -1 "$OUTROOT/ATLAS_meta"/*_res-1mm.nii.gz 2>/dev/null | wc -l || true)"
echo "  disconnect.: $(ls -1 "$OUTROOT/ATLAS_disconnectome"/*_res-1mm.nii.gz 2>/dev/null | wc -l || true)"
echo "  components:  $(ls -1 "$OUTROOT/ATLAS_components"/*_res-1mm.nii.gz 2>/dev/null | wc -l || true)"