#!/usr/bin/env bash
ATLAS_DIR=/home/sskgroup/Documents/atlases
OUT_DIR=/home/sskgroup/Documents/atlases/custom_qsirecon
REF_TPL=/home/sskgroup/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz  # or _res-01_

mkdir -p "$OUT_DIR"

for atlas in "$ATLAS_DIR"/*.nii.gz; do
  base=$(basename "$atlas")
  out="$OUT_DIR/${base%.nii.gz}_MNI2009cAsym.nii.gz"

  antsApplyTransforms -d 3 -i "$atlas" -r "$REF_TPL" \
    -n NearestNeighbor -t identity -o "$out"

  c3d "$out" -orient LPS -o "$out"
  nifti_tool -mod_hdr -mod_field sform_code 0 -infiles "$out" -overwrite
  fslorient -setsform 0 0 0 0 0 0 0 0 0 0 0 0 "$out"

  echo "✓ $base → $(basename "$out")"
done
