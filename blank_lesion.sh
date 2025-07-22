#!/bin/bash

# === CONFIGURATION ===
BIDS_ROOT="/media/sskgroup/Thesis/stroke_test/input"  # e.g., /home/sskgroup/Documents/stroke_test

# === SCRIPT START ===
echo "🔍 Scanning BIDS directory for control subjects..."

for subj_dir in "${BIDS_ROOT}"/sub-*; do
  subj=$(basename "$subj_dir")

  T1W_IMG="${subj_dir}/anat/${subj}_T1w.nii.gz"
  LESION_MASK="${subj_dir}/anat/${subj}_label-lesion_roi.nii.gz"

  # Check if lesion mask already exists
  if [[ -f "$LESION_MASK" ]]; then
    echo "✅ $subj already has lesion mask — skipping"
    continue
  fi

  # Check if T1w image exists
  if [[ ! -f "$T1W_IMG" ]]; then
    echo "❌ Missing T1w image for $subj — skipping"
    continue
  fi

  echo "🧠 Creating blank lesion mask for $subj in T1w space..."
  fslmaths "$T1W_IMG" -mul 0 "$LESION_MASK"

  if [[ -f "$LESION_MASK" ]]; then
    echo "✅ Blank lesion mask created: $LESION_MASK"
  else
    echo "⚠️ Failed to create lesion mask for $subj"
  fi
done

echo "🏁 Done."
