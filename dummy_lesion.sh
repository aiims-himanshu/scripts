#!/bin/bash

# === CONFIGURATION ===
BIDS_ROOT="/media/sskgroup/Thesis/stroke_test/input"

echo "🔍 Scanning BIDS directory for control subjects..."

for subj_dir in "${BIDS_ROOT}"/sub-*; do
  subj=$(basename "$subj_dir")

  T1W_IMG="${subj_dir}/anat/${subj}_T1w.nii.gz"
  LESION_MASK="${subj_dir}/anat/${subj}_label-lesion_roi.nii.gz"

  if [[ -f "$LESION_MASK" ]]; then
    echo "✅ $subj already has lesion mask — skipping"
    continue
  fi

  if [[ ! -f "$T1W_IMG" ]]; then
    echo "❌ Missing T1w image for $subj — skipping"
    continue
  fi

  echo "🧠 Creating visible dummy lesion mask for $subj..."

  # Step 1: Zero-filled template
  temp_zero="${subj_dir}/anat/${subj}_temp_zero.nii.gz"
  fslmaths "$T1W_IMG" -mul 0 "$temp_zero"

  # Step 2: Add a small visible cube in deep white matter (e.g., x=40:50, y=60:70, z=30:40)
  # Adjust ROI to fit in most brains with standard orientation
  fslmaths "$temp_zero" -roi 40 10 60 10 30 10 0 1 "$LESION_MASK"

  # Step 3: Clean up
  rm "$temp_zero"

  if [[ -f "$LESION_MASK" ]]; then
    echo "✅ Dummy lesion mask created: $LESION_MASK"
  else
    echo "⚠️ Failed to create dummy lesion mask for $subj"
  fi
done

echo "🏁 All dummy lesions created."
