#!/bin/bash

# ------------------- Configuration -------------------

BIDS_ROOT="/media/sskgroup/Thesis/Lesion_BIDS/BIDS"
DATASETS=("Primary")

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=16
APPLY_N4=1

ANTS_REG="antsRegistration"
ANTS_APPLY="antsApplyTransforms"
N4CORRECT="N4BiasFieldCorrection"

TEMP_DIR="$BIDS_ROOT/temp"
mkdir -p "$TEMP_DIR"

# ------------------- Processing Loop -------------------

for dataset in "${DATASETS[@]}"; do
  echo "📂 Processing dataset: $dataset"
  for SUBJ_DIR in $BIDS_ROOT/$dataset/sub-*/anat; do
    subj=$(basename "$(dirname "$SUBJ_DIR")")
    echo "🔁 Processing $subj"

    T1="$SUBJ_DIR/${subj}_T1w.nii.gz"
    FLAIR="$SUBJ_DIR/${subj}_FLAIR.nii.gz"
    T2="$SUBJ_DIR/${subj}_T2w.nii.gz"
    TEMP_SUBJ="$TEMP_DIR/$subj"
    mkdir -p "$TEMP_SUBJ"

    [[ ! -f $T1 ]] && echo "❌ Missing T1w for $subj, skipping" && continue
    [[ ! -f $FLAIR && ! -f $T2 ]] && echo "❌ Missing FLAIR and T2w for $subj, skipping" && continue

    # ---------- FLAIR ----------
    if [[ -f $FLAIR ]]; then
      echo "🧠 Registering FLAIR to T1w (Affine only)..."
      $ANTS_REG -d 3 \
        -r [$T1,$FLAIR,1] \
        -m mattes[$T1,$FLAIR,1,32,regular,0.2] \
        -t affine[0.1] \
        -c [1000x500x250x0,1e-6,10] \
        -s 4x2x1x0vox \
        -f 6x4x2x1 \
        -o "$TEMP_SUBJ/${subj}_FLAIR_to_T1_affine_"

      echo "🎯 Resampling FLAIR to T1w..."
      $ANTS_APPLY -d 3 \
        -i "$FLAIR" \
        -r "$T1" \
        -t "$TEMP_SUBJ/${subj}_FLAIR_to_T1_affine_0GenericAffine.mat" \
        -n BSpline \
        -o "$SUBJ_DIR/${subj}_FLAIR_resampled.nii.gz"

      if [[ "$APPLY_N4" -eq 1 ]]; then
        echo "✨ N4BiasFieldCorrection on FLAIR..."
        $N4CORRECT -d 3 \
          -i "$SUBJ_DIR/${subj}_FLAIR_resampled.nii.gz" \
          -o "$SUBJ_DIR/${subj}_FLAIR_resampled_n4.nii.gz"
      fi
    fi

    # ---------- T2w ----------
    if [[ -f $T2 ]]; then
      echo "🧠 Registering T2w to T1w (Affine only)..."
      $ANTS_REG -d 3 \
        -r [$T1,$T2,1] \
        -m mattes[$T1,$T2,1,32,regular,0.2] \
        -t affine[0.1] \
        -c [1000x500x250x0,1e-6,10] \
        -s 4x2x1x0vox \
        -f 6x4x2x1 \
        -o "$TEMP_SUBJ/${subj}_T2w_to_T1_affine_"

      echo "🎯 Resampling T2w to T1w..."
      $ANTS_APPLY -d 3 \
        -i "$T2" \
        -r "$T1" \
        -t "$TEMP_SUBJ/${subj}_T2w_to_T1_affine_0GenericAffine.mat" \
        -n BSpline \
        -o "$SUBJ_DIR/${subj}_T2w_resampled.nii.gz"
    fi

    echo "✅ Finished $subj"
    echo "-------------------------------"
  done
done
