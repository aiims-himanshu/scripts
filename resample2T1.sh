#!/bin/bash

# ------------------- Configuration -------------------

# Hardcoded base directory (adjust if needed)
BIDS_DIR="/media/sskgroup/Thesis/Lesion_BIDS/BIDS"
PRIMARY_DIR="/media/sskgroup/Thesis/Lesion_BIDS/BIDS/Primary"

# Enable or disable N4BiasFieldCorrection for FLAIR
APPLY_N4=1  # set to 0 to skip

# Max threads for ANTs
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$(nproc)

# ANTs path (optional if antsRegistration, antsApplyTransforms, N4BiasFieldCorrection in $PATH)
ANTS_REG="antsRegistration"
ANTS_APPLY="antsApplyTransforms"
N4CORRECT="N4BiasFieldCorrection"

# ------------------- Processing Loop -------------------

for SUBJ_DIR in $PRIMARY_DIR/sub-*/anat; do
    subj=$(basename "$(dirname "$SUBJ_DIR")")
    echo "🔁 Processing $subj"

    T1="$SUBJ_DIR/${subj}_T1w.nii.gz"
    FLAIR="$SUBJ_DIR/${subj}_FLAIR.nii.gz"
    T2="$SUBJ_DIR/${subj}_T2w.nii.gz"

    # Check existence
    [[ ! -f $T1 ]] && echo "❌ Missing T1w for $subj, skipping" && continue
    [[ ! -f $FLAIR && ! -f $T2 ]] && echo "❌ Missing FLAIR and T2w for $subj, skipping" && continue

    # ----- Register and Resample FLAIR -----
    if [[ -f $FLAIR ]]; then
        echo "🧠 Registering FLAIR to T1w..."

        $ANTS_REG -d 3 \
          -r [$T1,$FLAIR,1] \
          -m mattes[$T1,$FLAIR,1,32,regular,0.25] \
          -t affine[0.1] \
          -c [1000x500x250x0,1e-6,10] \
          -s 4x2x1x0vox \
          -f 6x4x2x1 \
          -o "$SUBJ_DIR/${subj}_FLAIR_to_T1_"

        echo "🎯 Resampling FLAIR to T1w space..."
        $ANTS_APPLY -d 3 \
          -i $FLAIR \
          -r $T1 \
          -t "${SUBJ_DIR}/${subj}_FLAIR_to_T1_0GenericAffine.mat" \
          -n BSpline \
          -o "${SUBJ_DIR}/${subj}_FLAIR_resampled.nii.gz"

        # Optional N4 correction
        if [[ "$APPLY_N4" -eq 1 ]]; then
            echo "✨ Applying N4BiasFieldCorrection to FLAIR..."
            $N4CORRECT -d 3 \
              -i "${SUBJ_DIR}/${subj}_FLAIR_resampled.nii.gz" \
              -o "${SUBJ_DIR}/${subj}_FLAIR_resampled_n4.nii.gz"
        fi
    fi

    # ----- Register and Resample T2w -----
    if [[ -f $T2 ]]; then
        echo "🧠 Registering T2w to T1w..."

        $ANTS_REG -d 3 \
          -r [$T1,$T2,1] \
          -m mattes[$T1,$T2,1,32,regular,0.25] \
          -t affine[0.1] \
          -c [1000x500x250x0,1e-6,10] \
          -s 4x2x1x0vox \
          -f 6x4x2x1 \
          -o "$SUBJ_DIR/${subj}_T2w_to_T1_"

        echo "🎯 Resampling T2w to T1w space..."
        $ANTS_APPLY -d 3 \
          -i $T2 \
          -r $T1 \
          -t "${SUBJ_DIR}/${subj}_T2w_to_T1_0GenericAffine.mat" \
          -n BSpline \
          -o "${SUBJ_DIR}/${subj}_T2w_resampled.nii.gz"
    fi

    echo "✅ Finished $subj"
    echo "-------------------------------"
done
