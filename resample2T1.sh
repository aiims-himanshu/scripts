#!/bin/bash

# ----------- Config -----------

BIDS_ROOT="/media/sskgroup/Thesis/Lesion_BIDS/BIDS"
DATASETS=("Primary" "Secondary")  # Process both

ITK_THREADS=4  # threads per subject
PARALLEL_SUBJECTS=4  # number of subjects in parallel

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$ITK_THREADS
APPLY_N4=1

ANTS_REG="antsRegistration"
ANTS_APPLY="antsApplyTransforms"
N4CORRECT="N4BiasFieldCorrection"

# ----------- Subject Function -----------

process_subject() {
    SUBJ_DIR="$1"
    subj=$(basename "$(dirname "$SUBJ_DIR")")

    echo "🔁 Processing $subj"
    T1="$SUBJ_DIR/${subj}_T1w.nii.gz"
    FLAIR="$SUBJ_DIR/${subj}_FLAIR.nii.gz"
    T2="$SUBJ_DIR/${subj}_T2w.nii.gz"

    [[ ! -f $T1 ]] && echo "❌ Missing T1w for $subj, skipping" && return
    [[ ! -f $FLAIR && ! -f $T2 ]] && echo "❌ Missing FLAIR and T2w for $subj, skipping" && return

    # ----- FLAIR -----
    if [[ -f $FLAIR ]]; then
        echo "🧠 Registering FLAIR to T1w for $subj..."
        $ANTS_REG -d 3 \
          -r [$T1,$FLAIR,1] \
          -m mattes[$T1,$FLAIR,1,32,regular,0.2] \
          -t affine[0.1] \
          -c [1000x500x250x0,1e-6,10] \
          -s 4x2x1x0vox \
          -f 6x4x2x1 \
          -m CC[$T1,$FLAIR,1,4] \
          -t SyN[0.1,3,0] \
          -c [100x70x50x0,1e-6,10] \
          -s 2x1x0x0vox \
          -f 4x2x1x1 \
          -o "$SUBJ_DIR/${subj}_FLAIR_to_T1_"

        echo "🎯 Resampling FLAIR to T1w..."
        $ANTS_APPLY -d 3 \
          -i $FLAIR \
          -r $T1 \
          -t "${SUBJ_DIR}/${subj}_FLAIR_to_T1_1Warp.nii.gz" \
          -t "${SUBJ_DIR}/${subj}_FLAIR_to_T1_0GenericAffine.mat" \
          -n BSpline \
          -o "${SUBJ_DIR}/${subj}_FLAIR_resampled.nii.gz"

        if [[ "$APPLY_N4" -eq 1 ]]; then
            echo "✨ Applying N4BiasFieldCorrection on FLAIR..."
            $N4CORRECT -d 3 \
              -i "${SUBJ_DIR}/${subj}_FLAIR_resampled.nii.gz" \
              -o "${SUBJ_DIR}/${subj}_FLAIR_resampled_n4.nii.gz"
        fi
    fi

    # ----- T2w -----
    if [[ -f $T2 ]]; then
        echo "🧠 Registering T2w to T1w for $subj..."
        $ANTS_REG -d 3 \
          -r [$T1,$T2,1] \
          -m mattes[$T1,$T2,1,32,regular,0.2] \
          -t affine[0.1] \
          -c [1000x500x250x0,1e-6,10] \
          -s 4x2x1x0vox \
          -f 6x4x2x1 \
          -m CC[$T1,$T2,1,4] \
          -t SyN[0.1,3,0] \
          -c [100x70x50x0,1e-6,10] \
          -s 2x1x0x0vox \
          -f 4x2x1x1 \
          -o "$SUBJ_DIR/${subj}_T2w_to_T1_"

        echo "🎯 Resampling T2w to T1w..."
        $ANTS_APPLY -d 3 \
          -i $T2 \
          -r $T1 \
          -t "${SUBJ_DIR}/${subj}_T2w_to_T1_1Warp.nii.gz" \
          -t "${SUBJ_DIR}/${subj}_T2w_to_T1_0GenericAffine.mat" \
          -n BSpline \
          -o "${SUBJ_DIR}/${subj}_T2w_resampled.nii.gz"
    fi

    echo "✅ Finished $subj"
}

export -f process_subject
export ANTS_REG ANTS_APPLY N4CORRECT APPLY_N4 ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

# ----------- Run Parallel on Primary + Secondary -----------

for dataset in "${DATASETS[@]}"; do
    echo "📂 Processing $dataset dataset..."
    find "$BIDS_ROOT/$dataset" -type d -name anat | sort | \
        xargs -n 1 -P $PARALLEL_SUBJECTS -I{} bash -c 'process_subject "$@"' _ {}
done
