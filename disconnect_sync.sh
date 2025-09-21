#!/usr/bin/env bash
set -euo pipefail


#usage note bash ./disconnect_sync.sh
# ====== USER CONFIGURATION ======
OUTROOT="/home/speech/Desktop/Primary_MNI"
TMPROOT="/home/speech/Desktop/Primary_MNI_TMP"
BCBTOOLS="/home/speech/Desktop/BCBToolKit/Tools"
MNI_TEMPLATE="${BCBTOOLS}/extraFiles/MNI152.nii.gz"
THREADS_PER_SYN=4
NTHREADS=44
MAX_PARALLEL_SYN=$((NTHREADS/THREADS_PER_SYN))

export FSLOUTPUTTYPE="NIFTI_GZ"

logmsg() { echo "[$(date +"%F %T")] $*"; }

logmsg "RESUME: Starting SyN registration from TMP for all preprocessed subjects."

# Find all TMP folders with finished preprocessing
SUBJECTS=()
for TMPD in "$TMPROOT"/sub-*; do
    [[ -d "$TMPD" ]] || continue
    # Only add subjects with t1_brain_final.txt
    [[ -f "$TMPD/t1_brain_final.txt" ]] && SUBJECTS+=("$(basename "$TMPD")")
done

#########################################
# 1. SyN Registration (Parallel)        #
#########################################
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$THREADS_PER_SYN
export OMP_NUM_THREADS=$THREADS_PER_SYN

syn_registration() {
    SUB="$1"
    TMPD="$TMPROOT/$SUB"
    OD="$OUTROOT/$SUB/norm_final"
    MNI_TEMPLATE="$2"

    METHOD=$(cat "$TMPD/method.txt")
    T1_BRAIN=$(cat "$TMPD/t1_brain_final.txt")
    COSTMASK=$(cat "$TMPD/costmask.txt")
    ANTSOUT="$TMPD/${SUB}_xfm_"
    
    # Skip if already done
    if [[ -f "${ANTSOUT}1Warp.nii.gz" && -f "${ANTSOUT}0GenericAffine.mat" ]]; then
        echo "[$SUB] SyN outputs already present, skipping..."
        return
    fi

    ANTS_CMD="antsRegistration -d 3 --float 1 \
        --output [$ANTSOUT,${ANTSOUT}Warped.nii.gz,${ANTSOUT}InverseWarped.nii.gz] \
        --interpolation Linear \
        --winsorize-image-intensities [0.005,0.995] \
        --use-histogram-matching 1 \
        --initial-moving-transform [$MNI_TEMPLATE,$T1_BRAIN,1] \
        --transform Rigid[0.1] \
        --metric MI[$MNI_TEMPLATE,$T1_BRAIN,1,32,Regular,0.25] \
        --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 \
        --smoothing-sigmas 3x2x1x0vox \
        --transform Affine[0.1] \
        --metric MI[$MNI_TEMPLATE,$T1_BRAIN,1,32,Regular,0.25] \
        --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 8x4x2x1 \
        --smoothing-sigmas 3x2x1x0vox \
        --transform SyN[0.1,3,0] \
        --metric CC[$MNI_TEMPLATE,$T1_BRAIN,1,4] \
        --convergence [100x70x50x20,1e-6,10] \
        --shrink-factors 8x4x2x1 \
        --smoothing-sigmas 3x2x1x0vox"
    if [[ "$METHOD" == "CLASSIC" ]]; then
        ANTS_CMD="$ANTS_CMD --masks [$COSTMASK,NULL]"
    fi
    mkdir -p "$OD"
    logmsg "[$SUB] Starting ANTs registration..." | tee -a "$OD/ants_registration.log"
    eval $ANTS_CMD | tee -a "$OD/ants_registration.log"
    logmsg "[$SUB] Finished ANTs registration." | tee -a "$OD/ants_registration.log"
}
export -f syn_registration
export TMPROOT
export OUTROOT
export MNI_TEMPLATE

parallel -j $MAX_PARALLEL_SYN syn_registration ::: "${SUBJECTS[@]}" ::: "$MNI_TEMPLATE"

#########################################
# 2. Post-SyN Loop (Apply, QC, Output)  #
#########################################
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NTHREADS
export OMP_NUM_THREADS=$NTHREADS

logmsg "Step 2: Post-SyN (apply transforms, QC, etc.)"
for SUB in "${SUBJECTS[@]}"; do
    TMPD="$TMPROOT/$SUB"
    OD="$OUTROOT/$SUB/norm_final"
    mkdir -p "$OD" "$OD/qc"
    MNI_TEMPLATE="$MNI_TEMPLATE"

    T1_BRAIN=$(cat "$TMPD/t1_brain_final.txt")
    LES_REO="$TMPD/${SUB}_lesion_reo.nii.gz"
    ANTSOUT="$TMPD/${SUB}_xfm_"

    # Transform T1 brain
    if [[ ! -f "$OD/${SUB}_T1_MNI.nii.gz" ]]; then
        antsApplyTransforms -d 3 \
            -i "$T1_BRAIN" \
            -r "$MNI_TEMPLATE" \
            -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
            -o "$OD/${SUB}_T1_MNI.nii.gz" \
            -n Linear
        fslcpgeom "$MNI_TEMPLATE" "$OD/${SUB}_T1_MNI.nii.gz"
    fi

    # Transform lesion mask
    if [[ ! -f "$OD/${SUB}_lesion_MNI.nii.gz" ]]; then
        antsApplyTransforms -d 3 \
            -i "$LES_REO" \
            -r "$MNI_TEMPLATE" \
            -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
            -o "$OD/${SUB}_lesion_MNI.nii.gz" \
            -n NearestNeighbor
        fslcpgeom "$MNI_TEMPLATE" "$OD/${SUB}_lesion_MNI.nii.gz"
    fi

    # QC PNG
    if [[ ! -f "$OD/qc/${SUB}_lesion_on_mni.png" ]]; then
        overlay 1 1 "$MNI_TEMPLATE" -a "$OD/${SUB}_lesion_MNI.nii.gz" 0.5 1.0 "$OD/qc/${SUB}_lesion_on_mni.nii.gz"
        slices "$OD/qc/${SUB}_lesion_on_mni.nii.gz" -o "$OD/qc/${SUB}_lesion_on_mni.png"
        rm -f "$OD/qc/${SUB}_lesion_on_mni.nii.gz"
    fi

    logmsg "[OK] [$SUB] Lesion/T1 normalized to MNI and QC done." | tee -a "$OUTROOT/logs/normalize.log"
done

logmsg "All subjects post-SyN steps completed! Outputs in $OUTROOT and intermediates in $TMPROOT."

