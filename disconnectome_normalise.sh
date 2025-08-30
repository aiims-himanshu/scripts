#!/usr/bin/env bash
set -euo pipefail

# ╭─────────────────────────────────────────────────────────────╮
# │     BCB-Style Lesion/T1 Normalization: Loop–Parallel–Loop  │
# ╰─────────────────────────────────────────────────────────────╯
# run command ====   bash ./disconnectome_normalise.sh
# ====== USER CONFIGURATION ======
INROOT="/media/sskgroup/fast/Lesion/batch1"
OUTROOT="/media/sskgroup/fast/Lesion/Primary_MNI"
TMPROOT="/media/sskgroup/fast/Lesion/Primary_MNI_TMP"
BCBTOOLS="/home/sskgroup/Documents/BCBToolKit/Tools"

MNI_TEMPLATE="${BCBTOOLS}/extraFiles/MNI152.nii.gz"
MNI_TEMPLATE_WSKULL="${BCBTOOLS}/extraFiles/MNI152_wskull.nii.gz"
BRAIN_EXTRACTION_PRIOR="${BCBTOOLS}/extraFiles/Priors/brainPrior.nii.gz"
BRAIN_EXTRACTION_TEMPLATE="${BCBTOOLS}/extraFiles/Priors/brainWithSkullTemplate.nii.gz"

NTHREADS=16                # System total
THREADS_PER_SYN=2          # Threads per SyN process (for ITK/OMP)
MAX_PARALLEL_SYN=$((NTHREADS/THREADS_PER_SYN))
LESION_VOX_THRESHOLD=5000

export FSLOUTPUTTYPE="NIFTI_GZ"
# ===============================

mkdir -p "$OUTROOT/logs"
mkdir -p "$TMPROOT"

SUBJECTS=($(cd "$INROOT"; ls sub-*_*T1w.nii.gz | sed 's/_T1w.nii.gz//'))

logmsg() { echo "[$(date +"%F %T")] $*"; }

#####################################
# 1. Preprocess Loop (Multi-thread) #
#####################################

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NTHREADS
export OMP_NUM_THREADS=$NTHREADS

logmsg "Step 1: Preprocessing all subjects (serial, max threads)"
for SUB in "${SUBJECTS[@]}"; do
    T1="$INROOT/${SUB}_T1w.nii.gz"
    LES="$INROOT/${SUB}_label-lesion_roi.nii.gz"
    if [[ ! -f "$LES" ]]; then
        logmsg "[WARN] Skipping $SUB: No lesion mask found." | tee -a "$OUTROOT/logs/normalize.log"
        continue
    fi
    OD="$OUTROOT/$SUB/norm_final"
    TMPD="$TMPROOT/$SUB"
    mkdir -p "$OD/qc" "$TMPD"
    
    # Reorient inputs, save to TMP
    T1_REO="$TMPD/${SUB}_T1w_reo.nii.gz"
    LES_REO="$TMPD/${SUB}_lesion_reo.nii.gz"
    fslreorient2std "$T1" "$T1_REO"
    fslreorient2std "$LES" "$LES_REO"

    # Lesion size & basic bilateral check (adjust if needed for your data)
    nvox=$(fslstats "$LES_REO" -V | awk '{print $1+0}')
    right_lesion=$(fslmaths "$LES_REO" -roi 0 -1 0 -1 0 -1 0 -40 "$TMPD/right_lesion.nii.gz"; fslstats "$TMPD/right_lesion.nii.gz" -V | awk '{print $1+0}')
    left_lesion=$(fslmaths "$LES_REO" -roi 40 -1 0 -1 0 -1 0 -1 "$TMPD/left_lesion.nii.gz"; fslstats "$TMPD/left_lesion.nii.gz" -V | awk '{print $1+0}')
    
    # Choose method: Classic (large, bilateral) or Enantiomorphic (else)
    if [[ "$nvox" -gt $LESION_VOX_THRESHOLD && "$right_lesion" -gt 100 && "$left_lesion" -gt 100 ]]; then
        echo "CLASSIC" > "$TMPD/method.txt"
        logmsg "[$SUB] Classic method (large, bilateral: $nvox vox)" | tee -a "$OUTROOT/logs/normalize.log"
        # Classic cost-function mask: all brain except lesion
        COSTMASK="$TMPD/${SUB}_costmask.nii.gz"
        fslmaths "$LES_REO" -bin -mul -1 -add 1 "$COSTMASK"
        echo "$COSTMASK" > "$TMPD/costmask.txt"
        echo "$T1_REO" > "$TMPD/t1reg.txt"
    else
        echo "ENANTIO" > "$TMPD/method.txt"
        logmsg "[$SUB] Enantiomorphic method (small/unilateral: $nvox vox)" | tee -a "$OUTROOT/logs/normalize.log"
        # Enantiomorphic fill
        BASE=$(basename "$T1_REO" .nii.gz)
        # --- Inline enantiomorphic fill (as BCB) ---
        flirt -in "$T1_REO" -ref "$MNI_TEMPLATE_WSKULL" -omat "$TMPD/affine.mat" -out "$TMPD/output${BASE}.nii.gz"
        flirt -in "$LES_REO" -ref "$MNI_TEMPLATE_WSKULL" -applyxfm -init "$TMPD/affine.mat" -out "$TMPD/affineLesion${BASE}.nii.gz"
        fslswapdim "$TMPD/affineLesion${BASE}.nii.gz" -x y z "$TMPD/flippedaffine${BASE}.nii.gz"
        fslmaths "$TMPD/output${BASE}.nii.gz" -mas "$TMPD/flippedaffine${BASE}.nii.gz" "$TMPD/healthytissue${BASE}.nii.gz"
        fslswapdim "$TMPD/healthytissue${BASE}.nii.gz" -x y z "$TMPD/flippedhealthytissue${BASE}.nii.gz"
        convert_xfm -omat "$TMPD/inverseAffine.mat" -inverse "$TMPD/affine.mat"
        flirt -in "$TMPD/flippedhealthytissue${BASE}.nii.gz" -ref "$T1_REO" -applyxfm -init "$TMPD/inverseAffine.mat" -out "$TMPD/nativeflippedhealthytissue${BASE}.nii.gz"
        fslmaths "$TMPD/nativeflippedhealthytissue${BASE}.nii.gz" -mas "$LES_REO" "$TMPD/mnativeflippedhealthytissue${BASE}.nii.gz"
        fslmaths "$LES_REO" -add 1 -uthr 1 -bin "$TMPD/lesionedMask.nii.gz"
        fslmaths "$T1_REO" -mul "$TMPD/lesionedMask.nii.gz" "$TMPD/T1pitted.nii.gz"
        fslmaths "$TMPD/T1pitted.nii.gz" -add "$TMPD/mnativeflippedhealthytissue${BASE}.nii.gz" "$TMPD/Enantiomorphic${BASE}.nii.gz"
        ENAN_T1="$TMPD/Enantiomorphic${BASE}.nii.gz"
        echo "$ENAN_T1" > "$TMPD/t1reg.txt"
        echo "" > "$TMPD/costmask.txt"
    fi

    # Brain Extraction (ANTs)
    T1REG=$(cat "$TMPD/t1reg.txt")
    T1_BRAIN="$TMPD/${SUB}_T1_brain.nii.gz"
    antsBrainExtraction.sh -d 3 \
        -a "$T1REG" \
        -e "$BRAIN_EXTRACTION_TEMPLATE" \
        -m "$BRAIN_EXTRACTION_PRIOR" \
        -o "$TMPD/${SUB}_be_"
    mv "$TMPD/${SUB}_be_BrainExtractionBrain.nii.gz" "$T1_BRAIN"
    echo "$T1_BRAIN" > "$TMPD/t1_brain_final.txt"
    # Costmask/costfunc passed via txt file for later use
done

#########################################
# 2. Parallel SyN Registration (ANTS)   #
#########################################
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$THREADS_PER_SYN
export OMP_NUM_THREADS=$THREADS_PER_SYN

logmsg "Step 2: Parallel SyN registration for all subjects (${MAX_PARALLEL_SYN} jobs @ ${THREADS_PER_SYN} threads each)"

syn_registration() {
    SUB="$1"
    TMPD="$TMPROOT/$SUB"
    OD="$OUTROOT/$SUB/norm_final"
    MNI_TEMPLATE="$2"

    METHOD=$(cat "$TMPD/method.txt")
    T1_BRAIN=$(cat "$TMPD/t1_brain_final.txt")
    COSTMASK=$(cat "$TMPD/costmask.txt")
    ANTSOUT="$TMPD/${SUB}_xfm_"
    
    # Build ANTs registration command (BCB style)
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
# 3. Post-SyN Loop (Multi-thread again) #
#########################################
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NTHREADS
export OMP_NUM_THREADS=$NTHREADS

logmsg "Step 3: Post-SyN (apply transforms, QC, etc.)"
for SUB in "${SUBJECTS[@]}"; do
    TMPD="$TMPROOT/$SUB"
    OD="$OUTROOT/$SUB/norm_final"
    mkdir -p "$OD" "$OD/qc"
    MNI_TEMPLATE="$MNI_TEMPLATE"

    T1_BRAIN=$(cat "$TMPD/t1_brain_final.txt")
    LES_REO="$TMPD/${SUB}_lesion_reo.nii.gz"
    ANTSOUT="$TMPD/${SUB}_xfm_"

    # Transform T1 brain
    antsApplyTransforms -d 3 \
        -i "$T1_BRAIN" \
        -r "$MNI_TEMPLATE" \
        -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
        -o "$OD/${SUB}_T1_MNI.nii.gz" \
        -n Linear
    fslcpgeom "$MNI_TEMPLATE" "$OD/${SUB}_T1_MNI.nii.gz"

    # Transform lesion mask
    antsApplyTransforms -d 3 \
        -i "$LES_REO" \
        -r "$MNI_TEMPLATE" \
        -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
        -o "$OD/${SUB}_lesion_MNI.nii.gz" \
        -n NearestNeighbor
    fslcpgeom "$MNI_TEMPLATE" "$OD/${SUB}_lesion_MNI.nii.gz"

    # QC PNG
    overlay 1 1 "$MNI_TEMPLATE" -a "$OD/${SUB}_lesion_MNI.nii.gz" 0.5 1.0 "$OD/qc/${SUB}_lesion_on_mni.nii.gz"
    slices "$OD/qc/${SUB}_lesion_on_mni.nii.gz" -o "$OD/qc/${SUB}_lesion_on_mni.png"
    rm -f "$OD/qc/${SUB}_lesion_on_mni.nii.gz"

    logmsg "[OK] [$SUB] Lesion/T1 normalized to MNI and QC done." | tee -a "$OUTROOT/logs/normalize.log"
done

logmsg "All subjects processed! All intermediates are in $TMPROOT for backup/restore."

