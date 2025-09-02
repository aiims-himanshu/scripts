#!/usr/bin/env bash
set -euo pipefail

# ╭─────────────────────────────────────────────────────────────╮
# │  normalize_patients_and_controls_to_MNI.sh                  │
# │  Unified BCB-style normalization for PATIENTS & CONTROLS    │
# │  - Patients: lesion-aware (cost-function mask / enantiomorphic)
# │  - Controls: standard (no lesion), robust ANTs pipeline     │
# │  Parallel SyN with fixed threads per job (e.g., 4 cores)    │
# ╰─────────────────────────────────────────────────────────────╯
# USAGE:
#   bash ./normalize_patients_and_controls_to_MNI.sh
#
# BEFORE YOU RUN:
#   1) Edit the USER CONFIGURATION section below.
#   2) Ensure FSL, ANTs, GNU parallel are installed & on PATH.
#   3) Input naming: expects files like sub-XXX_T1w.nii.gz .
#      For PATIENTS it expects sub-XXX_label-lesion_roi.nii.gz next to T1.
#
# OUTPUT LAYOUT:
#   ${OUTROOT}/patients/<SUB>/norm_final/
#   ${OUTROOT}/controls/<SUB>/norm_final/
#   with QC PNGs under .../qc/
#
# THREADING MODEL:
#   - System cores: NTHREADS (e.g., 32)
#   - Each SyN job: THREADS_PER_SYN (e.g., 4)
#   - Max parallel SyN jobs = NTHREADS / THREADS_PER_SYN
#   - ITK/OMP threads are set appropriately per stage.

# ======================= USER CONFIGURATION =======================
# Patients (with lesion masks)
PATIENT_INROOT="/media/sskgroup/fast/Lesion/batch1"
# Controls (T1 only). Point this to the folder that contains T1s (e.g., /path/Control/T1)
CONTROL_INROOT="/media/sskgroup/fast/Controls/T1"

# Output & temp roots (separate parents for clarity)
OUTROOT="/media/sskgroup/fast/Lesion/Primary_MNI"
TMPROOT="/media/sskgroup/fast/Lesion/Primary_MNI_TMP"

# BCBToolKit resources
BCBTOOLS="/home/sskgroup/Documents/BCBToolKit/Tools"
MNI_TEMPLATE="${BCBTOOLS}/extraFiles/MNI152.nii.gz"
MNI_TEMPLATE_WSKULL="${BCBTOOLS}/extraFiles/MNI152_wskull.nii.gz"
BRAIN_EXTRACTION_PRIOR="${BCBTOOLS}/extraFiles/Priors/brainPrior.nii.gz"
BRAIN_EXTRACTION_TEMPLATE="${BCBTOOLS}/extraFiles/Priors/brainWithSkullTemplate.nii.gz"

# Threads (32-core system; 4 cores per SyN job)
NTHREADS=${NTHREADS:-32}
THREADS_PER_SYN=${THREADS_PER_SYN:-4}
MAX_PARALLEL_SYN=$(( NTHREADS / THREADS_PER_SYN ))

# Lesion logic
LESION_VOX_THRESHOLD=${LESION_VOX_THRESHOLD:-5000}

# Filename patterns (adjust only if your naming differs)
T1_PATTERN="sub-*_*T1w.nii.gz"
LESION_SUFFIX="label-lesion_roi.nii.gz"

# FSL NIfTI gz output by default
export FSLOUTPUTTYPE="NIFTI_GZ"
# =================================================================

# --------------- sanity checks ---------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "[FATAL] Missing: $1"; exit 1; }; }
need fslmaths; need fslstats; need fslreorient2std; need fslcpgeom
need antsRegistration; need antsApplyTransforms; need antsBrainExtraction.sh
need flirt; need fslswapdim; need convert_xfm; need overlay; need slices
need parallel

mkdir -p "${OUTROOT}/patients" "${OUTROOT}/controls" "${OUTROOT}/logs"
mkdir -p "${TMPROOT}/patients" "${TMPROOT}/controls"

logmsg() { echo "[$(date +"%F %T")] $*"; }

# --------------- subject discovery ---------------
collect_subjects() {
  local root="$1"; local kind="$2"  # kind ∈ {patients,controls}
  local -a subs=()
  if [[ -d "$root" ]]; then
    while IFS= read -r -d '' t1; do
      local base
      base=$(basename "$t1")
      base=${base%_T1w.nii.gz}
      subs+=("$base")
    done < <(find "$root" -maxdepth 1 -type f -name "$T1_PATTERN" -print0 | sort -z)
  fi
  echo "${subs[@]:-}"
}

PATIENT_SUBJECTS=( $(collect_subjects "$PATIENT_INROOT" patients) )
CONTROL_SUBJECTS=( $(collect_subjects "$CONTROL_INROOT" controls) )

logmsg "Found ${#PATIENT_SUBJECTS[@]} patient T1s in $PATIENT_INROOT"
logmsg "Found ${#CONTROL_SUBJECTS[@]} control T1s in $CONTROL_INROOT"

# --------------- threading presets ---------------
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NTHREADS
export OMP_NUM_THREADS=$NTHREADS

# --------------- preprocessing: patients ---------------
logmsg "Step 1A: Preprocess PATIENTS (serial, max threads=$NTHREADS)"
for SUB in "${PATIENT_SUBJECTS[@]:-}"; do
  T1="${PATIENT_INROOT}/${SUB}_T1w.nii.gz"
  LES="${PATIENT_INROOT}/${SUB}_${LESION_SUFFIX}"
  if [[ ! -f "$LES" ]]; then
    logmsg "[WARN][PAT] Skipping $SUB: No lesion mask found." | tee -a "${OUTROOT}/logs/normalize_patients.log"
    continue
  fi
  OD="${OUTROOT}/patients/$SUB/norm_final"
  TMPD="${TMPROOT}/patients/$SUB"
  mkdir -p "$OD/qc" "$TMPD"

  # Reorient inputs → TMP
  T1_REO="$TMPD/${SUB}_T1w_reo.nii.gz"
  LES_REO="$TMPD/${SUB}_lesion_reo.nii.gz"
  fslreorient2std "$T1" "$T1_REO"
  fslreorient2std "$LES" "$LES_REO"

  # Lesion size & crude bilateral check
  nvox=$(fslstats "$LES_REO" -V | awk '{print $1+0}')
  fslmaths "$LES_REO" -roi 0 -1 0 -1 0 -1 0 -40 "$TMPD/right_lesion.nii.gz" >/dev/null 2>&1 || true
  right_lesion=$(fslstats "$TMPD/right_lesion.nii.gz" -V | awk '{print $1+0}')
  fslmaths "$LES_REO" -roi 40 -1 0 -1 0 -1 0 -1 "$TMPD/left_lesion.nii.gz" >/dev/null 2>&1 || true
  left_lesion=$(fslstats "$TMPD/left_lesion.nii.gz" -V | awk '{print $1+0}')

  if [[ "$nvox" -gt $LESION_VOX_THRESHOLD && "$right_lesion" -gt 100 && "$left_lesion" -gt 100 ]]; then
    echo "CLASSIC" > "$TMPD/method.txt"
    logmsg "[PAT][$SUB] Classic method (large/bilateral: $nvox vox)" | tee -a "${OUTROOT}/logs/normalize_patients.log"
    COSTMASK="$TMPD/${SUB}_costmask.nii.gz"
    fslmaths "$LES_REO" -bin -mul -1 -add 1 "$COSTMASK"
    echo "$COSTMASK" > "$TMPD/costmask.txt"
    echo "$T1_REO" > "$TMPD/t1reg.txt"
  else
    echo "ENANTIO" > "$TMPD/method.txt"
    logmsg "[PAT][$SUB] Enantiomorphic method (small/unilateral: $nvox vox)" | tee -a "${OUTROOT}/logs/normalize_patients.log"
    BASE=$(basename "$T1_REO" .nii.gz)
    # --- Enantiomorphic fill (BCB style) ---
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

done

# --------------- preprocessing: controls ---------------
logmsg "Step 1B: Preprocess CONTROLS (serial, max threads=$NTHREADS)"
for SUB in "${CONTROL_SUBJECTS[@]:-}"; do
  T1="${CONTROL_INROOT}/${SUB}_T1w.nii.gz"
  OD="${OUTROOT}/controls/$SUB/norm_final"
  TMPD="${TMPROOT}/controls/$SUB"
  mkdir -p "$OD/qc" "$TMPD"

  # Reorient T1
  T1_REO="$TMPD/${SUB}_T1w_reo.nii.gz"
  fslreorient2std "$T1" "$T1_REO"

  # Controls: no lesion; mark method
  echo "CTRL" > "$TMPD/method.txt"
  echo "" > "$TMPD/costmask.txt"
  echo "$T1_REO" > "$TMPD/t1reg.txt"

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

done

# --------------- SyN registration (parallel) ---------------
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$THREADS_PER_SYN
export OMP_NUM_THREADS=$THREADS_PER_SYN

logmsg "Step 2: Parallel SyN registration (jobs=${MAX_PARALLEL_SYN}, threads/job=${THREADS_PER_SYN})"

syn_registration() {
  local kind="$1"   # patients|controls
  local SUB="$2"
  local TMPD
  local OD
  if [[ "$kind" == "patients" ]]; then
    TMPD="${TMPROOT}/patients/$SUB"
    OD="${OUTROOT}/patients/$SUB/norm_final"
  else
    TMPD="${TMPROOT}/controls/$SUB"
    OD="${OUTROOT}/controls/$SUB/norm_final"
  fi

  local METHOD T1_BRAIN COSTMASK ANTSOUT
  METHOD=$(cat "$TMPD/method.txt")
  T1_BRAIN=$(cat "$TMPD/t1_brain_final.txt")
  COSTMASK=$(cat "$TMPD/costmask.txt")
  ANTSOUT="$TMPD/${SUB}_xfm_"

  mkdir -p "$OD"

  local ANTS_CMD
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

  if [[ "$METHOD" == "CLASSIC" && -n "$COSTMASK" ]]; then
    ANTS_CMD+=" --masks [$COSTMASK,NULL]"
  fi

  echo "[INFO][$kind][$SUB] Starting ANTs registration..." | tee -a "$OD/ants_registration.log"
  eval $ANTS_CMD | tee -a "$OD/ants_registration.log"
  echo "[INFO][$kind][$SUB] Finished ANTs registration." | tee -a "$OD/ants_registration.log"
}
export -f syn_registration
export MNI_TEMPLATE MNI_TEMPLATE_WSKULL TMPROOT OUTROOT THREADS_PER_SYN

# Build job lists for parallel
P_JOBS=()
for s in "${PATIENT_SUBJECTS[@]:-}"; do P_JOBS+=("patients $s"); done
C_JOBS=()
for s in "${CONTROL_SUBJECTS[@]:-}"; do C_JOBS+=("controls $s"); done

parallel -j "$MAX_PARALLEL_SYN" --colsep ' ' syn_registration ::: "${P_JOBS[@]:-}" "${C_JOBS[@]:-}"

# --------------- Post-SyN (apply transforms & QC) ---------------
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NTHREADS
export OMP_NUM_THREADS=$NTHREADS

apply_and_qc_pat() {
  local SUB="$1"
  local TMPD="${TMPROOT}/patients/$SUB"
  local OD="${OUTROOT}/patients/$SUB/norm_final"
  mkdir -p "$OD" "$OD/qc"

  local T1_BRAIN LES_REO ANTSOUT
  T1_BRAIN=$(cat "$TMPD/t1_brain_final.txt")
  LES_REO="$TMPD/${SUB}_lesion_reo.nii.gz"
  ANTSOUT="$TMPD/${SUB}_xfm_"

  # T1 → MNI
  antsApplyTransforms -d 3 \
    -i "$T1_BRAIN" \
    -r "$MNI_TEMPLATE" \
    -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
    -o "$OD/${SUB}_T1_MNI.nii.gz" \
    -n Linear
  fslcpgeom "$MNI_TEMPLATE" "$OD/${SUB}_T1_MNI.nii.gz"

  # Lesion → MNI (NN)
  antsApplyTransforms -d 3 \
    -i "$LES_REO" \
    -r "$MNI_TEMPLATE" \
    -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
    -o "$OD/${SUB}_lesion_MNI.nii.gz" \
    -n NearestNeighbor
  fslcpgeom "$MNI_TEMPLATE" "$OD/${SUB}_lesion_MNI.nii.gz"

  # QC (lesion over MNI)
  overlay 1 1 "$MNI_TEMPLATE" -a "$OD/${SUB}_lesion_MNI.nii.gz" 0.5 1.0 "$OD/qc/${SUB}_lesion_on_mni.nii.gz"
  slices "$OD/qc/${SUB}_lesion_on_mni.nii.gz" -o "$OD/qc/${SUB}_lesion_on_mni.png"
  rm -f "$OD/qc/${SUB}_lesion_on_mni.nii.gz"
  echo "[OK][PAT][$SUB] T1 & lesion normalized + QC" | tee -a "${OUTROOT}/logs/post_patients.log"
}
export -f apply_and_qc_pat

apply_and_qc_ctrl() {
  local SUB="$1"
  local TMPD="${TMPROOT}/controls/$SUB"
  local OD="${OUTROOT}/controls/$SUB/norm_final"
  mkdir -p "$OD" "$OD/qc"

  local T1_BRAIN ANTSOUT
  T1_BRAIN=$(cat "$TMPD/t1_brain_final.txt")
  ANTSOUT="$TMPD/${SUB}_xfm_"

  # T1 → MNI
  antsApplyTransforms -d 3 \
    -i "$T1_BRAIN" \
    -r "$MNI_TEMPLATE" \
    -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
    -o "$OD/${SUB}_T1_MNI.nii.gz" \
    -n Linear
  fslcpgeom "$MNI_TEMPLATE" "$OD/${SUB}_T1_MNI.nii.gz"

  # QC (T1 over MNI): generate a quick montage of the normalized T1 alone
  slices "$OD/${SUB}_T1_MNI.nii.gz" -o "$OD/qc/${SUB}_t1_on_mni.png"
  echo "[OK][CTRL][$SUB] T1 normalized + QC" | tee -a "${OUTROOT}/logs/post_controls.log"
}
export -f apply_and_qc_ctrl

logmsg "Step 3A: Post-SyN (patients)"
if [[ ${#PATIENT_SUBJECTS[@]:-0} -gt 0 ]]; then
  parallel -j "$MAX_PARALLEL_SYN" apply_and_qc_pat ::: "${PATIENT_SUBJECTS[@]}"
fi

logmsg "Step 3B: Post-SyN (controls)"
if [[ ${#CONTROL_SUBJECTS[@]:-0} -gt 0 ]]; then
  parallel -j "$MAX_PARALLEL_SYN" apply_and_qc_ctrl ::: "${CONTROL_SUBJECTS[@]}"
fi

logmsg "All done! Intermediates preserved under $TMPROOT (patients/controls)."
