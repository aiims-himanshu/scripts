#!/usr/bin/env bash
set -euo pipefail

# ╭─────────────────────────────────────────────────────────────╮
# │    MASTER LESION NORMALIZATION PIPELINE FOR BCBTOOLKIT       │
# │  Extension-Agnostic with Automated Resume Capabilities      │
# ╰─────────────────────────────────────────────────────────────╯

# ======= USER CONFIGURATION =======
INROOT="/home/sskgroup/Documents/Lesion/Secondary"  
OUTROOT="/home/sskgroup/Documents/Lesion/Secondary_MNI"
BCBTOOLS="/home/sskgroup/Documents/BCBToolKit/Tools"

# Core BCBtoolkit Reference Files
MNI_TEMPLATE="${BCBTOOLS}/extraFiles/MNI152.nii.gz"
REFMASK="${BCBTOOLS}/extraFiles/MNI152_for_antsBrainExtractionMask.nii.gz"
BRAINTPL="${BCBTOOLS}/extraFiles/MNI152_for_antsBrainExtractionBrain.nii.gz"
ENFILL_SCRIPT="/home/sskgroup/Git/scripts/enantiomorphic-fill-v1.sh"

N_THREADS=16
LESION_VOX_THRESHOLD=1000  
export FSLOUTPUTTYPE="NIFTI_GZ" # Keeps outputs strictly compressed to save space
# ==================================

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${N_THREADS}
export OMP_NUM_THREADS=${N_THREADS}

mkdir -p "$OUTROOT"

cd "$INROOT"

# FIX 1: Dynamic glob matching both .nii and .nii.gz extensions for T1 scans
for T1 in sub-*_*T1w.nii sub-*_*T1w.nii.gz; do
  # Break safely if no matching files exist in the directory
  [[ -e "$T1" ]] || continue

  # FIX 2: Dynamic regex pattern extraction to isolate the clean subject ID string
  SUB=$(basename "$T1" | sed -E 's/_T1w\.nii(\.gz)?$//')
  
  # FIX 3: Locate the lesion mask whether it ends in .nii or .nii.gz
  LES=""
  if [[ -f "${INROOT}/${SUB}_label-lesion_roi.nii.gz" ]]; then
    LES="${INROOT}/${SUB}_label-lesion_roi.nii.gz"
  elif [[ -f "${INROOT}/${SUB}_label-lesion_roi.nii" ]]; then
    LES="${INROOT}/${SUB}_label-lesion_roi.nii"
  else
    echo "[WARN] Skipping $SUB: No lesion mask found (.nii or .nii.gz)."
    continue
  fi

  OD="${OUTROOT}/${SUB}/norm_final"
  TMPD="${OUTROOT}/${SUB}/tmp_norm"

  # --- CRITICAL ADDITION: SMART TRACKING RESUME CHECK ---
  FINAL_LESION_OUT="${OD}/${SUB}_lesion_MNI1mm.nii.gz"
  FINAL_T1_OUT="${OD}/${SUB}_T1_MNI1mm.nii.gz"
  FINAL_QC_OUT="${OD}/qc/${SUB}_lesion_on_mni.png"

  if [[ -f "$FINAL_LESION_OUT" && -f "$FINAL_T1_OUT" && -f "$FINAL_QC_OUT" ]]; then
    echo "[INFO] Skip [$SUB]: Process outputs are already fully complete."
    continue
  fi

  # If any output was missing, ensure directories exist and begin calculation
  mkdir -p "$OD" "$OD/logs" "$OD/qc" "$TMPD"

  echo "[*] [$SUB] Reorienting headers to standard MNI coordinates..."
  T1_REO="$TMPD/${SUB}_T1w_reo.nii.gz"
  LES_REO="$TMPD/${SUB}_lesion_reo.nii.gz"
  
  # FSL tools handle mixed input extensions natively and convert them to .nii.gz cleanly here
  fslreorient2std "$T1" "$T1_REO"
  fslreorient2std "$LES" "$LES_REO"

  # Measure lesion properties
  nvox=$(fslstats "$LES_REO" -V | awk '{print $1+0}')
  
  FINAL_ANATOMY_BRAIN="${OD}/${SUB}_T1_brain.nii.gz"
  USE_COSTMASK=0

  # --- DECISION ROUTINE ---
  if [[ "$nvox" -le $LESION_VOX_THRESHOLD ]]; then
    echo "[*] [$SUB] Target fits unilateral profile ($nvox vox). Launching Enantiomorphic fill..."
    
    if [[ ! -f "$ENFILL_SCRIPT" ]]; then
      echo "[ERR] Enantiomorphic fill script not found at: $ENFILL_SCRIPT"
      exit 1
    fi
    
    # Execute structural patching
    bash "$ENFILL_SCRIPT" "$SUB" "$INROOT" "$OUTROOT" 1 $N_THREADS
    FILLED="${OUTROOT}/${SUB}/enfill_v1/${SUB}_T1w_enantiomorphic.nii.gz"
    
    if [[ -f "$FILLED" ]]; then
      echo "[*] [$SUB] Patch generated. Running brain extraction on filled structural image..."
      antsBrainExtraction.sh -d 3 -a "$FILLED" -e "$BRAINTPL" -m "$REFMASK" -o "${OD}/${SUB}_filled_be_"
      mv "${OD}/${SUB}_filled_be_BrainExtractionBrain.nii.gz" "$FINAL_ANATOMY_BRAIN"
      rm -f ${OD}/${SUB}_filled_be_*
    else
      echo "[WARN] [$SUB] Patch tracking failed; defaulting to raw structural extraction."
      antsBrainExtraction.sh -d 3 -a "$T1_REO" -e "$BRAINTPL" -m "$REFMASK" -o "${OD}/${SUB}_be_"
      mv "${OD}/${SUB}_be_BrainExtractionBrain.nii.gz" "$FINAL_ANATOMY_BRAIN"
      rm -f ${OD}/${SUB}_be_*
    fi
  else
    echo "[*] [$SUB] Target exceeds localized threshold ($nvox vox). Initializing Cost-Masking configuration..."
    COSTMASK="$TMPD/${SUB}_costmask.nii.gz"
    fslmaths "$LES_REO" -bin -mul -1 -add 1 "$COSTMASK"
    fslmaths "$COSTMASK" -dilM -dilM "$TMPD/${SUB}_costmask_dil.nii.gz"
    COSTMASK_FINAL="$TMPD/${SUB}_costmask_dil.nii.gz"
    USE_COSTMASK=1

    # Extract regular native brain
    antsBrainExtraction.sh -d 3 -a "$T1_REO" -e "$BRAINTPL" -m "$REFMASK" -o "${OD}/${SUB}_be_"
    mv "${OD}/${SUB}_be_BrainExtractionBrain.nii.gz" "$FINAL_ANATOMY_BRAIN"
    rm -f ${OD}/${SUB}_be_*
  fi

  # --- CREATE DILATED COST MASK FOR FINAL SyN ---
  DILM="${OD}/${SUB}_lesion_dil.nii.gz"
  fslmaths "$LES_REO" -bin -dilM -dilM "$DILM"

  # --- REGISTER CLEAN ANATOMICAL BRAIN TO MNI TEMPLATE ---
  ANTSOUT="${OD}/${SUB}_xfm_"
  echo "[*] [$SUB] Launching symmetric diffeomorphic mapping (SyN)..."
  
  if [[ "$USE_COSTMASK" -eq 1 ]]; then
    antsRegistrationSyN.sh -d 3 -f "$BRAINTPL" -m "$FINAL_ANATOMY_BRAIN" -o "$ANTSOUT" -n "$N_THREADS" -x "$COSTMASK_FINAL" -t s > "$OD/logs/${SUB}_ants_reg.log" 2>&1
  else
    antsRegistrationSyN.sh -d 3 -f "$BRAINTPL" -m "$FINAL_ANATOMY_BRAIN" -o "$ANTSOUT" -n "$N_THREADS" -x "$DILM" -t s > "$OD/logs/${SUB}_ants_reg.log" 2>&1
  fi

  # --- WARP LESION MASK TO MNI (CRITICAL: NEAREST NEIGHBOR & MATRIX COPY) ---
  echo "[*] [$SUB] Resampling lesion mask to standard space grid..."
  antsApplyTransforms -d 3 \
    -i "$LES_REO" \
    -r "$MNI_TEMPLATE" \
    -t "${ANTSOUT}1Warp.nii.gz" \
    -t "${ANTSOUT}0GenericAffine.mat" \
    -n NearestNeighbor \
    -o "$FINAL_LESION_OUT"

  # Force exact geometric boundary template dimensions matching BCBtoolkit
  fslcpgeom "$MNI_TEMPLATE" "$FINAL_LESION_OUT"

  # --- WARP T1 ANATOMY TO MNI (LINEAR) ---
  antsApplyTransforms -d 3 -i "$FINAL_ANATOMY_BRAIN" -r "$MNI_TEMPLATE" -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" -n Linear -o "$FINAL_T1_OUT"

  # --- VISUAL VERIFICATION AND CLEANUP ---
  echo "[*] [$SUB] Exporting Quality Control overlays..."
  overlay 1 1 "$MNI_TEMPLATE" -a "$FINAL_LESION_OUT" 0.5 1.0 "${TMPD}/qc_render.nii.gz"
  slices "${TMPD}/qc_render.nii.gz" -o "$FINAL_QC_OUT"
  
  rm -rf "$TMPD"
  echo "[OK] Disconnectome file generation complete: $FINAL_LESION_OUT"
done
