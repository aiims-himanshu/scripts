#!/usr/bin/env bash
set -euo pipefail

# ╭─────────────────────────────────────────────────────────────╮
# │ Batch normalize all lesion masks in a directory to MNI1mm   │
# │ including enantiomorphic fill when needed.                  │
# ╰─────────────────────────────────────────────────────────────╯

# ======= USER: EDIT THESE =======
INROOT="/media/sskgroup/Thesis/lesion/Secondary"  # T1, lesion in native space
OUTROOT="/media/sskgroup/Thesis/lesion/Secondary_MNI"
MNI1="/home/sskgroup/Documents/BCBToolKit/Tools/extraFiles/MNI152.nii.gz"
REFMASK="/home/sskgroup/Documents/BCBToolKit/Tools/extraFiles/MNI152_for_antsBrainExtractionMask.nii.gz"
BRAINTPL="/home/sskgroup/Documents/BCBToolKit/Tools/extraFiles/MNI152_for_antsBrainExtractionBrain.nii.gz"
ENFILL_SCRIPT="/home/sskgroup/Git/scripts/enantiomorphic-fill-v1.sh"
N_THREADS=12
LESION_VOX_THRESHOLD=5000
# ================================

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${N_THREADS}
export OMP_NUM_THREADS=${N_THREADS}

cd "$INROOT"
for T1 in sub-*_*T1w.nii.gz; do
  SUB=$(basename "$T1" | sed 's/_T1w.nii.gz//')
  LES="${INROOT}/${SUB}_label-lesion_roi.nii.gz"

  # Check if lesion mask exists
  if [[ ! -f "$LES" ]]; then
    echo "[WARN] Skipping $SUB: No lesion mask found."
    continue
  fi

  OD="${OUTROOT}/${SUB}/norm_final"
  mkdir -p "$OD" "$OD/logs" "$OD/qc"

  # --- STEP 0: Native space QC (manual step, but print notice) ---
  echo "[INFO] [$SUB] Native-space QC recommended for $T1 and $LES before proceeding."

  # --- STEP 1: Brain Extraction (ANTs) ---
  if [[ ! -f "${OD}/${SUB}_T1_brain.nii.gz" ]]; then
    echo "[*] [$SUB] Brain extraction via ANTs..."
    antsBrainExtraction.sh \
      -d 3 \
      -a "$T1" \
      -e "$BRAINTPL" \
      -m "$REFMASK" \
      -o "${OD}/${SUB}_be_"
    mv "${OD}/${SUB}_be_BrainExtractionBrain.nii.gz" "${OD}/${SUB}_T1_brain.nii.gz"
    mv "${OD}/${SUB}_be_BrainExtractionMask.nii.gz"  "${OD}/${SUB}_brainmask.nii.gz"
  fi

  # --- STEP 2: Lesion size check for enantiomorphic fill ---
  nvox=$(fslstats "$LES" -V | awk '{print $1+0}')
  USE_FILL=0
  if [[ "$nvox" -gt $LESION_VOX_THRESHOLD ]]; then
    echo "[*] [$SUB] Lesion > $LESION_VOX_THRESHOLD voxels ($nvox) -- will use enantiomorphic fill."
    USE_FILL=1
  else
    echo "[*] [$SUB] Lesion <= $LESION_VOX_THRESHOLD voxels ($nvox) -- using original T1."
  fi

  MOVIMG="${OD}/${SUB}_T1_brain.nii.gz"
  if [[ "$USE_FILL" -eq 1 ]]; then
    if [[ ! -f "$ENFILL_SCRIPT" ]]; then
      echo "[ERR] Enantiomorphic fill script not found: $ENFILL_SCRIPT"
      exit 1
    fi
    bash "$ENFILL_SCRIPT" "$SUB" "$INROOT" "$OUTROOT" 1 $N_THREADS
    FILLED="${OUTROOT}/${SUB}/enfill_v1/${SUB}_T1w_enantiomorphic.nii.gz"
    if [[ -f "$FILLED" ]]; then
      MOVIMG="$FILLED"
      echo "[*] [$SUB] Using filled T1 for registration: $FILLED"
    else
      echo "[WARN] [$SUB] Enantiomorphic fill failed; falling back to brain-extracted T1."
    fi
  fi

  # --- STEP 3: Create cost-function mask (dilated lesion) ---
  DILM="${OD}/${SUB}_lesion_dil.nii.gz"
  fslmaths "$LES" -bin -dilM -dilM "$DILM"

  # --- STEP 4: Register T1 (or filled) to MNI152 1mm with SyN, cost-func mask ---
  ANTSOUT="${OD}/${SUB}_xfm_"
  echo "[*] [$SUB] Running ANTs SyN registration (max threads)..."
  antsRegistrationSyN.sh \
    -d 3 \
    -f "$MNI1" \
    -m "$MOVIMG" \
    -o "$ANTSOUT" \
    -n "$N_THREADS" \
    -x "$DILM" \
    -t s \
    > "$OD/logs/${SUB}_ants_reg.log" 2>&1

  # --- STEP 5: Apply transform to lesion mask (NEAREST NEIGHBOR) ---
  echo "[*] [$SUB] Applying transform to lesion mask (NN)..."
  antsApplyTransforms -d 3 -i "$LES" -r "$MNI1" \
    -t "${ANTSOUT}1Warp.nii.gz" -t "${ANTSOUT}0GenericAffine.mat" \
    -n NearestNeighbor \
    -o "${OD}/${SUB}_lesion_MNI1mm.nii.gz"
  fslcpgeom "$MNI1" "${OD}/${SUB}_lesion_MNI1mm.nii.gz"

  # --- STEP 6: QC PNGs ---
  echo "[*] [$SUB] Making QC overlays..."
  overlay 1 1 "$MNI1" -a "${OD}/${SUB}_lesion_MNI1mm.nii.gz" 0.5 1.0 "${OD}/qc/${SUB}_lesion_on_mni.nii.gz"
  slices "${OD}/qc/${SUB}_lesion_on_mni.nii.gz" -o "${OD}/qc/${SUB}_lesion_on_mni.png"
  rm -f "${OD}/qc/${SUB}_lesion_on_mni.nii.gz"

  echo "[OK] [$SUB] Lesion normalized to MNI1mm: ${OD}/${SUB}_lesion_MNI1mm.nii.gz"
  echo "     QC PNG: ${OD}/qc/${SUB}_lesion_on_mni.png"
done
