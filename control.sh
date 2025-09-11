#!/usr/bin/env bash
set -euo pipefail

# ╭──────────────────────────────────────────────────────────────╮
# │  fetch_controls_to_BCB_bridge.sh (ONE BATCH MODE, HARD-CODED)│
# │  Convert CONTROL T1s from fMRIPrep batch → BCB MNI space     │
# │  Strategy: Template bridge (TF MNI2009cAsym → BCB MNI152)    │
# │  Process all subjects inside a single batch dir              │
# ╰──────────────────────────────────────────────────────────────╯

# ===================== HARD-CODED PATHS =====================
BATCH_DIR="/media/sskgroup/Thesis/functional/fmriprep_24/batch_6"
OUTROOT="/media/sskgroup/Thesis/functional/Controls_BCB"
TMPROOT="/media/sskgroup/Thesis/functional/Controls_BCB_TMP"
BCBTOOLS="/home/sskgroup/Documents/BCBToolKit/Tools"

# Precomputed bridge transforms
TF2BCB_AFF="/media/sskgroup/Thesis/functional/Controls_BCB/bridge/tflow_to_bcb_0GenericAffine.mat"
TF2BCB_WARP="/media/sskgroup/Thesis/functional/Controls_BCB/bridge/tflow_to_bcb_1Warp.nii.gz"

# Templates
TFLOW_T1_PATTERN="*_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz"
MNI_BCB_T1="${BCBTOOLS}/extraFiles/MNI152.nii.gz"

# Threading (not used in loop mode; kept for reference)
NTHREADS=16
THREADS_PER_SYN=2
# ============================================================

export FSLOUTPUTTYPE=NIFTI_GZ

need(){ command -v "$1" >/dev/null 2>&1 || { echo "[FATAL] Missing: $1"; exit 1; }; }
need fslmaths; need fslstats; need fslreorient2std; need fslcpgeom
need antsApplyTransforms; need slices

log(){ echo "[$(date +"%F %T")] $*"; }

mkdir -p "$OUTROOT/controls" "$OUTROOT/logs" "$TMPROOT"

# ---------------- Check bridge files ----------------
if [[ ! -f "$TF2BCB_AFF" || ! -f "$TF2BCB_WARP" ]]; then
  log "[FATAL] Missing bridge transforms in $(dirname $TF2BCB_AFF)"; exit 1
fi

# ---------------- Worker ----------------
process_control_bridge(){ # $1: subject ID (e.g., sub-01)
  local SUB="$1"
  local AD="$BATCH_DIR/$SUB/anat"
  local OUTD="$OUTROOT/controls/$SUB/norm_final"; mkdir -p "$OUTD/qc"

  local TF_T1
  TF_T1=$(find "$AD" -maxdepth 1 -type f -name "$TFLOW_T1_PATTERN" | head -n 1 || true)
  if [[ -z "${TF_T1:-}" ]]; then
    log "[WARN][$SUB] No ${TFLOW_T1_PATTERN} in $AD — skipping" | tee -a "$OUTROOT/logs/controls_bridge.log"
    return
  fi

  antsApplyTransforms -d 3 \
    -i "$TF_T1" \
    -r "$MNI_BCB_T1" \
    -t "$TF2BCB_WARP" -t "$TF2BCB_AFF" \
    -o "$OUTD/${SUB}_T1_BCBMNI.nii.gz" \
    -n Linear
  fslcpgeom "$MNI_BCB_T1" "$OUTD/${SUB}_T1_BCBMNI.nii.gz"
  slices "$OUTD/${SUB}_T1_BCBMNI.nii.gz" -o "$OUTD/qc/${SUB}_t1_on_bcb.png"
  log "[OK][BRIDGE][$SUB] → $OUTD/${SUB}_T1_BCBMNI.nii.gz" | tee -a "$OUTROOT/logs/controls_bridge.log"
}

# ---------------- Collect subjects from batch ----------------
mapfile -t SUBS < <(find "$BATCH_DIR" -maxdepth 1 -type d -name 'sub-*' -printf '%f\n' | sort)

log "Planned ${#SUBS[@]} control conversions in $BATCH_DIR"

for SUB in "${SUBS[@]}"; do
  process_control_bridge "$SUB"
done

log "All controls in $BATCH_DIR processed. Outputs under $OUTROOT/controls/<SUB>/norm_final/"
