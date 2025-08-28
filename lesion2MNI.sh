#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────
# HARD-CODED CONFIG — EDIT THESE IF PATHS CHANGE
# ─────────────────────────────────────────────────────────────
INROOT="/media/sskgroup/fast/Lesion/Secondary"                        # folder with sub-XX_T1w.nii.gz + sub-XX_label-lesion_roi.nii.gz
OUTROOT="/media/sskgroup/fast/Lesion/Secondary_MNI"                         # where outputs go
REF="/home/sskgroup/Documents/BCBToolKit/Tools/extraFiles/restingState/MNI152_T1_2mm_brain.nii.gz"  # MNI2009cAsym 2mm brain
REFMASK="/home/sskgroup/Documents/BCBToolKit/Tools/extraFiles/MNI152_for_antsBrainExtractionMask_2mm.nii.gz"  # resampled to 2mm
THREADS=16                                  # how many CPU threads for ANTs
ALL=1                                      # 1=batch all subjects, 0=single subject
SUB="sub-"                               # subject ID if ALL=0
# ─────────────────────────────────────────────────────────────

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS

die(){ echo "[ERR] $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }

same_grid(){
  local a="$1" b="$2"
  local ad bd ap bp
  ad=$(fslhd "$a" | awk '/^dim[123]/ {printf "%s ",$2}'); ad="${ad% }"
  bd=$(fslhd "$b" | awk '/^dim[123]/ {printf "%s ",$2}'); bd="${bd% }"
  ap=$(fslhd "$a" | awk '/^pixdim[123]/ {printf "%s ",$2}'); ap="${ap% }"
  bp=$(fslhd "$b" | awk '/^pixdim[123]/ {printf "%s ",$2}'); bp="${bp% }"
  [[ "$ad" == "$bd" && "$ap" == "$bp" ]]
}

resample_to_ref(){
  antsApplyTransforms -d 3 -i "$1" -r "$2" -t identity -n "$4" -o "$3"
  fslcpgeom "$2" "$3" || true
}

normalize_subject(){
  local SUB="$1"
  local T1="${INROOT}/${SUB}_T1w.nii.gz"
  local LES="${INROOT}/${SUB}_label-lesion_roi.nii.gz"
  [[ -f "$T1" ]] || die "Missing $T1"
  [[ -f "$LES" ]] || die "Missing $LES"

  local OUTDIR="${OUTROOT}/${SUB}/norm_v3"
  local LOGDIR="${OUTDIR}/logs"
  mkdir -p "$OUTDIR" "$LOGDIR"

  # N4 bias correction
  if [[ ! -f "${OUTDIR}/${SUB}_T1w_N4.nii.gz" ]]; then
    echo "[*] $SUB: N4 bias correction"
    N4BiasFieldCorrection -d 3 -i "$T1" -o "${OUTDIR}/${SUB}_T1w_N4.nii.gz"
  fi
  local T1N4="${OUTDIR}/${SUB}_T1w_N4.nii.gz"

  # lesion + moving mask
  fslmaths "$LES" -bin "${OUTDIR}/${SUB}_lesion_bin_native.nii.gz"
  if [[ -f "${INROOT}/${SUB}_brainmask.nii.gz" ]]; then
    fslmaths "${INROOT}/${SUB}_brainmask.nii.gz" -sub "${OUTDIR}/${SUB}_lesion_bin_native.nii.gz" -thr 0 -bin "${OUTDIR}/${SUB}_moving_mask.nii.gz"
  else
    fslmaths "${OUTDIR}/${SUB}_lesion_bin_native.nii.gz" -mul -1 -add 1 -thr 0 -bin "${OUTDIR}/${SUB}_moving_mask.nii.gz"
  fi

  # fixed mask argument
  local FIXMASK_ARG=()
  if [[ -n "$REFMASK" && -f "$REFMASK" ]]; then
    if same_grid "$REF" "$REFMASK"; then
      FIXMASK_ARG=(-x "${OUTDIR}/${SUB}_moving_mask.nii.gz,${REFMASK}")
    else
      echo "[*] $SUB: resampling fixed mask"
      local RM="${OUTDIR}/refmask_2mm.nii.gz"
      resample_to_ref "$REFMASK" "$REF" "$RM" NearestNeighbor
      FIXMASK_ARG=(-x "${OUTDIR}/${SUB}_moving_mask.nii.gz,${RM}")
    fi
  else
    FIXMASK_ARG=(-x "${OUTDIR}/${SUB}_moving_mask.nii.gz")
  fi

  # ANTs registration
  echo "[*] $SUB: ANTs SyN → MNI"
  set +e
  antsRegistrationSyN.sh -d 3 -f "$REF" -m "$T1N4" -o "${OUTDIR}/${SUB}_xfm_" -n "$THREADS" "${FIXMASK_ARG[@]}" \
    > "${LOGDIR}/${SUB}_ants_reg.log" 2>&1
  RC=$?
  set -e
  if [[ $RC -ne 0 || ! -f "${OUTDIR}/${SUB}_xfm_1Warp.nii.gz" ]]; then
    echo "[WARN] $SUB: SyN failed, fallback to affine"
    antsRegistrationSyN.sh -d 3 -f "$REF" -m "$T1N4" -o "${OUTDIR}/${SUB}_xfm_" -n "$THREADS" -t a "${FIXMASK_ARG[@]}" \
      >> "${LOGDIR}/${SUB}_ants_reg.log" 2>&1
  fi

  # apply to T1
  antsApplyTransforms -d 3 -i "$T1N4" -r "$REF" \
    -t "${OUTDIR}/${SUB}_xfm_1Warp.nii.gz" -t "${OUTDIR}/${SUB}_xfm_0GenericAffine.mat" \
    -o "${OUTDIR}/${SUB}_T1_MNI2mm.nii.gz"
  fslcpgeom "$REF" "${OUTDIR}/${SUB}_T1_MNI2mm.nii.gz" || true

  # apply to lesion
  antsApplyTransforms -d 3 -i "${OUTDIR}/${SUB}_lesion_bin_native.nii.gz" -r "$REF" \
    -t "${OUTDIR}/${SUB}_xfm_1Warp.nii.gz" -t "${OUTDIR}/${SUB}_xfm_0GenericAffine.mat" \
    -n NearestNeighbor -o "${OUTDIR}/${SUB}_lesion_MNI2mm.nii.gz"
  fslcpgeom "$REF" "${OUTDIR}/${SUB}_lesion_MNI2mm.nii.gz" || true
  fslmaths "${OUTDIR}/${SUB}_lesion_MNI2mm.nii.gz" -bin "${OUTDIR}/${SUB}_lesion_MNI2mm.nii.gz"

  # QC
  local V
  V=$(fslstats "${OUTDIR}/${SUB}_lesion_MNI2mm.nii.gz" -V)
  echo "[QC] $SUB lesion voxels=$(echo $V | awk '{print $1}')  volume_mm3=$(echo $V | awk '{print $2}')"
  echo "[OK] $SUB done"
}

# ---- main run
need antsRegistrationSyN.sh; need antsApplyTransforms; need N4BiasFieldCorrection
need fslmaths; need fslstats; need fslcpgeom; need fslhd

if [[ $ALL -eq 1 ]]; then
  mapfile -t SUBS < <(find "$INROOT" -maxdepth 1 -type f -name "sub-*_T1w.nii.gz" -printf "%f\n" | sed -E 's/_T1w\.nii\.gz$//' | sort)
  [[ ${#SUBS[@]} -gt 0 ]] || die "No subjects found in $INROOT"
  for S in "${SUBS[@]}"; do
    echo "========== $S =========="
    normalize_subject "$S"
  done
  echo "[OK] batch complete"
else
  normalize_subject "$SUB"
fi