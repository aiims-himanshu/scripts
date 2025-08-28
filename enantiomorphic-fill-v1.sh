#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────
# enantiomorphic-fill-v1.sh
# Create an enantiomorphically-filled T1 in native space.
#
# USAGE:
#   bash enantiomorphic-fill-v1.sh <SUBJECT> <INROOT> <OUTROOT> <DO_RIGID> <THREADS>
#     SUBJECT   : subject ID (e.g., sub-93)
#     INROOT    : input dir with ${SUBJECT}_T1w.nii.gz and ${SUBJECT}_label-lesion_roi.nii.gz
#     OUTROOT   : output root (writes ${OUTROOT}/${SUBJECT}/enfill_v1/)
#     DO_RIGID  : 1 = use rigid alignment (recommended), 0 = skip (fast)
#     THREADS   : number of ANTs threads to use
#
# DEPENDS: FSL (fslmaths, fslcpgeom, fslreorient2std, fslswapdim, overlay, slices)
#          ANTs (antsRegistrationSyN.sh, antsApplyTransforms) [if DO_RIGID=1]
# ─────────────────────────────────────────────────────────────

SUB="${1:?sub-XX}"
INROOT="${2:?/path/to/inroot}"
OUTROOT="${3:?/path/to/outroot}"
DO_RIGID="${4:-1}"
THREADS="${5:-8}"

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS

need(){ command -v "$1" >/dev/null 2>&1 || { echo "[ERR] Missing: $1" >&2; exit 1; }; }
need fslmaths; need fslcpgeom; need fslreorient2std; need fslswapdim; need overlay; need slices
if [[ "$DO_RIGID" -eq 1 ]]; then
  need antsRegistrationSyN.sh; need antsApplyTransforms
fi

T1="${INROOT}/${SUB}_T1w.nii.gz"
LES="${INROOT}/${SUB}_label-lesion_roi.nii.gz"
[[ -f "$T1"  ]] || { echo "[ERR] Missing $T1"; exit 1; }
[[ -f "$LES" ]] || { echo "[ERR] Missing $LES"; exit 1; }

OD="${OUTROOT}/${SUB}/enfill_v1"
LD="${OD}/logs"; QCD="${OD}/qc"
mkdir -p "$OD" "$LD" "$QCD"

echo "[*] $SUB: reorienting headers to std"
fslreorient2std "$T1"  "${OD}/${SUB}_T1w_std.nii.gz"
fslreorient2std "$LES" "${OD}/${SUB}_lesion_std.nii.gz"

# Bin and (lightly) dilate lesion to ensure full coverage during replacement
echo "[*] $SUB: lesion bin + 2× dilation"
fslmaths "${OD}/${SUB}_lesion_std.nii.gz" -bin "${OD}/${SUB}_lesion_bin.nii.gz"
fslmaths "${OD}/${SUB}_lesion_bin.nii.gz" -dilM -dilM "${OD}/${SUB}_lesion_bin_dil.nii.gz"

# Create inverse mask (1 - lesion)
fslmaths "${OD}/${SUB}_lesion_bin_dil.nii.gz" -mul -1 -add 1 -bin "${OD}/${SUB}_lesion_inv.nii.gz"

# Mirror the T1 in native grid (left-right flip around image center)
echo "[*] $SUB: generating mirrored T1"
fslswapdim "${OD}/${SUB}_T1w_std.nii.gz" -x y z "${OD}/${SUB}_T1w_mirror_raw.nii.gz"
fslcpgeom "${OD}/${SUB}_T1w_std.nii.gz" "${OD}/${SUB}_T1w_mirror_raw.nii.gz"

MIR_ALIGNED="${OD}/${SUB}_T1w_mirror_on_orig.nii.gz"

if [[ "$DO_RIGID" -eq 1 ]]; then
  echo "[*] $SUB: rigid alignment of mirrored→original with lesion-excluded masks"

  # Make a simple heuristic brainmask (or supply your own)
  if [[ ! -f "${OD}/${SUB}_brainmask_std.nii.gz" ]]; then
    fslmaths "${OD}/${SUB}_T1w_std.nii.gz" -bin -dilM -s 2 -thr 0.1 -bin "${OD}/${SUB}_brainmask_std.nii.gz"
  fi

  # Moving (mirrored) masks: flip brainmask & lesion to mirrored space, copy geometry
  fslswapdim "${OD}/${SUB}_brainmask_std.nii.gz" -x y z "${OD}/${SUB}_brainmask_mirror.nii.gz"
  fslswapdim "${OD}/${SUB}_lesion_bin_dil.nii.gz" -x y z "${OD}/${SUB}_lesion_mirror_dil.nii.gz"
  fslcpgeom "${OD}/${SUB}_T1w_mirror_raw.nii.gz" "${OD}/${SUB}_brainmask_mirror.nii.gz"
  fslcpgeom "${OD}/${SUB}_T1w_mirror_raw.nii.gz" "${OD}/${SUB}_lesion_mirror_dil.nii.gz"

  # Fixed and moving cost-function masks (exclude lesions)
  fslmaths "${OD}/${SUB}_brainmask_std.nii.gz"   -sub "${OD}/${SUB}_lesion_bin_dil.nii.gz"   -thr 0 -bin "${OD}/${SUB}_fixed_mask.nii.gz"
  fslmaths "${OD}/${SUB}_brainmask_mirror.nii.gz" -sub "${OD}/${SUB}_lesion_mirror_dil.nii.gz" -thr 0 -bin "${OD}/${SUB}_moving_mask.nii.gz"

  # Rigid registration (mirrored → original)
  set +e
  antsRegistrationSyN.sh -d 3 \
    -f "${OD}/${SUB}_T1w_std.nii.gz" \
    -m "${OD}/${SUB}_T1w_mirror_raw.nii.gz" \
    -o "${OD}/${SUB}_rigid_" \
    -t r -n "$THREADS" \
    -x "${OD}/${SUB}_moving_mask.nii.gz,${OD}/${SUB}_fixed_mask.nii.gz" \
    >> "${LD}/${SUB}_rigid.log" 2>&1
  RC=$?
  set -e

  if [[ $RC -ne 0 || ! -f "${OD}/${SUB}_rigid_0GenericAffine.mat" ]]; then
    echo "[WARN] $SUB: rigid alignment failed; falling back to no-alignment"
    cp "${OD}/${SUB}_T1w_mirror_raw.nii.gz" "$MIR_ALIGNED"
  else
    antsApplyTransforms -d 3 -i "${OD}/${SUB}_T1w_mirror_raw.nii.gz" -r "${OD}/${SUB}_T1w_std.nii.gz" \
      -t "${OD}/${SUB}_rigid_0GenericAffine.mat" -o "$MIR_ALIGNED"
    fslcpgeom "${OD}/${SUB}_T1w_std.nii.gz" "$MIR_ALIGNED" || true
  fi
else
  echo "[i] $SUB: DO_RIGID=0 → skipping rigid alignment"
  cp "${OD}/${SUB}_T1w_mirror_raw.nii.gz" "$MIR_ALIGNED"
fi

# Compose the filled T1: original outside lesion, mirrored inside lesion
echo "[*] $SUB: composing enantiomorphic fill"
fslmaths "${OD}/${SUB}_T1w_std.nii.gz" -mas "${OD}/${SUB}_lesion_inv.nii.gz" "${OD}/${SUB}_T1_bg.nii.gz"
fslmaths "$MIR_ALIGNED" -mas "${OD}/${SUB}_lesion_bin_dil.nii.gz" "${OD}/${SUB}_T1_fillpart.nii.gz"
fslmaths "${OD}/${SUB}_T1_bg.nii.gz" -add "${OD}/${SUB}_T1_fillpart.nii.gz" "${OD}/${SUB}_T1w_enantiomorphic.nii.gz"

# QC: overlays
overlay 1 1 "${OD}/${SUB}_T1w_std.nii.gz" -a "${OD}/${SUB}_lesion_bin_dil.nii.gz" 0.5 1.0 "${QCD}/${SUB}_lesion_on_T1.nii.gz"
slices  "${QCD}/${SUB}_lesion_on_T1.nii.gz" -o "${QCD}/${SUB}_lesion_on_T1.png"
rm -f   "${QCD}/${SUB}_lesion_on_T1.nii.gz"

overlay 1 1 "${OD}/${SUB}_T1w_enantiomorphic.nii.gz" -a "${OD}/${SUB}_lesion_bin_dil.nii.gz" 0.5 1.0 "${QCD}/${SUB}_lesion_on_T1filled.nii.gz"
slices  "${QCD}/${SUB}_lesion_on_T1filled.nii.gz" -o "${QCD}/${SUB}_lesion_on_T1filled.png"
rm -f   "${QCD}/${SUB}_lesion_on_T1filled.nii.gz"

echo "[OK] $SUB → ${OD}/${SUB}_T1w_enantiomorphic.nii.gz"
