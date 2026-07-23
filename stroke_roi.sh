#!/usr/bin/env bash
set -euo pipefail

# ╭──────────────────────────────────────────────────────────────╮
# │ stroke_roi.sh                                                │
# │ Build CONN-compatible lesion, disconnectome, and rings      │
# │ for Primary-style ROI trees and CONN/BIDS datasets          │
# ╰──────────────────────────────────────────────────────────────╯
#
# Default layout:
#   ROI_ROOT/
#     T1/                sub-XXX_T1.nii.gz
#     Lesion/            sub-XXX_lesion.nii.gz
#     Disconnectome/     sub-XXX_discon.nii.gz
#
#   CONN_ROOT/
#     sub-*/anat/*_T1w*.nii[.gz]
#     sub-*/func/*_run-01_bold*.nii[.gz]
#
# Example:
#   bash stroke_roi.sh \
#     --roi-root  /media/sskgroup/Thesis/Primary/roi \
#     --conn-root /media/sskgroup/Thesis/Primary/standard_lesion/data/BIDS/dataset \
#     --out-root  /media/sskgroup/Thesis/Primary/standard_lesion/conn_rois \
#     --threads 16
#
# Optional:
#   --id-map /path/to/idmap.csv   # CSV columns: conn_id,source_id
#
# Requirements: ANTs (antsRegistration, antsApplyTransforms), FSL (fslmaths, fslval, imcp optional)

# ---------- defaults ----------
ROI_ROOT="/media/sskgroup/Thesis/Primary/roi"
T1_ROOT=""
LESION_ROOT=""
DISC_ROOT=""
CONN_ROOT="/media/sskgroup/Thesis/Primary/standard_lesion/data/BIDS/dataset"
OUT_ROOT="/media/sskgroup/Thesis/Primary/conn_rois"
THREADS="${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS:-16}"
IDMAP=""

# peri-lesional ring radii (mm)
RING_MM_A=2
RING_MM_B=4
RING_MM_C=8

# ---------- helpers ----------
die(){ echo "ERROR: $*" >&2; exit 1; }
log(){ echo "[`date +%F' '%T`] $*"; }

have(){ command -v "$1" >/dev/null 2>&1; }
have_imcp(){ have imcp; }

pick_first(){ # glob and echo first existing match
  local pattern="$1"
  local f
  shopt -s nullglob
  for f in $pattern; do
    echo "$f"
    return 0
  done
  return 1
}

id_from_csv(){ # map conn_id -> source_id via CSV "conn_id,source_id"
  local csv="$1"
  local key="$2"
  awk -F, -v k="$key" 'NR>1 && $1==k {print $2; exit}' "$csv"
}

id_variants(){
  local src="$1" core
  if [[ "$src" =~ ^sub-([0-9]+)$ ]]; then
    core="${BASH_REMATCH[1]}"
  else
    core="$(echo "$src" | grep -oE '[0-9]+$' || true)"
    [[ -z "$core" ]] && { echo "$src"; return; }
  fi
  local n=$((10#$core))
  echo "sub-$n"
  printf "sub-%03d\n" "$n"
  printf "sub-%04d\n" "$n"
  printf "sub-%05d\n" "$n"
}

find_exact_source(){ # ROOT, ID, SUFFIX -> first exact file match
  local root="$1"
  local id="$2"
  local suffix="$3"
  local idv cand
  while read -r idv; do
    for cand in \
      "$root/${idv}${suffix}.nii.gz" \
      "$root/${idv}${suffix}.nii"; do
      [[ -f "$cand" ]] && { echo "$cand"; return 0; }
    done
    cand="$(find "$root" -maxdepth 2 -type f \( -name "${idv}${suffix}.nii.gz" -o -name "${idv}${suffix}.nii" \) 2>/dev/null | sort | head -n1 || true)"
    [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  done < <(id_variants "$id")
  return 1
}

vox2dils(){ # args: mm voxel_mm -> integer dilations (>=1)
  python3 - "$1" "$2" <<'PY'
import sys
mm=float(sys.argv[1]); vox=float(sys.argv[2])
print(max(1, int(round(mm/vox))))
PY
}

# ---------- parse CLI ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --roi-root)    ROI_ROOT="$2"; shift 2;;
    --t1-root)     T1_ROOT="$2"; shift 2;;
    --lesion-root) LESION_ROOT="$2"; shift 2;;
    --disc-root)   DISC_ROOT="$2"; shift 2;;
    --conn-root)   CONN_ROOT="$2"; shift 2;;
    --out-root)    OUT_ROOT="$2"; shift 2;;
    --threads)     THREADS="$2"; shift 2;;
    --id-map)      IDMAP="$2"; shift 2;;
    -h|--help)     sed -n '1,220p' "$0"; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -n "$T1_ROOT" ]]     || T1_ROOT="$ROI_ROOT/T1"
[[ -n "$LESION_ROOT" ]] || LESION_ROOT="$ROI_ROOT/Lesion"
[[ -n "$DISC_ROOT" ]]   || DISC_ROOT="$ROI_ROOT/Disconnectome"

[[ -d "$CONN_ROOT" ]]   || die "--conn-root not found"
[[ -d "$T1_ROOT" ]]     || die "--t1-root not found"
[[ -d "$LESION_ROOT" ]] || die "--lesion-root not found"
[[ -n "$OUT_ROOT" ]]    || die "--out-root is required"
mkdir -p "$OUT_ROOT"

have antsRegistration || die "antsRegistration not in PATH"
have antsApplyTransforms || die "antsApplyTransforms not in PATH"
have fslmaths || die "fslmaths not in PATH"
have fslval   || die "fslval not in PATH"

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"

# ---------- subjects ----------
mapfile -t CONN_SUBS < <(find "$CONN_ROOT" -maxdepth 1 -type d -name "sub-*" -printf "%f\n" | sort)
[[ ${#CONN_SUBS[@]} -gt 0 ]] || die "No sub-* in $CONN_ROOT"

# manifest
MANIFEST="$OUT_ROOT/manifest_$(date +%Y%m%d_%H%M%S).csv"
echo "conn_id,status,anat_ref,func_ref,src_id,t1_src,lesion_src,disc_src,out_dir" > "$MANIFEST"

for CONN_ID in "${CONN_SUBS[@]}"; do
  log "Subject: $CONN_ID"

  SUBJ_DIR="$CONN_ROOT/$CONN_ID"
  OUT_SUB="$OUT_ROOT/$CONN_ID"
  TX="$OUT_SUB/tx"
  mkdir -p "$OUT_SUB" "$TX"

  # CONN anatomical reference
  ANAT_REF="$(pick_first "$SUBJ_DIR/anat/*_run-01_T1w*.nii.gz" || true)"
  [[ -n "$ANAT_REF" ]] || ANAT_REF="$(pick_first "$SUBJ_DIR/anat/*_run-01_T1w*.nii" || true)"
  [[ -n "$ANAT_REF" ]] || ANAT_REF="$(pick_first "$SUBJ_DIR/anat/*_T1w*.nii.gz" || true)"
  [[ -n "$ANAT_REF" ]] || ANAT_REF="$(pick_first "$SUBJ_DIR/anat/*_T1w*.nii" || true)"
  [[ -n "$ANAT_REF" ]] || { log "  No CONN ANAT T1w - skip"; echo "$CONN_ID,fail_no_anat,,,,,,,$OUT_SUB" >> "$MANIFEST"; continue; }

  # CONN functional reference
  FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_run-01_bold*.nii.gz" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_run-01_bold*.nii" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_bold*.nii.gz" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_bold*.nii" || true)"
  [[ -n "$FUNC_REF" ]] || { log "  No CONN FUNC bold - skip"; echo "$CONN_ID,fail_no_func,$ANAT_REF,,,,,,$OUT_SUB" >> "$MANIFEST"; continue; }

  FUNC_REF3D="$OUT_SUB/func_mean_ref.nii.gz"
  fslmaths "$FUNC_REF" -Tmean "$FUNC_REF3D" >/dev/null

  # conn_id -> source_id
  SRC_ID="$CONN_ID"
  if [[ -n "$IDMAP" && -f "$IDMAP" ]]; then
    MAPPED="$(id_from_csv "$IDMAP" "$CONN_ID" || true)"
    [[ -n "$MAPPED" ]] && SRC_ID="$MAPPED"
  fi

  # source files
  T1_SRC="$(find_exact_source "$T1_ROOT" "$SRC_ID" "_T1" || true)"
  [[ -z "$T1_SRC" ]] && T1_SRC="$(find_exact_source "$T1_ROOT" "$SRC_ID" "_T1w" || true)"
  LES_SRC="$(find_exact_source "$LESION_ROOT" "$SRC_ID" "_lesion" || true)"
  DISC_SRC=""
  if [[ -n "$DISC_ROOT" && -d "$DISC_ROOT" ]]; then
    DISC_SRC="$(find_exact_source "$DISC_ROOT" "$SRC_ID" "_discon" || true)"
  fi

  if [[ -n "$LES_SRC" && -f "$LES_SRC" && -n "$T1_SRC" && -f "$T1_SRC" ]]; then
    status="patient"
    log "  Patient: compute ANAT + FUNC ROIs"

    # Source T1 -> CONN anatomy
    antsRegistration \
      -d 3 \
      -r ["$ANAT_REF","$T1_SRC",1] \
      -m mattes["$ANAT_REF","$T1_SRC",1,64,regular,0.2] \
      -t translation[0.1] -c 1000x500x250 -s 4x2x1 -f 4x2x1 \
      -m mattes["$ANAT_REF","$T1_SRC",1,64,regular,0.2] \
      -t rigid[0.1]       -c 1000x500x250 -s 4x2x1 -f 4x2x1 \
      -m mattes["$ANAT_REF","$T1_SRC",1,64,regular,0.2] \
      -t affine[0.1]      -c 1000x500x250 -s 4x2x1 -f 4x2x1 \
      -m cc["$ANAT_REF","$T1_SRC",1,4] \
      -t SyN[0.1,3,0]     -c 100x70x50   -s 2x1x0 -f 4x2x1 \
      -o "$TX/t1SRC_to_CONNanat_"

    antsApplyTransforms -d 3 \
      -i "$LES_SRC" -r "$ANAT_REF" -n NearestNeighbor \
      -o "$OUT_SUB/lesion_in_CONNanat.nii.gz" \
      -t "$TX/t1SRC_to_CONNanat_1Warp.nii.gz" \
      -t "$TX/t1SRC_to_CONNanat_0GenericAffine.mat"
    fslmaths "$OUT_SUB/lesion_in_CONNanat.nii.gz" -bin "$OUT_SUB/lesion_in_CONNanat.nii.gz"

    if [[ -n "$DISC_SRC" && -f "$DISC_SRC" ]]; then
      antsApplyTransforms -d 3 \
        -i "$DISC_SRC" -r "$ANAT_REF" -n Linear \
        -o "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" \
        -t "$TX/t1SRC_to_CONNanat_1Warp.nii.gz" \
        -t "$TX/t1SRC_to_CONNanat_0GenericAffine.mat"
      for thr in 0.3 0.5 0.7; do
        fslmaths "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" -thr "$thr" -bin \
          "$OUT_SUB/disconnectome_thr${thr}_in_CONNanat.nii.gz"
      done
    fi

    voxmmA=$(fslval "$ANAT_REF" pixdim1)
    nA=$(vox2dils "$RING_MM_A" "$voxmmA")
    nB=$(vox2dils "$RING_MM_B" "$voxmmA")
    nC=$(vox2dils "$RING_MM_C" "$voxmmA")

    fslmaths "$OUT_SUB/lesion_in_CONNanat.nii.gz" -dilM -kernel 3D -dilM "$TX/dilA"
    for i in $(seq 2 "$nA"); do fslmaths "$TX/dilA" -dilM "$TX/dilA"; done
    fslmaths "$TX/dilA" -sub "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_0to2mm_in_CONNanat.nii.gz"
    fslmaths "$OUT_SUB/ring_0to2mm_in_CONNanat.nii.gz" -bin "$OUT_SUB/ring_0to2mm_in_CONNanat.nii.gz"

    if have_imcp; then imcp "$TX/dilA" "$TX/dilB"; else cp "$TX/dilA.nii.gz" "$TX/dilB.nii.gz"; fi
    for i in $(seq $((nA+1)) "$nB"); do fslmaths "$TX/dilB" -dilM "$TX/dilB"; done
    fslmaths "$TX/dilB" -sub "$TX/dilA" "$OUT_SUB/ring_2to4mm_in_CONNanat.nii.gz"
    fslmaths "$OUT_SUB/ring_2to4mm_in_CONNanat.nii.gz" -bin "$OUT_SUB/ring_2to4mm_in_CONNanat.nii.gz"

    if have_imcp; then imcp "$TX/dilB" "$TX/dilC"; else cp "$TX/dilB.nii.gz" "$TX/dilC.nii.gz"; fi
    for i in $(seq $((nB+1)) "$nC"); do fslmaths "$TX/dilC" -dilM "$TX/dilC"; done
    fslmaths "$TX/dilC" -sub "$TX/dilB" "$OUT_SUB/ring_4to8mm_in_CONNanat.nii.gz"
    fslmaths "$OUT_SUB/ring_4to8mm_in_CONNanat.nii.gz" -bin "$OUT_SUB/ring_4to8mm_in_CONNanat.nii.gz"

    rm -f "$TX"/dilA* "$TX"/dilB* "$TX"/dilC* 2>/dev/null || true

    # CONN anatomy -> CONN functional
    antsRegistration \
      -d 3 \
      -r ["$FUNC_REF3D","$ANAT_REF",1] \
      -m mattes["$FUNC_REF3D","$ANAT_REF",1,64,regular,0.2] \
      -t rigid[0.1]  -c 1000x500x250 -s 4x2x1 -f 4x2x1 \
      -m mattes["$FUNC_REF3D","$ANAT_REF",1,64,regular,0.2] \
      -t affine[0.1] -c 1000x500x250 -s 4x2x1 -f 4x2x1 \
      -o "$TX/CONNanat_to_CONNfunc_"

    antsApplyTransforms -d 3 \
      -i "$OUT_SUB/lesion_in_CONNanat.nii.gz" \
      -r "$FUNC_REF3D" -n NearestNeighbor \
      -o "$OUT_SUB/lesion_in_CONNfunc.nii.gz" \
      -t "$TX/CONNanat_to_CONNfunc_0GenericAffine.mat"
    fslmaths "$OUT_SUB/lesion_in_CONNfunc.nii.gz" -bin "$OUT_SUB/lesion_in_CONNfunc.nii.gz"

    if [[ -f "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" ]]; then
      antsApplyTransforms -d 3 \
        -i "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" \
        -r "$FUNC_REF3D" -n Linear \
        -o "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
        -t "$TX/CONNanat_to_CONNfunc_0GenericAffine.mat"
      for thr in 0.3 0.5 0.7; do
        fslmaths "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" -thr "$thr" -bin \
          "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
      done
    else
      fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz"
      for thr in 0.3 0.5 0.7; do
        cp "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
           "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
      done
    fi

    voxmmF=$(fslval "$FUNC_REF3D" pixdim1)
    nAf=$(vox2dils "$RING_MM_A" "$voxmmF")
    nBf=$(vox2dils "$RING_MM_B" "$voxmmF")
    nCf=$(vox2dils "$RING_MM_C" "$voxmmF")

    fslmaths "$OUT_SUB/lesion_in_CONNfunc.nii.gz" -dilM -kernel 3D -dilM "$TX/dilA"
    [[ -f "$TX/dilA.nii.gz" ]] || { log "  Failed to create $TX/dilA.nii.gz"; echo "$CONN_ID,fail_dilate,$ANAT_REF,$FUNC_REF,$SRC_ID,$T1_SRC,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"; continue; }
    for i in $(seq 2 "$nAf"); do fslmaths "$TX/dilA" -dilM "$TX/dilA"; done
    fslmaths "$TX/dilA" -sub "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_0to2mm_in_CONNfunc.nii.gz"
    fslmaths "$OUT_SUB/ring_0to2mm_in_CONNfunc.nii.gz" -bin "$OUT_SUB/ring_0to2mm_in_CONNfunc.nii.gz"

    if have_imcp; then imcp "$TX/dilA" "$TX/dilB"; else cp "$TX/dilA.nii.gz" "$TX/dilB.nii.gz"; fi
    for i in $(seq $((nAf+1)) "$nBf"); do fslmaths "$TX/dilB" -dilM "$TX/dilB"; done
    fslmaths "$TX/dilB" -sub "$TX/dilA" "$OUT_SUB/ring_2to4mm_in_CONNfunc.nii.gz"
    fslmaths "$OUT_SUB/ring_2to4mm_in_CONNfunc.nii.gz" -bin "$OUT_SUB/ring_2to4mm_in_CONNfunc.nii.gz"

    if have_imcp; then imcp "$TX/dilB" "$TX/dilC"; else cp "$TX/dilB.nii.gz" "$TX/dilC.nii.gz"; fi
    for i in $(seq $((nBf+1)) "$nCf"); do fslmaths "$TX/dilC" -dilM "$TX/dilC"; done
    fslmaths "$TX/dilC" -sub "$TX/dilB" "$OUT_SUB/ring_4to8mm_in_CONNfunc.nii.gz"
    fslmaths "$OUT_SUB/ring_4to8mm_in_CONNfunc.nii.gz" -bin "$OUT_SUB/ring_4to8mm_in_CONNfunc.nii.gz"

    rm -f "$TX"/dilA* "$TX"/dilB* "$TX"/dilC* 2>/dev/null || true
  else
    status="control"
    log "  Control: create blanks in ANAT and FUNC"

    fslmaths "$ANAT_REF" -mul 0 "$OUT_SUB/lesion_in_CONNanat.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_0to2mm_in_CONNanat.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_2to4mm_in_CONNanat.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_4to8mm_in_CONNanat.nii.gz"
    fslmaths "$ANAT_REF" -mul 0 "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz"
    for thr in 0.3 0.5 0.7; do
      cp "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" \
         "$OUT_SUB/disconnectome_thr${thr}_in_CONNanat.nii.gz"
    done

    fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/lesion_in_CONNfunc.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_0to2mm_in_CONNfunc.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_2to4mm_in_CONNfunc.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_4to8mm_in_CONNfunc.nii.gz"
    fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz"
    for thr in 0.3 0.5 0.7; do
      cp "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
         "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
    done
  fi

  echo "$CONN_ID,$status,$ANAT_REF,$FUNC_REF,$SRC_ID,$T1_SRC,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"
done

log "All subjects processed."
log "Outputs -> $OUT_ROOT"
log "Manifest -> $MANIFEST"
