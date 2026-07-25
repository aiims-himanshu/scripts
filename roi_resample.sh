#!/usr/bin/env bash
set -euo pipefail

# ╭──────────────────────────────────────────────────────────────╮
# │ roi_resample.sh                                              │
# │ Resample MNI152 lesion/disconnectome ROIs into CONN         │
# │ MNI152NLin2009cAsym anat and func spaces using one bridge   │
# ╰──────────────────────────────────────────────────────────────╯
#
# Workflow:
#   source ROIs in MNI152
#   -> reusable affine bridge to MNI152NLin2009cAsym
#   -> resample directly into CONN anat grid
#   -> resample directly into CONN func grid
#   -> build peri-lesional rings in each target grid
#
# Default layout:
#   ROI_ROOT/
#     Lesion/            sub-XXX_lesion.nii.gz
#     Disconnectome/     sub-XXX_discon.nii.gz
#
#   CONN_ROOT/
#     sub-*/anat/*_T1w*.nii[.gz]
#     sub-*/func/*_bold*.nii[.gz]
#
# Example:
#   bash roi_resample.sh \
#     --roi-root  /media/sskgroup/Thesis/Primary/roi \
#     --conn-root /media/sskgroup/Thesis/Primary/trial/conn_project01/data/BIDS/dataset \
#     --out-root  /media/sskgroup/Thesis/Primary/trial/conn_rois_mni_bridge \
#     --bridge-affine /home/sskgroup/Git/MNI152/bridge/MNI152_to_MNI152NLin2009cAsym_0GenericAffine.mat
#
# Optional:
#   --id-map /path/to/idmap.csv   # CSV columns: conn_id,source_id
#
# Requirements: antsApplyTransforms, FSL (fslmaths, fslval, imcp optional)

# ---------- defaults ----------
ROI_ROOT="/media/sskgroup/Thesis/Primary/trial/roi"
LESION_ROOT=""
DISC_ROOT=""
CONN_ROOT="/media/sskgroup/Thesis/Primary/trial/conn_project01/data/BIDS/dataset"
OUT_ROOT="/media/sskgroup/Thesis/Primary/trial/conn_rois_mni_bridge"
BRIDGE_AFFINE="/home/sskgroup/Git/MNI152/bridge/MNI152_to_MNI152NLin2009cAsym_0GenericAffine.mat"
THREADS="${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS:-16}"
IDMAP=""

# peri-lesional ring radii (mm)
RING_MM_A=3
RING_MM_B=6
RING_MM_C=12

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

resample_nii(){ # src ref interp out [bridge_affine]
  local src="$1"
  local ref="$2"
  local interp="$3"
  local out="$4"
  local xf="${5:-}"
  if [[ -n "$xf" ]]; then
    antsApplyTransforms -d 3 \
      -i "$src" -r "$ref" -n "$interp" \
      -o "$out" \
      -t "$xf"
  else
    antsApplyTransforms -d 3 \
      -i "$src" -r "$ref" -n "$interp" \
      -o "$out"
  fi
}

voxcount(){ # binary/thresholded nifti -> voxel count
  local file="$1"
  fslstats "$file" -V | awk '{print $1+0}'
}

# ---------- parse CLI ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --roi-root)      ROI_ROOT="$2"; shift 2;;
    --lesion-root)   LESION_ROOT="$2"; shift 2;;
    --disc-root)     DISC_ROOT="$2"; shift 2;;
    --conn-root)     CONN_ROOT="$2"; shift 2;;
    --out-root)      OUT_ROOT="$2"; shift 2;;
    --bridge-affine) BRIDGE_AFFINE="$2"; shift 2;;
    --threads)       THREADS="$2"; shift 2;;
    --id-map)        IDMAP="$2"; shift 2;;
    -h|--help)
      sed -n '1,220p' "$0"
      exit 0
      ;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -n "$LESION_ROOT" ]] || LESION_ROOT="$ROI_ROOT/Lesion"
[[ -n "$DISC_ROOT" ]]   || DISC_ROOT="$ROI_ROOT/Disconnectome"

[[ -d "$CONN_ROOT" ]]   || die "--conn-root not found"
[[ -d "$LESION_ROOT" ]] || die "--lesion-root not found"
[[ -n "$OUT_ROOT" ]]    || die "--out-root is required"
[[ -f "$BRIDGE_AFFINE" ]] || die "--bridge-affine not found: $BRIDGE_AFFINE"
mkdir -p "$OUT_ROOT"

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
echo "conn_id,status,anat_ref,func_ref,src_id,lesion_src,disc_src,out_dir" > "$MANIFEST"
QC_TSV="$OUT_ROOT/qc_counts_$(date +%Y%m%d_%H%M%S).tsv"
cat > "$QC_TSV" <<'EOF'
conn_id	status	src_id	lesion_src_vox	lesion_anat_vox	lesion_func_vox	disc_src_vox	disc_src_thr0.3_vox	disc_src_thr0.5_vox	disc_src_thr0.7_vox	disc_anat_thr0.3_vox	disc_anat_thr0.5_vox	disc_anat_thr0.7_vox	disc_func_thr0.3_vox	disc_func_thr0.5_vox	disc_func_thr0.7_vox	ring0to3_anat_vox	ring3to6_anat_vox	ring6to12_anat_vox	ring0to3_func_vox	ring3to6_func_vox	ring6to12_func_vox	anat_ref	func_ref	lesion_src	disc_src	out_dir
EOF

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
  [[ -n "$ANAT_REF" ]] || { log "  No CONN ANAT T1w - skip"; echo "$CONN_ID,fail_no_anat,,,,,,$OUT_SUB" >> "$MANIFEST"; continue; }

  # CONN functional reference
  FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_run-01_bold*.nii.gz" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_run-01_bold*.nii" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_bold*.nii.gz" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_bold*.nii" || true)"
  [[ -n "$FUNC_REF" ]] || { log "  No CONN FUNC bold - skip"; echo "$CONN_ID,fail_no_func,$ANAT_REF,,,,,$OUT_SUB" >> "$MANIFEST"; continue; }

  FUNC_REF3D="$OUT_SUB/func_mean_ref.nii.gz"
  fslmaths "$FUNC_REF" -Tmean "$FUNC_REF3D" >/dev/null

  # conn_id -> source_id
  SRC_ID="$CONN_ID"
  if [[ -n "$IDMAP" && -f "$IDMAP" ]]; then
    MAPPED="$(id_from_csv "$IDMAP" "$CONN_ID" || true)"
    [[ -n "$MAPPED" ]] && SRC_ID="$MAPPED"
  fi

  # source files in MNI152 space
  LES_SRC="$(find_exact_source "$LESION_ROOT" "$SRC_ID" "_lesion" || true)"
  DISC_SRC=""
  if [[ -n "$DISC_ROOT" && -d "$DISC_ROOT" ]]; then
    DISC_SRC="$(find_exact_source "$DISC_ROOT" "$SRC_ID" "_discon" || true)"
  fi

  if [[ -n "$LES_SRC" && -f "$LES_SRC" ]]; then
    status="patient"
    log "  Patient: bridge + resample to ANAT/FUNC"
    LES_SRC_VOX=$(voxcount "$LES_SRC")
    DISC_SRC_VOX=0
    DISC_SRC_T03=0
    DISC_SRC_T05=0
    DISC_SRC_T07=0

    # Source lesion -> CONN ANAT/FUNC
    resample_nii "$LES_SRC" "$ANAT_REF" NearestNeighbor "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$BRIDGE_AFFINE"
    fslmaths "$OUT_SUB/lesion_in_CONNanat.nii.gz" -bin "$OUT_SUB/lesion_in_CONNanat.nii.gz"

    resample_nii "$LES_SRC" "$FUNC_REF3D" NearestNeighbor "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$BRIDGE_AFFINE"
    fslmaths "$OUT_SUB/lesion_in_CONNfunc.nii.gz" -bin "$OUT_SUB/lesion_in_CONNfunc.nii.gz"
    LES_ANAT_VOX=$(voxcount "$OUT_SUB/lesion_in_CONNanat.nii.gz")
    LES_FUNC_VOX=$(voxcount "$OUT_SUB/lesion_in_CONNfunc.nii.gz")

    # Disconnectome probability maps
    if [[ -n "$DISC_SRC" && -f "$DISC_SRC" ]]; then
      DISC_SRC_VOX=$(voxcount "$DISC_SRC")
      fslmaths "$DISC_SRC" -thr 0.3 -bin "$TX/disc_src_thr0p3.nii.gz"
      fslmaths "$DISC_SRC" -thr 0.5 -bin "$TX/disc_src_thr0p5.nii.gz"
      fslmaths "$DISC_SRC" -thr 0.7 -bin "$TX/disc_src_thr0p7.nii.gz"
      DISC_SRC_T03=$(voxcount "$TX/disc_src_thr0p3.nii.gz")
      DISC_SRC_T05=$(voxcount "$TX/disc_src_thr0p5.nii.gz")
      DISC_SRC_T07=$(voxcount "$TX/disc_src_thr0p7.nii.gz")

      resample_nii "$DISC_SRC" "$ANAT_REF" Linear "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" "$BRIDGE_AFFINE"
      for thr in 0.3 0.5 0.7; do
        fslmaths "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" -thr "$thr" -bin \
          "$OUT_SUB/disconnectome_thr${thr}_in_CONNanat.nii.gz"
      done
      DISC_ANAT_T03=$(voxcount "$OUT_SUB/disconnectome_thr0.3_in_CONNanat.nii.gz")
      DISC_ANAT_T05=$(voxcount "$OUT_SUB/disconnectome_thr0.5_in_CONNanat.nii.gz")
      DISC_ANAT_T07=$(voxcount "$OUT_SUB/disconnectome_thr0.7_in_CONNanat.nii.gz")

      resample_nii "$DISC_SRC" "$FUNC_REF3D" Linear "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" "$BRIDGE_AFFINE"
      for thr in 0.3 0.5 0.7; do
        fslmaths "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" -thr "$thr" -bin \
          "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
      done
      DISC_FUNC_T03=$(voxcount "$OUT_SUB/disconnectome_thr0.3_in_CONNfunc.nii.gz")
      DISC_FUNC_T05=$(voxcount "$OUT_SUB/disconnectome_thr0.5_in_CONNfunc.nii.gz")
      DISC_FUNC_T07=$(voxcount "$OUT_SUB/disconnectome_thr0.7_in_CONNfunc.nii.gz")
    else
      DISC_SRC_VOX=0
      DISC_SRC_T03=0
      DISC_SRC_T05=0
      DISC_SRC_T07=0
      fslmaths "$ANAT_REF" -mul 0 "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz"
      for thr in 0.3 0.5 0.7; do
        cp "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" \
           "$OUT_SUB/disconnectome_thr${thr}_in_CONNanat.nii.gz"
      done
      DISC_ANAT_T03=0
      DISC_ANAT_T05=0
      DISC_ANAT_T07=0

      fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz"
      for thr in 0.3 0.5 0.7; do
        cp "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
           "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
      done
      DISC_FUNC_T03=0
      DISC_FUNC_T05=0
      DISC_FUNC_T07=0
    fi

    # Rings in ANAT space
    voxmmA=$(fslval "$ANAT_REF" pixdim1)
    nA=$(vox2dils "$RING_MM_A" "$voxmmA")
    nB=$(vox2dils "$RING_MM_B" "$voxmmA")
    nC=$(vox2dils "$RING_MM_C" "$voxmmA")

    fslmaths "$OUT_SUB/lesion_in_CONNanat.nii.gz" -dilM -kernel 3D -dilM "$TX/dilA"
    for i in $(seq 2 "$nA"); do fslmaths "$TX/dilA" -dilM "$TX/dilA"; done
    fslmaths "$TX/dilA" -sub "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_0to3mm_in_CONNanat.nii.gz"
    fslmaths "$OUT_SUB/ring_0to3mm_in_CONNanat.nii.gz" -bin "$OUT_SUB/ring_0to3mm_in_CONNanat.nii.gz"

    if have_imcp; then imcp "$TX/dilA" "$TX/dilB"; else cp "$TX/dilA.nii.gz" "$TX/dilB.nii.gz"; fi
    for i in $(seq $((nA+1)) "$nB"); do fslmaths "$TX/dilB" -dilM "$TX/dilB"; done
    fslmaths "$TX/dilB" -sub "$TX/dilA" "$OUT_SUB/ring_3to6mm_in_CONNanat.nii.gz"
    fslmaths "$OUT_SUB/ring_3to6mm_in_CONNanat.nii.gz" -bin "$OUT_SUB/ring_3to6mm_in_CONNanat.nii.gz"

    if have_imcp; then imcp "$TX/dilB" "$TX/dilC"; else cp "$TX/dilB.nii.gz" "$TX/dilC.nii.gz"; fi
    for i in $(seq $((nB+1)) "$nC"); do fslmaths "$TX/dilC" -dilM "$TX/dilC"; done
    fslmaths "$TX/dilC" -sub "$TX/dilB" "$OUT_SUB/ring_6to12mm_in_CONNanat.nii.gz"
    fslmaths "$OUT_SUB/ring_6to12mm_in_CONNanat.nii.gz" -bin "$OUT_SUB/ring_6to12mm_in_CONNanat.nii.gz"

    rm -f "$TX"/dilA* "$TX"/dilB* "$TX"/dilC* 2>/dev/null || true

    RING0A=$(voxcount "$OUT_SUB/ring_0to3mm_in_CONNanat.nii.gz")
    RING3A=$(voxcount "$OUT_SUB/ring_3to6mm_in_CONNanat.nii.gz")
    RING6A=$(voxcount "$OUT_SUB/ring_6to12mm_in_CONNanat.nii.gz")

    # Rings in FUNC space
    voxmmF=$(fslval "$FUNC_REF3D" pixdim1)
    nAf=$(vox2dils "$RING_MM_A" "$voxmmF")
    nBf=$(vox2dils "$RING_MM_B" "$voxmmF")
    nCf=$(vox2dils "$RING_MM_C" "$voxmmF")

    fslmaths "$OUT_SUB/lesion_in_CONNfunc.nii.gz" -dilM -kernel 3D -dilM "$TX/dilA"
    [[ -f "$TX/dilA.nii.gz" ]] || { log "  Failed to create $TX/dilA.nii.gz"; echo "$CONN_ID,fail_dilate,$ANAT_REF,$FUNC_REF,$SRC_ID,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"; continue; }
    for i in $(seq 2 "$nAf"); do fslmaths "$TX/dilA" -dilM "$TX/dilA"; done
    fslmaths "$TX/dilA" -sub "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_0to3mm_in_CONNfunc.nii.gz"
    fslmaths "$OUT_SUB/ring_0to3mm_in_CONNfunc.nii.gz" -bin "$OUT_SUB/ring_0to3mm_in_CONNfunc.nii.gz"

    if have_imcp; then imcp "$TX/dilA" "$TX/dilB"; else cp "$TX/dilA.nii.gz" "$TX/dilB.nii.gz"; fi
    for i in $(seq $((nAf+1)) "$nBf"); do fslmaths "$TX/dilB" -dilM "$TX/dilB"; done
    fslmaths "$TX/dilB" -sub "$TX/dilA" "$OUT_SUB/ring_3to6mm_in_CONNfunc.nii.gz"
    fslmaths "$OUT_SUB/ring_3to6mm_in_CONNfunc.nii.gz" -bin "$OUT_SUB/ring_3to6mm_in_CONNfunc.nii.gz"

    if have_imcp; then imcp "$TX/dilB" "$TX/dilC"; else cp "$TX/dilB.nii.gz" "$TX/dilC.nii.gz"; fi
    for i in $(seq $((nBf+1)) "$nCf"); do fslmaths "$TX/dilC" -dilM "$TX/dilC"; done
    fslmaths "$TX/dilC" -sub "$TX/dilB" "$OUT_SUB/ring_6to12mm_in_CONNfunc.nii.gz"
    fslmaths "$OUT_SUB/ring_6to12mm_in_CONNfunc.nii.gz" -bin "$OUT_SUB/ring_6to12mm_in_CONNfunc.nii.gz"

    rm -f "$TX"/dilA* "$TX"/dilB* "$TX"/dilC* 2>/dev/null || true

    RING0F=$(voxcount "$OUT_SUB/ring_0to3mm_in_CONNfunc.nii.gz")
    RING3F=$(voxcount "$OUT_SUB/ring_3to6mm_in_CONNfunc.nii.gz")
    RING6F=$(voxcount "$OUT_SUB/ring_6to12mm_in_CONNfunc.nii.gz")
  else
    status="control"
    log "  Control: create blanks in ANAT and FUNC"

    fslmaths "$ANAT_REF" -mul 0 "$OUT_SUB/lesion_in_CONNanat.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_0to3mm_in_CONNanat.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_3to6mm_in_CONNanat.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNanat.nii.gz" "$OUT_SUB/ring_6to12mm_in_CONNanat.nii.gz"
    fslmaths "$ANAT_REF" -mul 0 "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz"
    for thr in 0.3 0.5 0.7; do
      cp "$OUT_SUB/disconnectome_in_CONNanat_prob.nii.gz" \
         "$OUT_SUB/disconnectome_thr${thr}_in_CONNanat.nii.gz"
    done

    fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/lesion_in_CONNfunc.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_0to3mm_in_CONNfunc.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_3to6mm_in_CONNfunc.nii.gz"
    cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_6to12mm_in_CONNfunc.nii.gz"
    fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz"
    for thr in 0.3 0.5 0.7; do
      cp "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
         "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
    done

    LES_SRC_VOX=0
    LES_ANAT_VOX=0
    LES_FUNC_VOX=0
    DISC_SRC_VOX=0
    DISC_SRC_T03=0
    DISC_SRC_T05=0
    DISC_SRC_T07=0
    DISC_ANAT_T03=0
    DISC_ANAT_T05=0
    DISC_ANAT_T07=0
    DISC_FUNC_T03=0
    DISC_FUNC_T05=0
    DISC_FUNC_T07=0
    RING0A=0
    RING3A=0
    RING6A=0
    RING0F=0
    RING3F=0
    RING6F=0
  fi

  if [[ "$status" == "patient" ]]; then
    echo "[`date +%F' '%T`]   lesion voxels: src=$LES_SRC_VOX anat=$LES_ANAT_VOX func=$LES_FUNC_VOX"
    echo "[`date +%F' '%T`]   discon thr voxels: src(0.3/0.5/0.7)=$DISC_SRC_T03/$DISC_SRC_T05/$DISC_SRC_T07 anat=$DISC_ANAT_T03/$DISC_ANAT_T05/$DISC_ANAT_T07 func=$DISC_FUNC_T03/$DISC_FUNC_T05/$DISC_FUNC_T07"
    echo "[`date +%F' '%T`]   rings voxels anat(0-3/3-6/6-12)=$RING0A/$RING3A/$RING6A func=$RING0F/$RING3F/$RING6F"
  fi

  echo -e "${CONN_ID}\t${status}\t${SRC_ID}\t${LES_SRC_VOX}\t${LES_ANAT_VOX}\t${LES_FUNC_VOX}\t${DISC_SRC_VOX}\t${DISC_SRC_T03}\t${DISC_SRC_T05}\t${DISC_SRC_T07}\t${DISC_ANAT_T03}\t${DISC_ANAT_T05}\t${DISC_ANAT_T07}\t${DISC_FUNC_T03}\t${DISC_FUNC_T05}\t${DISC_FUNC_T07}\t${RING0A}\t${RING3A}\t${RING6A}\t${RING0F}\t${RING3F}\t${RING6F}\t${ANAT_REF}\t${FUNC_REF}\t${LES_SRC}\t${DISC_SRC}\t${OUT_SUB}" >> "$QC_TSV"
  echo "$CONN_ID,$status,$ANAT_REF,$FUNC_REF,$SRC_ID,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"
done

log "All subjects processed."
log "Outputs -> $OUT_ROOT"
log "Manifest -> $MANIFEST"
log "QC TSV -> $QC_TSV"
