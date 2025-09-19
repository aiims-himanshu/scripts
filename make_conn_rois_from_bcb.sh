#!/usr/bin/env bash
set -euo pipefail

# ╭──────────────────────────────────────────────────────────────╮
# │ make_conn_rois_from_bcb.sh                                   │
# │ Lesion/Perilesional/Disconnectome ROIs in CONN space         │
# │ • New: "Repair mode" fixes ANAT→FUNC for existing patients   │
# │ • Skips finished controls                                    │
# ╰──────────────────────────────────────────────────────────────╯
# Example:
#   bash make_conn_rois_from_bcb.sh \
#     --conn-root   /path/to/CONN_dataset_root \
#     --t1-root     /path/to/T1 \
#     --lesion-root /path/to/Lesion \
#     --disc-root   /path/to/Disconnectome \
#     --out-root    /path/to/derivatives/conn_rois \
#     --threads 16
#
# CONN dataset layout matches your tree (anat/*_T1w*, func/*_run-01_bold*). :contentReference[oaicite:3]{index=3}
# Lesion/Disconnectome names match your sources (e.g., sub-52_lesion_MNI.nii.gz). :contentReference[oaicite:4]{index=4}
# This builds on your previous script with fixes (skip, ID variants, imcp, 3D func ref). :contentReference[oaicite:5]{index=5}

# ---------- defaults (override with CLI) ----------
CONN_ROOT="/media/sskgroup/Thesis/functional/primary/data/BIDS/dataset"
T1_ROOT="/media/sskgroup/Thesis/functional/Primary_Tractron/T1"
LESION_ROOT="/media/sskgroup/Thesis/functional/Primary_Tractron/Lesion"
DISC_ROOT="/media/sskgroup/Thesis/functional/Primary_Tractron/Disconnectome"   # optional
OUT_ROOT="/media/sskgroup/Thesis/functional/Primary_Tractron/derivatives/conn_rois"
THREADS=16
IDMAP=""   # optional CSV: conn_id,source_id

# ring radii in mm
RING_MM_A=2
RING_MM_B=4
RING_MM_C=8

# ---------- helpers ----------
die(){ echo "ERROR: $*" >&2; exit 1; }
log(){ echo "[`date +%F' '%T`] $*"; }

have_imcp(){ command -v imcp >/dev/null 2>&1; }

vox2dils(){ # args: mm voxel_mm  -> integer dilations (≥1)
  python3 - "$1" "$2" <<'PY'
import sys, math
mm=float(sys.argv[1]); vox=float(sys.argv[2])
print(max(1, int(round(mm/vox))))
PY
}

pick_first(){ # glob and echo first existing match
  local pattern="$1"; shopt -s nullglob
  for f in $pattern; do echo "$f"; return 0; done
  return 1
}

id_from_csv(){ # map conn_id → source_id via CSV "conn_id,source_id"
  local csv="$1"; local key="$2"
  awk -F, -v k="$key" 'NR>1 && $1==k {print $2; exit}' "$csv"
}

# emit variants for sub-51 ↔ sub-0051 without colliding with sub-510
id_variants(){
  local src="$1" core
  if [[ "$src" =~ ^sub-([0-9]+)$ ]]; then core="${BASH_REMATCH[1]}"
  else core="$(echo "$src" | grep -oE '[0-9]+$' || true)"; [[ -z "$core" ]] && { echo "$src"; return; }
  fi
  local n=$((10#$core))
  echo "sub-$n"
  printf "sub-%03d\n" "$n"
  printf "sub-%04d\n" "$n"
  printf "sub-%05d\n" "$n"
}

find_with_variants(){ # ROOT, ID, MID, SUFFIX -> first match
  local root="$1"; local id="$2"; local mid="$3"; local tail="$4"; local cand idv
  while read -r idv; do
    cand="$(pick_first "$root/${idv}${mid}${tail}" || true)"
    [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  done < <(id_variants "$id")
  return 1
}

already_done_func(){ # out_sub -> 0 if minimal functional outputs exist (skip)
  local out="$1"
  [[ -f "$out/lesion_in_CONNfunc.nii.gz" ]] && return 0
  [[ -f "$out/lesion_in_CONNfunc.nii"    ]] && return 0
  return 1
}

same_grid(){ # a b -> 0 if dim/pixdim 1..3 match exactly
  local a="$1" b="$2" d va vb pa pb
  for d in 1 2 3; do va=$(fslval "$a" dim$d); vb=$(fslval "$b" dim$d); [[ "$va" == "$vb" ]] || return 1; done
  for d in 1 2 3; do pa=$(fslval "$a" pixdim$d); pb=$(fslval "$b" pixdim$d); [[ "$pa" == "$pb" ]] || return 1; done
  return 0
}

needs_repair_func(){ # OUT_SUB FUNC_REF3D -> 0 if OK, 1 if needs repair
  local out="$1" func3d="$2"
  local les="$out/lesion_in_CONNfunc.nii.gz"
  [[ -f "$les" ]] || les="$out/lesion_in_CONNfunc.nii"
  [[ -f "$les" ]] || return 0  # nothing to check; treat as "new run"
  # 1) same grid?
  if ! same_grid "$les" "$func3d"; then return 1; fi
  # 2) any nonzero voxels?
  local nz; nz=$(fslstats "$les" -V | awk '{print $1}')
  [[ "${nz:-0}" -gt 0 ]] || return 1
  # 3) overlap with a functional brain mask?
  local m="$out/tx/funcmask"; mkdir -p "$out/tx"
  fslmaths "$func3d" -bin -dilF "$m" >/dev/null
  local ov; ov=$(fslstats "$les" -k "${m}.nii.gz" -V | awk '{print $1}')
  [[ "${ov:-0}" -gt 0 ]] || return 1
  return 0
}

# ---------- parse CLI ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --conn-root)   CONN_ROOT="$2"; shift 2;;
    --t1-root)     T1_ROOT="$2"; shift 2;;
    --lesion-root) LESION_ROOT="$2"; shift 2;;
    --disc-root)   DISC_ROOT="$2"; shift 2;;
    --out-root)    OUT_ROOT="$2"; shift 2;;
    --threads)     THREADS="$2"; shift 2;;
    --id-map)      IDMAP="$2"; shift 2;;
    -h|--help) sed -n '1,200p' "$0"; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -d "$CONN_ROOT" ]]   || die "--conn-root not found"
[[ -d "$T1_ROOT" ]]     || die "--t1-root not found"
[[ -d "$LESION_ROOT" ]] || die "--lesion-root not found"
mkdir -p "$OUT_ROOT"

command -v antsRegistration >/dev/null || die "antsRegistration not in PATH"
command -v antsApplyTransforms >/dev/null || die "antsApplyTransforms not in PATH"
command -v fslmaths >/dev/null || die "fslmaths not in PATH"
command -v fslval >/dev/null   || die "fslval not in PATH"

# ---------- discover subjects from CONN dataset ----------
mapfile -t CONN_SUBS < <(find "$CONN_ROOT" -maxdepth 1 -type d -name "sub-*" -printf "%f\n" | sort)
[[ ${#CONN_SUBS[@]} -gt 0 ]] || die "No sub-* in $CONN_ROOT"

# manifest
MANIFEST="$OUT_ROOT/manifest_$(date +%Y%m%d_%H%M%S).csv"
echo "conn_id,mode,status,anat_ref,func_ref,src_id,t1_src,lesion_src,disc_src,out_dir" > "$MANIFEST"

# ---------- main loop ----------
for CONN_ID in "${CONN_SUBS[@]}"; do
  log "Subject: $CONN_ID"
  SUBJ_DIR="$CONN_ROOT/$CONN_ID"
  OUT_SUB="$OUT_ROOT/$CONN_ID"
  TX="$OUT_SUB/tx"; mkdir -p "$OUT_SUB" "$TX"

  # CONN references (from your dataset tree) :contentReference[oaicite:6]{index=6}
  ANAT_REF="$(pick_first "$SUBJ_DIR/anat/*_T1w*.nii.gz" || true)"
  [[ -n "$ANAT_REF" ]] || ANAT_REF="$(pick_first "$SUBJ_DIR/anat/*_T1w*.nii" || true)"
  [[ -n "$ANAT_REF" ]] || { log "❌ $CONN_ID no CONN ANAT T1w"; echo "$CONN_ID,,fail_no_anat,,,,,,,$OUT_SUB" >> "$MANIFEST"; continue; }

  FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_run-01_bold*.nii.gz" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_run-01_bold*.nii" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_bold*.nii.gz" || true)"
  [[ -n "$FUNC_REF" ]] || FUNC_REF="$(pick_first "$SUBJ_DIR/func/*_bold*.nii" || true)"
  [[ -n "$FUNC_REF" ]] || { log "❌ $CONN_ID no CONN FUNC bold"; echo "$CONN_ID,,fail_no_func,$ANAT_REF,,,,,,$OUT_SUB" >> "$MANIFEST"; continue; }

  # 3D functional reference (mean of time)
  FUNC_REF3D="$OUT_SUB/func_mean_ref.nii.gz"
  if [[ ! -f "$FUNC_REF3D" ]]; then
    fslmaths "$FUNC_REF" -Tmean "$FUNC_REF3D"
  fi

  # conn_id → source_id (if mapping provided)
  SRC_ID="$CONN_ID"
  if [[ -n "$IDMAP" ]]; then
    MAPPED="$(id_from_csv "$IDMAP" "$CONN_ID" || true)"; [[ -n "$MAPPED" ]] && SRC_ID="$MAPPED"
  fi

  # source files (variants; require '_' after ID to avoid sub-510 collisions)
  T1_SRC="$(find_with_variants "$T1_ROOT"     "$SRC_ID" "_T1_"     "*.nii.gz" || true)"
  [[ -z "$T1_SRC" ]] && T1_SRC="$(find_with_variants "$T1_ROOT"    "$SRC_ID" "_T1"      "*.nii.gz" || true)"
  [[ -z "$T1_SRC" ]] && T1_SRC="$(find_with_variants "$T1_ROOT"    "$SRC_ID" "_T1w"     "*.nii.gz" || true)"

  LES_SRC="$(find_with_variants "$LESION_ROOT" "$SRC_ID" "_lesion_" "*.nii.gz" || true)"

  DISC_SRC=""
  if [[ -n "$DISC_ROOT" && -d "$DISC_ROOT" ]]; then
    DISC_SRC="$(find_with_variants "$DISC_ROOT" "$SRC_ID" "_lesion_" "*.nii.gz" | head -n1 || true)"
  fi

  # --- REPAIR MODE: if lesion_in_CONNanat exists, only fix functional side if needed ---
  LES_ANAT="$OUT_SUB/lesion_in_CONNanat.nii.gz"
  if [[ -f "$LES_ANAT" ]]; then
    mode="repair"
    if needs_repair_func "$OUT_SUB" "$FUNC_REF3D"; then
      log "  → Repairing functional deliverables (ANAT→FUNC alignment)"
      # Estimate CONN ANAT → CONN FUNC (rigid/affine)
      antsRegistration \
        -d 3 \
        -r ["$FUNC_REF3D","$ANAT_REF",1] \
        -m mattes["$FUNC_REF3D","$ANAT_REF",1,64,regular,0.2] \
        -t rigid[0.1]  -c 1000x500x250 -s 4x2x1 -f 4x2x1 \
        -m mattes["$FUNC_REF3D","$ANAT_REF",1,64,regular,0.2] \
        -t affine[0.1] -c 1000x500x250 -s 4x2x1 -f 4x2x1 \
        -o "$TX/CONNanat_to_CONNfunc_"

      # lesion (ANAT) → FUNC
      antsApplyTransforms -d 3 \
        -i "$LES_ANAT" \
        -r "$FUNC_REF3D" -n NearestNeighbor \
        -o "$OUT_SUB/lesion_in_CONNfunc.nii.gz" \
        -t "$TX/CONNanat_to_CONNfunc_0GenericAffine.mat"
      fslmaths "$OUT_SUB/lesion_in_CONNfunc.nii.gz" -bin "$OUT_SUB/lesion_in_CONNfunc.nii.gz"

      # disconnectome (ANAT) → FUNC (prob) if present
      if [[ -f "$OUT_SUB/disconnectome_in_CONNanat.nii.gz" ]]; then
        antsApplyTransforms -d 3 \
          -i "$OUT_SUB/disconnectome_in_CONNanat.nii.gz" \
          -r "$FUNC_REF3D" -n Linear \
          -o "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
          -t "$TX/CONNanat_to_CONNfunc_0GenericAffine.mat"
        for thr in 0.3 0.5 0.7; do
          fslmaths "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" -thr $thr -bin \
            "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
        done
      fi

      # rings in corrected FUNC space
      voxmm=$(fslval "$FUNC_REF3D" pixdim1)
      nA=$(vox2dils "$RING_MM_A" "$voxmm")
      nB=$(vox2dils "$RING_MM_B" "$voxmm")
      nC=$(vox2dils "$RING_MM_C" "$voxmm")

      fslmaths "$OUT_SUB/lesion_in_CONNfunc.nii.gz" -dilM -kernel 3D -dilM "$TX/dilA"
      [[ -f "$TX/dilA.nii.gz" ]] || { log "  ❌ failed to create $TX/dilA.nii.gz"; echo "$CONN_ID,$mode,fail_dilate,$ANAT_REF,$FUNC_REF,$SRC_ID,$T1_SRC,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"; continue; }
      for i in $(seq 2 $nA); do fslmaths "$TX/dilA" -dilM "$TX/dilA"; done
      fslmaths "$TX/dilA" -sub "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_0to2mm.nii.gz"
      fslmaths "$OUT_SUB/ring_0to2mm.nii.gz" -bin "$OUT_SUB/ring_0to2mm.nii.gz"

      if have_imcp; then imcp "$TX/dilA" "$TX/dilB"; else cp "$TX/dilA.nii.gz" "$TX/dilB.nii.gz"; fi
      for i in $(seq $((nA+1)) $nB); do fslmaths "$TX/dilB" -dilM "$TX/dilB"; done
      fslmaths "$TX/dilB" -sub "$TX/dilA" "$OUT_SUB/ring_2to4mm.nii.gz"
      fslmaths "$OUT_SUB/ring_2to4mm.nii.gz" -bin "$OUT_SUB/ring_2to4mm.nii.gz"

      if have_imcp; then imcp "$TX/dilB" "$TX/dilC"; else cp "$TX/dilB.nii.gz" "$TX/dilC.nii.gz"; fi
      for i in $(seq $((nB+1)) $nC); do fslmaths "$TX/dilC" -dilM "$TX/dilC"; done
      fslmaths "$TX/dilC" -sub "$TX/dilB" "$OUT_SUB/ring_4to8mm.nii.gz"
      fslmaths "$OUT_SUB/ring_4to8mm.nii.gz" -bin "$OUT_SUB/ring_4to8mm.nii.gz"

      rm -f "$TX"/dilA* "$TX"/dilB* "$TX"/dilC* 2>/dev/null || true
      status="repaired"
    else
      log "  ✓ Functional deliverables already aligned — no repair needed"
      status="repair_not_needed"
    fi
    echo "$CONN_ID,$mode,$status,$ANAT_REF,$FUNC_REF,$SRC_ID,$T1_SRC,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"
    continue
  fi

  # --- NORMAL MODE (no ANAT outputs yet): produce from sources ---
  mode="normal"

  # Patient vs control decision by presence of lesion+T1
  if [[ -n "$LES_SRC" && -f "$LES_SRC" && -n "$T1_SRC" && -f "$T1_SRC" ]]; then
    log "  → PATIENT (lesion+T1 found)"

    # T1 (source) → CONN ANAT (nonlinear)
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

    # lesion (source) → CONN ANAT (NN)
    antsApplyTransforms -d 3 \
      -i "$LES_SRC" -r "$ANAT_REF" -n NearestNeighbor \
      -o "$OUT_SUB/lesion_in_CONNanat.nii.gz" \
      -t "$TX/t1SRC_to_CONNanat_1Warp.nii.gz" \
      -t "$TX/t1SRC_to_CONNanat_0GenericAffine.mat"

    # disconnectome (source) → CONN ANAT (prob), if present
    if [[ -n "$DISC_SRC" && -f "$DISC_SRC" ]]; then
      antsApplyTransforms -d 3 \
        -i "$DISC_SRC" -r "$ANAT_REF" -n Linear \
        -o "$OUT_SUB/disconnectome_in_CONNanat.nii.gz" \
        -t "$TX/t1SRC_to_CONNanat_1Warp.nii.gz" \
        -t "$TX/t1SRC_to_CONNanat_0GenericAffine.mat"
    fi

    # Now compute CONN ANAT → CONN FUNC (rigid/affine) and finish FUNC side
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

    if [[ -f "$OUT_SUB/disconnectome_in_CONNanat.nii.gz" ]]; then
      antsApplyTransforms -d 3 \
        -i "$OUT_SUB/disconnectome_in_CONNanat.nii.gz" \
        -r "$FUNC_REF3D" -n Linear \
        -o "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
        -t "$TX/CONNanat_to_CONNfunc_0GenericAffine.mat"
      for thr in 0.3 0.5 0.7; do
        fslmaths "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" -thr $thr -bin \
          "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
      done
    else
      # No disconnectome provided → patient gets blanks too
      fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz"
      for thr in 0.3 0.5 0.7; do
        cp "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
           "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
      done
    fi

    # Rings in FUNC space
    voxmm=$(fslval "$FUNC_REF3D" pixdim1)
    nA=$(vox2dils "$RING_MM_A" "$voxmm")
    nB=$(vox2dils "$RING_MM_B" "$voxmm")
    nC=$(vox2dils "$RING_MM_C" "$voxmm")

    fslmaths "$OUT_SUB/lesion_in_CONNfunc.nii.gz" -dilM -kernel 3D -dilM "$TX/dilA"
    [[ -f "$TX/dilA.nii.gz" ]] || { log "  ❌ failed to create $TX/dilA.nii.gz"; echo "$CONN_ID,$mode,fail_dilate,$ANAT_REF,$FUNC_REF,$SRC_ID,$T1_SRC,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"; continue; }
    for i in $(seq 2 $nA); do fslmaths "$TX/dilA" -dilM "$TX/dilA"; done
    fslmaths "$TX/dilA" -sub "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_0to2mm.nii.gz"
    fslmaths "$OUT_SUB/ring_0to2mm.nii.gz" -bin "$OUT_SUB/ring_0to2mm.nii.gz"

    if have_imcp; then imcp "$TX/dilA" "$TX/dilB"; else cp "$TX/dilA.nii.gz" "$TX/dilB.nii.gz"; fi
    for i in $(seq $((nA+1)) $nB); do fslmaths "$TX/dilB" -dilM "$TX/dilB"; done
    fslmaths "$TX/dilB" -sub "$TX/dilA" "$OUT_SUB/ring_2to4mm.nii.gz"
    fslmaths "$OUT_SUB/ring_2to4mm.nii.gz" -bin "$OUT_SUB/ring_2to4mm.nii.gz"

    if have_imcp; then imcp "$TX/dilB" "$TX/dilC"; else cp "$TX/dilB.nii.gz" "$TX/dilC.nii.gz"; fi
    for i in $(seq $((nB+1)) $nC); do fslmaths "$TX/dilC" -dilM "$TX/dilC"; done
    fslmaths "$TX/dilC" -sub "$TX/dilB" "$OUT_SUB/ring_4to8mm.nii.gz"
    fslmaths "$OUT_SUB/ring_4to8mm.nii.gz" -bin "$OUT_SUB/ring_4to8mm.nii.gz"

    rm -f "$TX"/dilA* "$TX"/dilB* "$TX"/dilC* 2>/dev/null || true
    status="done_patient"

  else
    # CONTROL (or missing lesion/T1) → write blanks (in FUNC grid) once
    if already_done_func "$OUT_SUB"; then
      mode="normal"
      log "  → CONTROL/MISSING → create blank ROIs"
      fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/lesion_in_CONNfunc.nii.gz"
      cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_0to2mm.nii.gz"
      cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_2to4mm.nii.gz"
      cp "$OUT_SUB/lesion_in_CONNfunc.nii.gz" "$OUT_SUB/ring_4to8mm.nii.gz"
      fslmaths "$FUNC_REF3D" -mul 0 "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz"
      for thr in 0.3 0.5 0.7; do
        cp "$OUT_SUB/disconnectome_in_CONNfunc_prob.nii.gz" \
           "$OUT_SUB/disconnectome_thr${thr}_in_CONNfunc.nii.gz"
      done
      status="done_control"
    else
      log "  ✓ CONTROL already has outputs — skip"
      status="skip_already_done_control"
    fi
  fi

  echo "$CONN_ID,$mode,$status,$ANAT_REF,$FUNC_REF,$SRC_ID,$T1_SRC,$LES_SRC,$DISC_SRC,$OUT_SUB" >> "$MANIFEST"
done

log "All subjects processed. Outputs → $OUT_ROOT"
log "Manifest → $MANIFEST"
