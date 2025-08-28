#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────
# lesion-2mm-to-1mm-v1.sh
# Upsample MNI 2mm lesion masks → MNI 1mm using ANTs (NN) + FSL.
# Writes outputs to a separate OUTROOT, mirroring the input tree.
#
# REQUIREMENTS: antsApplyTransforms, fslreorient2std, fslcpgeom, awk, sed, realpath
# ─────────────────────────────────────────────────────────────

# === HARD-CODE THESE PATHS ===
INROOT="/media/sskgroup/fast/Lesion/Secondary/lesion_MNI"          # contains *_lesion_MNI2mm.nii.gz
OUTROOT="/media/sskgroup/fast/Lesion/Secondary/Lesion_MNI_1mm"     # where *_lesion_MNI1mm.nii.gz will be written
REF1="/media/sskgroup/fast/Lesion/BCBToolKit/Tools/extraFiles/MNI152.nii.gz"   # MNI152 1mm reference
THREADS=16  # ANTs threads (optional)
# =============================

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${THREADS}
export OMP_NUM_THREADS=${THREADS}

die(){ echo "[ERR] $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }

need antsApplyTransforms
need fslreorient2std
need fslcpgeom
need awk
need sed
need realpath

[[ -d "$INROOT"  ]] || die "INROOT not found: $INROOT"
mkdir -p "$OUTROOT"
[[ -f "$REF1"    ]] || die "REF1 not found: $REF1"

# Show reference grid (helps catch wrong ref)
echo "[i] Using 1mm reference: $REF1"
fslhd "$REF1" | egrep 'dim[123]|pixdim[123]'

# Walk all 2mm lesions and upsample to 1mm into OUTROOT, mirroring structure
echo "[*] Scanning for *_lesion_MNI2mm.nii.gz under $INROOT ..."
mapfile -t FILES < <(find "$INROOT" -type f -name "*.nii.gz" | sort)
[[ ${#FILES[@]} -gt 0 ]] || die "No lesion files matching *_lesion_MNI2mm.nii.gz found in $INROOT"

for L2 in "${FILES[@]}"; do
  # relative path of L2 w.r.t. INROOT
  rel="${L2#$INROOT/}"
  # output path under OUTROOT, rename _MNI2mm → _MNI1mm
  out_rel="${rel/_MNI2mm/_MNI1mm}"
  L1="${OUTROOT}/${out_rel}"

  # ensure output directory exists
  mkdir -p "$(dirname "$L1")"

  echo "[*] upsampling $(basename "$L2") → $(basename "$L1")"

  # 1) resample 2mm → 1mm (NearestNeighbor to preserve binary mask)
  antsApplyTransforms -d 3 -i "$L2" -r "$REF1" -t identity -n NearestNeighbor -o "$L1"

  # 2) standardize orientation & copy geometry to avoid pipeline warnings
  fslreorient2std "$L1" "${L1%.nii.gz}_reor.nii.gz"
  mv "${L1%.nii.gz}_reor.nii.gz" "$L1"
  fslcpgeom "$REF1" "$L1" || true

  # 3) quick header check
  dims_les=$(fslhd "$L1" | awk '/^dim1|^dim2|^dim3/ {printf "%s ",$2}')
  pixs_les=$(fslhd "$L1" | awk '/^pixdim1|^pixdim2|^pixdim3/ {printf "%s ",$2}')
  echo "    dims: $dims_les   pixdim: $pixs_les"
done

echo "[OK] Converted ${#FILES[@]} lesion(s) to 1mm → $OUTROOT"
