#!/usr/bin/env bash
set -euo pipefail

# ╭──────────────────────────────────────────────────────────────╮
# │ lesion_clean.sh                                              │
# │ Review and clean tiny disconnected lesion islands           │
# │ from binary lesion masks                                     │
# ╰──────────────────────────────────────────────────────────────╯
#
# Two-step workflow:
#   1) scan: create a review CSV of small connected components
#   2) apply: use the edited review CSV to remove only flagged components
#
# Default lesion layout:
#   /media/sskgroup/Thesis/Primary/roi/Lesion/sub-XXX_lesion.nii.gz
#
# Example:
#   bash lesion_clean.sh --scan  --lesion-root /media/sskgroup/Thesis/Primary/roi/Lesion --out-root /tmp/lesion_qc
#   edit the generated review CSV
#   bash lesion_clean.sh --apply --lesion-root /media/sskgroup/Thesis/Primary/roi/Lesion --out-root /tmp/lesion_qc --review-csv /tmp/lesion_qc/review_....csv
#
# Requirements:
#   python3 with numpy, scipy, nibabel

LESION_ROOT="/media/sskgroup/Thesis/Primary/roi/Lesion"
OUT_ROOT="/media/sskgroup/Thesis/Primary/roi/lesion_clean"
REVIEW_CSV="/media/sskgroup/Thesis/Primary/roi/lesion_clean/lesion_review_20260725_113858.csv"
MODE="apply" #scan or apply
MAX_COMPONENT_VOX=4
CONNECTIVITY=26
THREADS="${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS:-1}"

die(){ echo "ERROR: $*" >&2; exit 1; }
log(){ echo "[`date +%F' '%T`] $*"; }

have(){ command -v "$1" >/dev/null 2>&1; }

usage(){
  sed -n '1,220p' "$0"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --lesion-root)      LESION_ROOT="$2"; shift 2;;
    --out-root)         OUT_ROOT="$2"; shift 2;;
    --review-csv)       REVIEW_CSV="$2"; shift 2;;
    --max-component-vox) MAX_COMPONENT_VOX="$2"; shift 2;;
    --connectivity)     CONNECTIVITY="$2"; shift 2;;
    --threads)          THREADS="$2"; shift 2;;
    --scan)             MODE="scan"; shift;;
    --apply)            MODE="apply"; shift;;
    --both)             MODE="both"; shift;;
    -h|--help)          usage; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -d "$LESION_ROOT" ]] || die "--lesion-root not found"
mkdir -p "$OUT_ROOT"

python3 - <<'PY' >/dev/null 2>&1 || die "python3 needs numpy, scipy, and nibabel"
import numpy, scipy, nibabel
PY

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"

scan_py() {
python3 - "$LESION_ROOT" "$OUT_ROOT" "$MAX_COMPONENT_VOX" "$CONNECTIVITY" <<'PY'
import csv
import os
import sys
from datetime import datetime

import nibabel as nib
import numpy as np
from scipy import ndimage

lesion_root = sys.argv[1]
out_root = sys.argv[2]
max_component_vox = int(sys.argv[3])
connectivity = int(sys.argv[4])

if connectivity not in (6, 18, 26):
    raise SystemExit("connectivity must be 6, 18, or 26")

if connectivity == 6:
    structure = ndimage.generate_binary_structure(3, 1)
elif connectivity == 18:
    structure = ndimage.generate_binary_structure(3, 2)
else:
    structure = np.ones((3, 3, 3), dtype=bool)

ts = datetime.now().strftime("%Y%m%d_%H%M%S")
review_csv = os.path.join(out_root, f"lesion_review_{ts}.csv")
summary_csv = os.path.join(out_root, f"lesion_summary_{ts}.csv")

files = []
for name in sorted(os.listdir(lesion_root)):
    if name.startswith("sub-") and (name.endswith(".nii") or name.endswith(".nii.gz")):
        files.append(os.path.join(lesion_root, name))

if not files:
    raise SystemExit(f"No lesion files found in {lesion_root}")

with open(review_csv, "w", newline="") as review_fh, open(summary_csv, "w", newline="") as summary_fh:
    review = csv.writer(review_fh)
    summary = csv.writer(summary_fh)
    review.writerow([
        "subject","file","component_id","component_voxels",
        "remove_flag","bbox_ijk","centroid_ijk","note"
    ])
    summary.writerow([
        "subject","file","total_voxels","num_components",
        "num_small_components","small_component_voxels_total",
        "largest_component_voxels","flagged_smallest_voxels"
    ])

    for path in files:
        base = os.path.basename(path)
        subject = base.split("_")[0]

        img = nib.load(path)
        data = img.get_fdata()
        mask = data > 0

        if not np.any(mask):
            summary.writerow([subject, path, 0, 0, 0, 0, 0, "no_foreground"])
            continue

        lab, nlab = ndimage.label(mask, structure=structure)
        counts = np.bincount(lab.ravel())
        counts[0] = 0

        total_vox = int(mask.sum())
        largest = int(counts.max()) if counts.size > 1 else 0
        small_components = 0
        small_total = 0
        flagged_smallest = 0

        for comp_id in range(1, nlab + 1):
            vox = int(counts[comp_id])
            if vox <= 0:
                continue

            if vox <= max_component_vox:
                small_components += 1
                small_total += vox
                flagged_smallest = vox if flagged_smallest == 0 else min(flagged_smallest, vox)
                coords = np.argwhere(lab == comp_id)
                mins = coords.min(axis=0).tolist()
                maxs = coords.max(axis=0).tolist()
                centroid = coords.mean(axis=0)
                review.writerow([
                    subject,
                    path,
                    comp_id,
                    vox,
                    1,
                    f"{mins}-{maxs}",
                    f"{centroid[0]:.2f},{centroid[1]:.2f},{centroid[2]:.2f}",
                    "small_component_candidate"
                ])

        summary.writerow([
            subject,
            path,
            total_vox,
            int(nlab),
            small_components,
            small_total,
            largest,
            flagged_smallest
        ])

print(review_csv)
print(summary_csv)
PY
}

apply_py() {
python3 - "$LESION_ROOT" "$OUT_ROOT" "$REVIEW_CSV" "$CONNECTIVITY" <<'PY'
import csv
import os
import sys
from collections import defaultdict

import nibabel as nib
import numpy as np
from scipy import ndimage

lesion_root = sys.argv[1]
out_root = sys.argv[2]
review_csv = sys.argv[3]
connectivity = int(sys.argv[4])

if not os.path.isfile(review_csv):
    raise SystemExit(f"Review CSV not found: {review_csv}")

def parse_bool(val, default=True):
    if val is None or val == "":
        return default
    val = str(val).strip().lower()
    if val in {"1", "true", "t", "yes", "y", "remove"}:
        return True
    if val in {"0", "false", "f", "no", "n", "keep"}:
        return False
    return default

remove_map = defaultdict(set)
with open(review_csv, newline="") as fh:
    reader = csv.DictReader(fh)
    if not reader.fieldnames:
        raise SystemExit("Review CSV has no header")
    for row in reader:
        file_path = row.get("file") or row.get("path")
        comp_id = row.get("component_id") or row.get("component")
        remove_flag = parse_bool(row.get("remove_flag"), default=True)
        if file_path and comp_id and remove_flag:
            remove_map[file_path].add(int(comp_id))

if not os.path.isdir(out_root):
    os.makedirs(out_root, exist_ok=True)

subject_dirs = sorted(
    f for f in os.listdir(lesion_root)
    if f.startswith("sub-") and (f.endswith(".nii") or f.endswith(".nii.gz"))
)

if not subject_dirs:
    raise SystemExit(f"No lesion files found in {lesion_root}")

for name in subject_dirs:
    src = os.path.join(lesion_root, name)
    dst = os.path.join(out_root, name)

    img = nib.load(src)
    data = img.get_fdata()
    mask = data > 0

    if connectivity == 6:
        structure = ndimage.generate_binary_structure(3, 1)
    elif connectivity == 18:
        structure = ndimage.generate_binary_structure(3, 2)
    else:
        structure = np.ones((3, 3, 3), dtype=bool)

    if src in remove_map and np.any(mask):
        lab, nlab = ndimage.label(mask, structure=structure)
        remove_ids = remove_map[src]
        if remove_ids:
            cleaned = mask.copy()
            for comp_id in remove_ids:
                cleaned[lab == comp_id] = False
            data_out = cleaned.astype(np.float32)
        else:
            data_out = mask.astype(np.float32)
    else:
        data_out = mask.astype(np.float32)

    out_img = nib.Nifti1Image(data_out, img.affine, img.header)
    out_img.set_data_dtype(np.float32)
    nib.save(out_img, dst)

print(out_root)
PY
}

if [[ "$MODE" == "scan" ]]; then
  scan_py
  exit 0
fi

if [[ "$MODE" == "apply" ]]; then
  [[ -n "$REVIEW_CSV" ]] || die "--review-csv is required in apply mode"
  apply_py
  exit 0
fi

if [[ "$MODE" == "both" ]]; then
  mapfile -t _SCAN_OUT < <(scan_py)
  REVIEW_CSV="${_SCAN_OUT[0]}"
  apply_py
  exit 0
fi

die "Unknown mode: $MODE"
