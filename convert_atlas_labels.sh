#!/bin/bash

# === Setup ===
SRC_DIR="/home/sskgroup/Downloads/lesion_correction-main/space/"
DEST_DIR="/home/sskgroup/Documents/atlases"
mkdir -p "$DEST_DIR"

# === Atlas & Label File List ===

declare -A atlases=(
  ["Cerebellum"]="Cerebellum-MNIflirt-maxprob-thr25-2mm.nii.gz"
  ["HarvardOxford"]="HarvardOxford-cort-maxprob-thr25-2mm.nii.gz HarvardOxford-sub-maxprob-thr25-2mm.nii.gz"
  ["JHU"]="JHU-ICBM-labels-2mm.nii.gz JHU-ICBM-tracts-maxprob-thr25-2mm.nii.gz"
  ["Juelich"]="Juelich-maxprob-thr25-2mm.nii.gz"
  ["MNI"]="MNI-maxprob-thr25-2mm.nii.gz"
  ["MarsParietalParcellation"]="Parietal_thr25_2mm.nii.gz"
  ["MarsTPJParcellation"]="TPJ_thr25_2mm.nii.gz"
  ["NeubertVentralFrontalParcellation"]="VentralFrontal_thr25_2mm.nii.gz"
  ["SalletDorsalFrontalParcellation"]="DorsalFrontal_thr25_2mm.nii.gz"
  ["Thalamus"]="Thalamus-maxprob-thr25-2mm.nii.gz"
  ["Striatum"]="striatum-con-label-thr25-7sub-2mm.nii.gz"
  ["Talairach"]="Talairach-labels-2mm.nii.gz"
  ["XTRACT"]="xtract-tract-atlases-maxprob5-1mm.nii.gz"
  ["SMATT"]="SMATT-labels-1mm.nii.gz"
)

# === Process Each Atlas ===

for folder in "${!atlases[@]}"; do
  for atlas_file in ${atlases[$folder]}; do
    echo "📦 Copying $atlas_file from $folder"
    cp "$SRC_DIR/$folder/$atlas_file" "$DEST_DIR/" 2>/dev/null || echo "⚠️ File not found: $SRC_DIR/$folder/$atlas_file"
    
    xmlfile="$SRC_DIR/$folder.xml"
    if [[ ! -f "$xmlfile" ]]; then
      xmlfile="$SRC_DIR/$folder.xml"  # fallback: top-level
    fi

    # Extract base name for label
    base_name="${folder%%.*}"

    if [[ -f "$SRC_DIR/$base_name.xml" ]]; then
      echo "🔁 Converting $base_name.xml to $DEST_DIR/${base_name}.txt"
      xmlstarlet sel -t -m "//label" -v "@index" -o " " -v "." -n "$SRC_DIR/$base_name.xml" > "$DEST_DIR/${base_name}.txt"
    else
      echo "❌ No matching .xml found for $folder"
    fi

    echo ""
  done
done

echo "✅ Conversion complete."
