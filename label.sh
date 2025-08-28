#!/bin/bash

cd /home/sskgroup/Documents/atlases

declare -A rename_map=(
  ["Cerebellum_MNIflirt.txt"]="Cerebellum-MNIflirt-maxprob-thr25-2mm.txt"
  ["Cerebellum_MNIfnirt.txt"]="Cerebellum-MNIfnirt-maxprob-thr25-2mm.txt"
  ["HarvardOxford-Cortical.txt"]="HarvardOxford-cort-maxprob-thr25-2mm.txt"
  ["HarvardOxford-Subcortical.txt"]="HarvardOxford-sub-maxprob-thr25-2mm.txt"
  ["JHU-labels.txt"]="JHU-ICBM-labels-2mm.txt"
  ["JHU-tracts.txt"]="JHU-ICBM-tracts-maxprob-thr25-2mm.txt"
  ["Juelich.txt"]="Juelich-maxprob-thr25-2mm.txt"
  ["MarsParietalParcellation.txt"]="Parietal_thr25_2mm.txt"
  ["MarsTPJParcellation.txt"]="TPJ_thr25_2mm.txt"
  ["MNI.txt"]="MNI-maxprob-thr25-2mm.txt"
  ["NeubertVentralFrontalParcellation.txt"]="VentralFrontal_thr25_2mm.txt"
  ["SalletDorsalFrontalParcellation.txt"]="DorsalFrontal_thr25_2mm.txt"
  ["Striatum-Connectivity-3sub.txt"]="striatum-con-label-thr25-3sub-2mm.txt"
  ["Striatum-Connectivity-7sub.txt"]="striatum-con-label-thr25-7sub-2mm.txt"
  ["Striatum-Structural.txt"]="striatum-structural-2mm.txt"
  ["Talairach.txt"]="Talairach-labels-2mm.txt"
  ["Thalamus.txt"]="Thalamus-maxprob-thr25-2mm.txt"
  ["SMATT.txt"]="SMATT-labels-1mm.txt"
  ["XTRACT.txt"]="xtract-tract-atlases-maxprob5-1mm.txt"
)

for old in "${!rename_map[@]}"; do
  if [[ -f "$old" ]]; then
    new="${rename_map[$old]}"
    echo "🔁 Renaming $old → $new"
    mv "$old" "$new"
  else
    echo "⚠️  Missing: $old"
  fi
done
