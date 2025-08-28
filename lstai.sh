#!/bin/bash

# ------------------- Hardcoded Paths -------------------

BIDS_ROOT="/home/sskgroup/Documents/Lesion_BIDS/BIDS"
OUTPUT_ROOT="/home/sskgroup/Documents/Lesion_BIDS/LST_AI/Primary"
DATASET="Primary"
TEMP_ROOT="$OUTPUT_ROOT/lst_temp"

# Create required folders
mkdir -p "$OUTPUT_ROOT"
mkdir -p "$TEMP_ROOT"

# ------------------- Subject Loop -------------------

for SUBJ_DIR in "$BIDS_ROOT/$DATASET"/sub-*/anat; do
  SUB=$(basename "$(dirname "$SUBJ_DIR")")
  echo "🔁 Processing $SUB..."

  INPUT_DIR="$SUBJ_DIR"
  OUTPUT_DIR="$OUTPUT_ROOT/$SUB"
  mkdir -p "$OUTPUT_DIR"

  T1_FILE="$INPUT_DIR/${SUB}_T1w.nii.gz"
  FLAIR_FILE="$INPUT_DIR/${SUB}_FLAIR_resampled_n4.nii.gz"

  # Fallback to non-N4
  if [[ ! -f "$FLAIR_FILE" ]]; then
    FLAIR_FILE="$INPUT_DIR/${SUB}_FLAIR_resampled.nii.gz"
  fi

  if [[ ! -f "$T1_FILE" ]]; then
    echo "❌ T1w not found for $SUB"
    continue
  fi

  if [[ ! -f "$FLAIR_FILE" ]]; then
    echo "❌ FLAIR (resampled) not found for $SUB"
    continue
  fi

  echo "🧠 Running LST-AI for $SUB..."
  sudo docker run --rm --gpus all \
    -v "$INPUT_DIR":/custom_apps/lst_input \
    -v "$OUTPUT_DIR":/custom_apps/lst_output \
    -v "$TEMP_ROOT":/custom_apps/lst_temp \
    jqmcginnis/lst-ai:v1.2.0 \
    --t1 /custom_apps/lst_input/$(basename "$T1_FILE") \
    --flair /custom_apps/lst_input/$(basename "$FLAIR_FILE") \
    --output /custom_apps/lst_output \
    --temp /custom_apps/lst_temp

  # ---------------- Rename lesion output ----------------
  RAW_MASK=$(find "$OUTPUT_DIR" -iname "*lesion*.nii*" | head -n 1)
  if [[ -f "$RAW_MASK" ]]; then
    NEW_NAME="${SUB}_label-lesion_roi.nii.gz"
    cp "$RAW_MASK" "$SUBJ_DIR/$NEW_NAME"
    echo "✅ Saved lesion mask as: $NEW_NAME"
  else
    echo "⚠️ Lesion mask not found in output for $SUB"
  fi

  echo "✅ Finished $SUB"
  echo "-------------------------------"
done
