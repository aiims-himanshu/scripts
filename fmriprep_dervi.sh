#!/bin/bash

SOURCE_DIR="/home/sskgroup/Documents/stroke_test/output23.0.223/fmriprep"
DEST_DIR="/home/sskgroup/Documents/stroke_test/output24.0.1/fmriprep"

mkdir -p "$DEST_DIR"

for subj in "$SOURCE_DIR"/sub-*; do
  subj_id=$(basename "$subj")
  src_anat="$subj/anat"
  dest_anat="$DEST_DIR/$subj_id/anat"

  if [ -d "$src_anat" ]; then
    mkdir -p "$DEST_DIR/$subj_id"
    cp -r "$src_anat" "$dest_anat"
    echo "✅ Copied $subj_id"
  else
    echo "⚠️ Missing anat for $subj_id"
  fi
done
