# fMRIPrep anatomical outputs
SOURCE_DIR="/home/sskgroup/Documents/stroke_test/output23.0.223/"
DEST_DIR="/home/sskgroup/Documents/stroke_test/output24.0.1/"

# Original sourcedata location
SOURCE_SOURCEDATA="/home/sskgroup/Documents/stroke_test/output23.0.223/sourcedata"
DEST_SOURCEDATA="/home/sskgroup/Documents/stroke_test/output24.0.1/sourcedata"

# Create destination folders
mkdir -p "$DEST_DIR"

# Copy anat folders
for subj in "$SOURCE_DIR"/sub-*; do
  subj_id=$(basename "$subj")
  src_anat="$subj/anat"
  dest_anat="$DEST_DIR/$subj_id/anat"

  if [ -d "$src_anat" ]; then
    mkdir -p "$DEST_DIR/$subj_id"
    cp -r "$src_anat" "$dest_anat"
    echo "✅ Copied anat for $subj_id"
  else
    echo "⚠️ Missing anat for $subj_id"
  fi
done

# Copy sourcedata folder (single copy)
if [ -d "$SOURCE_SOURCEDATA" ]; then
  cp -r "$SOURCE_SOURCEDATA" "$DEST_SOURCEDATA"
  echo "✅ Copied sourcedata folder"
else
  echo "⚠️ sourcedata folder not found at $SOURCE_SOURCEDATA"
fi