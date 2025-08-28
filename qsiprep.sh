#!/bin/bash

# ------------ CONFIGURATION ------------
QSIPREP_IMAGE="pennlinc/qsiprep:1.0.1"
BIDS_DIR="/home/sskgroup/Documents/qsiprep/Primary"
OUTPUT_DIR="/home/sskgroup/Documents/qsiprep/qsiprep"
WORK_DIR="/home/sskgroup/Documents/qsiprep/tmpo"
FS_LICENSE="/home/sskgroup/Documents/license/license.txt"
LOGFILE="qsiprep_batch_log_$(date +%Y%m%d_%H%M%S).txt"

# Number of parallel subjects
BATCH_SIZE=5

# ------------ INIT ------------
mkdir -p "$OUTPUT_DIR" "$WORK_DIR"
echo "🔁 Starting QSIPrep batch run" > "$LOGFILE"

# ------------ FUNCTION ------------
run_subject() {
    subj_id=$1
    subj=sub-$(printf "%02d" "$subj_id")

    anat_file="${BIDS_DIR}/${subj}/anat/${subj}_T1w.nii.gz"
    dwi_file="${BIDS_DIR}/${subj}/dwi/${subj}_dwi.nii.gz"
    bval_file="${BIDS_DIR}/${subj}/dwi/${subj}_dwi.bval"
    bvec_file="${BIDS_DIR}/${subj}/dwi/${subj}_dwi.bvec"

    if [[ ! -f "$anat_file" || ! -f "$dwi_file" || ! -f "$bval_file" || ! -f "$bvec_file" ]]; then
        echo "⚠️ Skipping $subj — missing files" | tee -a "$LOGFILE"
        return
    fi

    echo "🚀 Running $subj..." | tee -a "$LOGFILE"

    time sudo docker run -ti --rm \
        -v "${BIDS_DIR}:/data:ro" \
        -v "${OUTPUT_DIR}:/out" \
        -v "${WORK_DIR}:/w" \
        -v "$(dirname "$FS_LICENSE"):/freesurfer" \
        "$QSIPREP_IMAGE" \
        /data /out \
        participant \
        -w /w \
        --omp-nthreads 3 \
        --fs-license-file /freesurfer/$(basename "$FS_LICENSE") \
        --anat-modality T1w \
        --output-resolution 2.5 \
        --anatomical-template MNI152NLin2009cAsym \
        --dwi-denoise-window auto \
        --denoise-method dwidenoise \
        --unringing-method mrdegibbs \
        --b1-biascorrect-stage none \
        --b0-to-t1w-transform Affine \
        --hmc-transform Affine \
        --hmc-model eddy \
        --shoreline-iters 5 \
        --b0-motion-corr-to iterative \
        --write-graph \
        --participant-label "$subj" \
        --notrack \
        >> "$LOGFILE" 2>&1

    if [[ $? -eq 0 ]]; then
        echo "✅ Finished $subj successfully" | tee -a "$LOGFILE"
    else
        echo "❌ Error in $subj — check logs" | tee -a "$LOGFILE"
    fi
}

# ------------ MAIN LOOP ------------
subject_ids=($(seq -w 1 51))  # sub-01 to sub-51

for ((i=0; i<${#subject_ids[@]}; i+=BATCH_SIZE)); do
    echo "🔁 Batch starting: ${subject_ids[@]:i:BATCH_SIZE}" | tee -a "$LOGFILE"

    for j in "${subject_ids[@]:i:BATCH_SIZE}"; do
        run_subject "$j" &
    done

    wait
    echo "✅ Batch completed." | tee -a "$LOGFILE"
done

echo "🏁 All batches processed. Log: $LOGFILE"
