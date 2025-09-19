#!/usr/bin/env python3
import os, glob, pandas as pd

# ----- HARD-CODED PATHS -----
FMRIPREP = "/media/sskgroup/noise/primary/lag_confounds/batch_7"
LAG_DIR  = "/media/sskgroup/project/Thesis/fmristroke/primary/confounds"
# ----------------------------

lag_files = sorted(glob.glob(os.path.join(LAG_DIR, "sub-*_*_desc-confounds_timeseries.tsv")))
n_merged = n_skip = 0

for lag_path in lag_files:
    base = os.path.basename(lag_path)  # already corrected name
    sub  = base.split("_")[0]          # e.g. sub-052
    func = os.path.join(FMRIPREP, sub, "func")
    dst  = os.path.join(func, base)

    if not os.path.exists(dst):
        print(f"[skip] missing fMRIPrep confounds: {dst}")
        n_skip += 1
        continue

    # read original and lag
    a = pd.read_table(dst)
    b = pd.read_table(lag_path).select_dtypes('number')

    if len(a) != len(b):
        print(f"[skip] row mismatch {base}: fmriprep={len(a)} lag={len(b)}")
        n_skip += 1
        continue

    # backup once
    bak = dst.replace("_desc-confounds_timeseries.tsv",
                      "_desc-confounds_timeseries.ORIG.tsv")
    if not os.path.exists(bak):
        a.to_csv(bak, sep='\t', index=False)

    # prefix lag columns
    b = b.add_prefix("lag_")
    merged = pd.concat([a, b], axis=1)
    merged.to_csv(dst, sep='\t', index=False)

    print(f"[merged] {dst}")
    n_merged += 1

print(f"done: merged={n_merged}, skipped={n_skip}")
