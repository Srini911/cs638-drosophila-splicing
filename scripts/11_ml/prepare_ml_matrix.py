#!/usr/bin/env python3
"""
Convert psi_matrix.tsv (events x samples) to psi_matrix_qc.csv (samples x events)
for ML. Also filters events with too many NAs.
Usage: prepare_ml_matrix.py
"""
import csv
from pathlib import Path

def main():
    base = Path.home() / "workspace" / "ML_Ozgun"
    infile = base / "data" / "intermediate" / "psi_matrix.tsv"
    outdir = base / "results" / "psi"
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "psi_matrix_qc.csv"

    # Read TSV
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        samples = [c for c in reader.fieldnames if c != "event_id"]
        events = []
        data = {}  # event_id -> {sample: psi}
        for row in reader:
            eid = row["event_id"]
            events.append(eid)
            data[eid] = {s: row[s] for s in samples}

    print(f"Input: {len(events)} events x {len(samples)} samples", flush=True)

    # Filter: keep events with PSI in ALL samples, and variance > 0
    # (Strict: no missing values, to avoid imputation complications with n=4)
    import numpy as np
    kept = []
    for eid in events:
        vals = []
        has_na = False
        for s in samples:
            v = data[eid][s]
            if v == "NA":
                has_na = True
                break
            vals.append(float(v))
        # Require no NAs and non-zero variance
        if not has_na and np.std(vals) > 1e-6:
            kept.append(eid)

    print(f"Kept {len(kept)} events (variance > 0, >=50% samples)", flush=True)

    # Write CSV (samples x events)
    with open(outfile, "w", newline="") as o:
        w = csv.writer(o)
        w.writerow(["sample"] + kept)
        for s in samples:
            row = [s]
            for eid in kept:
                v = data[eid][s]
                row.append(v if v != "NA" else "")
            w.writerow(row)

    print(f"Wrote {outfile}: {len(samples)} samples x {len(kept)} features", flush=True)

if __name__ == "__main__":
    main()
