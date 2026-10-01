#!/usr/bin/env python3
"""Phase 11 — build the ML feature matrix (PSI) from rMATS outputs.

For each age (3d, 6d, 9d), rMATS tested all annotated events and wrote per-sample
junction counts (IJC/SJC) and PSI (IncLevel) for the 6 pseudo-bulk samples of
that age (3 male + 3 female replicates).

This script:
  1. Parses the 5 event-type JC files per age.
  2. Unions events across ages (same GTF => same event universe).
  3. Builds an 18 x N_events PSI matrix (rows = pseudo-bulk samples).
  4. Applies QC filters: complete PSI across all 18 samples AND >=10 total
     supporting junction reads (IJC+SJC summed over samples).
  5. Writes results/psi/psi_matrix_all.csv, results/psi/psi_matrix_qc.csv,
     results/ml_feature_inventory.csv, and a feature dictionary.

Event ID format: <age>_<type>_<rMATS_ID> to keep age-specific definitions distinct
(the same genomic event tested at two ages is treated as one feature; PSI comes
from the age-appropriate rMATS run).
"""
import glob
import os
import re

import pandas as pd

BASE = os.path.expanduser("~/workspace/ML_Ozgun")
AGES = ["3", "6", "9"]
TYPES = ["SE", "A3SS", "A5SS", "MXE", "RI"]

# sample order per age: b1 = male reps, b2 = female reps (see run_rmats.sh)
def samples_for_age(a):
    return ([f"w1118_{a}d_r{r}_male" for r in (1, 2, 3)] +
            [f"w1118_{a}d_r{r}_female" for r in (1, 2, 3)])

def parse_jc(path, age, etype):
    df = pd.read_csv(path, sep="\t")
    # per-sample counts: IJC_SAMPLE_1 / SJC_SAMPLE_1 (b1=male), _2 (b2=female)
    recs = []
    for _, row in df.iterrows():
        ijc1 = str(row["IJC_SAMPLE_1"]).split(",")
        sjc1 = str(row["SJC_SAMPLE_1"]).split(",")
        ijc2 = str(row["IJC_SAMPLE_2"]).split(",")
        sjc2 = str(row["SJC_SAMPLE_2"]).split(",")
        inc1 = str(row["IncLevel1"]).split(",")
        inc2 = str(row["IncLevel2"]).split(",")
        counts = ijc1 + ijc2
        sjcs = sjc1 + sjc2
        incs = inc1 + inc2
        samp = samples_for_age(age)
        assert len(counts) == 6 and len(incs) == 6, f"unexpected replicate count in {path}"
        d = {"event_id": f"{age}d_{etype}_{row['ID']}",
             "age": age, "event_type": etype,
             "gene_id": row["GeneID"], "gene_symbol": row.get("geneSymbol", ""),
             "chr": row["chr"], "strand": row["strand"],
             "pvalue": row["PValue"], "fdr": row["FDR"],
             "delta_psi": row["IncLevelDifference"]}
        for s, ijc, sjc, psi in zip(samp, counts, sjcs, incs):
            ijc_i, sjc_i = int(ijc), int(sjc)
            try:
                psi_f = float(psi)
            except ValueError:
                psi_f = float("nan")
            d[f"psi_{s}"] = psi_f
            d[f"reads_{s}"] = ijc_i + sjc_i
        recs.append(d)
    return pd.DataFrame(recs)

all_events = []
for age in AGES:
    for et in TYPES:
        pat = os.path.join(BASE, "results/rmats", f"{age}d", f"{et}.MATS.JC.txt")
        if not os.path.exists(pat):
            print(f"MISSING {pat}")
            continue
        df = parse_jc(pat, age, et)
        print(f"{age}d {et}: {len(df)} tested events")
        all_events.append(df)

events = pd.concat(all_events, ignore_index=True)
psi_cols = [c for c in events.columns if c.startswith("psi_")]
read_cols = [c for c in events.columns if c.startswith("reads_")]

os.makedirs(os.path.join(BASE, "results/psi"), exist_ok=True)
events.to_csv(os.path.join(BASE, "results/psi/psi_matrix_all.csv"), index=False)
print(f"union: {len(events)} events x {len(psi_cols)} samples")

# QC: complete PSI + >=10 total supporting reads
events["n_missing"] = events[psi_cols].isna().sum(axis=1)
events["total_reads"] = events[read_cols].sum(axis=1)
qc = events[(events["n_missing"] == 0) & (events["total_reads"] >= 10)].copy()
print(f"QC-passed (complete PSI, >=10 reads): {len(qc)}")
print(qc["event_type"].value_counts().to_string())

# feature dictionary
feat = qc[["event_id", "age", "event_type", "gene_id", "gene_symbol",
           "chr", "strand", "pvalue", "fdr", "delta_psi", "total_reads"]].copy()
feat.to_csv(os.path.join(BASE, "results/ml_feature_inventory.csv"), index=False)

# ML matrix: rows=samples, cols=events (QC-passed only)
samples = [c.replace("psi_", "") for c in psi_cols]
ml = pd.DataFrame(index=samples)
for _, r in qc.iterrows():
    ml[r["event_id"]] = [r[c] for c in psi_cols]
ml.index.name = "sample"
ml.to_csv(os.path.join(BASE, "results/psi/psi_matrix_qc.csv"))
print(f"ML matrix: {ml.shape[0]} samples x {ml.shape[1]} features -> results/psi/psi_matrix_qc.csv")
