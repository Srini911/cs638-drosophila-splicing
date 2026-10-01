#!/usr/bin/env python3
"""Phase 4 — reproduce and validate roX2-based sex classification.

roX2 (FBgn0019660) is a male-specific lncRNA involved in X-chromosome dosage
compensation in Drosophila. It is effectively not expressed in females, so any
cell with roX2 UMIs is classified male.

This script:
  1. Extracts the roX2 count vector from the GEO 57k count matrix (streaming).
  2. Classifies each cell: roX2 UMIs > 0 -> male, else female.
  3. Compares against the GEO metadata `sex` column (treated as reference, NOT
     blindly inherited) and reports agreement / disagreement.
  4. Writes metadata/cell_metadata_with_roX2_expression.csv,
     metadata/female_barcodes.tsv, metadata/male_barcodes.tsv.

Usage: python3 validate_rox2.py
"""
import gzip
import os
import subprocess
import pandas as pd

BASE = os.path.expanduser("~/workspace/ML_Ozgun")
GEO = os.path.join(BASE, "data/raw/geo")
META = os.path.join(BASE, "metadata")
os.makedirs(META, exist_ok=True)

# 1. locate roX2 row (1-indexed) in genes.tsv
rox2_idx = None
with open(os.path.join(GEO, "matrix57k/genes.tsv")) as fh:
    for i, line in enumerate(fh, start=1):
        if line.startswith("FBgn0019660\t"):
            rox2_idx = i
            break
assert rox2_idx is not None, "roX2 (FBgn0019660) not found in genes.tsv"
print(f"roX2 matrix row: {rox2_idx}")

# 2. stream the MTX, keep only roX2 entries (col -> UMI count)
mtx = os.path.join(GEO, "matrix57k/matrix.mtx")
proc = subprocess.run(
    ["awk", f"$1=={rox2_idx} {{print $2, $3}}", mtx],
    capture_output=True, text=True, check=True,
)
rox2_counts = {}
for line in proc.stdout.splitlines():
    col, val = line.split()
    rox2_counts[int(col)] = rox2_counts.get(int(col), 0) + int(val)
print(f"cells with roX2 > 0 UMIs: {len(rox2_counts)}")

# 3. barcodes (1-indexed columns)
with open(os.path.join(GEO, "matrix57k/barcodes.tsv")) as fh:
    barcodes = [l.strip() for l in fh]
assert len(barcodes) == 56902

# 4. metadata sex reference
meta = pd.read_csv(os.path.join(GEO, "meta.tsv.gz"), sep="\t",
                   usecols=["new_barcode", "Age", "Genotype", "Replicate", "sex"])
meta["rox2_umi"] = 0
bc_to_col = {bc: i + 1 for i, bc in enumerate(barcodes)}
cols = meta["new_barcode"].map(bc_to_col)
meta["rox2_umi"] = cols.map(rox2_counts).fillna(0).astype(int)

# 5. classify and compare
meta["sex_rox2"] = (meta["rox2_umi"] > 0).map({True: "male", False: "female"})
ct = pd.crosstab(meta["sex"], meta["sex_rox2"], dropna=False)
print("\nmetadata sex (rows) vs roX2 classification (cols):")
print(ct)
agree = (meta["sex"] == meta["sex_rox2"]).mean()
print(f"\nagreement: {agree:.4f}")
print("\nroX2 UMI distribution by metadata sex:")
print(meta.groupby("sex")["rox2_umi"].describe()[["count", "mean", "max"]])
print("\nambiguous check — female cells with roX2>0:", ((meta.sex == "female") & (meta.rox2_umi > 0)).sum())
print("male cells with roX2==0:", ((meta.sex == "male") & (meta.rox2_umi == 0)).sum())

# 6. outputs
meta.to_csv(os.path.join(META, "cell_metadata_with_roX2_expression.csv"), index=False)
for sx in ("male", "female"):
    sub = meta[meta["sex_rox2"] == sx]
    sub[["new_barcode"]].to_csv(os.path.join(META, f"{sx}_barcodes.tsv"),
                                index=False, header=False)
    print(f"{sx}: {len(sub)} barcodes -> {sx}_barcodes.tsv")
print("\nWrote:", os.path.join(META, "cell_metadata_with_roX2_expression.csv"))
