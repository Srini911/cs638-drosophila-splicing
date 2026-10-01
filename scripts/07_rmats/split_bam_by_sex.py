#!/usr/bin/env python3
"""Phase 7 — split a STARsolo genomic BAM into sex-specific pseudo-bulk BAMs.

Reads are assigned by their CB tag (corrected cell barcode). Barcodes come from
the validated per-library sex lists (metadata/barcodes/<lib>_<sex>.txt), which
already exclude the roX2>=2 suspicious cells.

Usage: split_bam_by_sex.py <starsolo_bam> <lib_label> <outdir>
Writes: <outdir>/<lib_label>_male.bam(.bai), <outdir>/<lib_label>_female.bam(.bai)
and a small assignment summary.
"""
import os
import sys

import pysam

bam_in, lib, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)
base = os.path.expanduser("~/workspace/ML_Ozgun")

def load_bcs(sex):
    p = os.path.join(base, "metadata/barcodes", f"{lib}_{sex}.txt")
    return set(l.strip() for l in open(p) if l.strip())

male_bcs = load_bcs("male")
female_bcs = load_bcs("female")
print(f"{lib}: {len(male_bcs)} male, {len(female_bcs)} female barcodes", flush=True)

src = pysam.AlignmentFile(bam_in, "rb")
out_m = pysam.AlignmentFile(os.path.join(outdir, f"{lib}_male.bam"), "wb", template=src)
out_f = pysam.AlignmentFile(os.path.join(outdir, f"{lib}_female.bam"), "wb", template=src)
n_m = n_f = n_un = 0
for read in src:
    try:
        cb = read.get_tag("CB")
    except KeyError:
        n_un += 1
        continue
    # STARsolo may append "-1" suffix to CBs; strip it
    cb = cb.split("-")[0]
    if cb in male_bcs:
        out_m.write(read); n_m += 1
    elif cb in female_bcs:
        out_f.write(read); n_f += 1
    else:
        n_un += 1
out_m.close(); out_f.close(); src.close()
for sx in ("male", "female"):
    pysam.index(os.path.join(outdir, f"{lib}_{sx}.bam"))
total = n_m + n_f + n_un
print(f"{lib}: male_reads={n_m} female_reads={n_f} unassigned={n_un} "
      f"(assigned {100.0*(n_m+n_f)/max(total,1):.1f}%)", flush=True)
with open(os.path.join(outdir, f"{lib}_split_summary.txt"), "w") as fh:
    fh.write(f"male_reads={n_m}\nfemale_reads={n_f}\nunassigned={n_un}\n")
