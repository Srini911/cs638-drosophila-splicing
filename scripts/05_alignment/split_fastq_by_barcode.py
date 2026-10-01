#!/usr/bin/env python3
"""
Split 10x FASTQs by cell barcode (using whitelist + sex assignments).
Backup approach if STARsolo on full FASTQs is too slow.

Input:  _1.fastq (barcodes), _2.fastq (cDNA), barcode lists
Output: male_1.fastq, male_2.fastq, female_1.fastq, female_2.fastq (only whitelist barcodes)

Usage: split_fastq_by_barcode.py SRR LIB_LABEL
"""
import sys, gzip
from pathlib import Path

def open_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def main():
    srr, lib = sys.argv[1], sys.argv[2]
    base = Path.home() / "workspace" / "ML_Ozgun"
    fqdir = base / "data" / "raw" / "sra" / srr
    metadir = base / "metadata" / "barcodes"

    # Load barcode -> sex mapping
    bc_to_sex = {}
    for sex in ["male", "female"]:
        bc_file = metadir / f"{lib}_{sex}.txt"
        if bc_file.exists():
            for line in open(bc_file):
                bc = line.strip()
                if bc:
                    bc_to_sex[bc] = sex

    print(f"[{lib}] loaded {len(bc_to_sex)} barcodes", flush=True)

    # Find input FASTQs
    r1 = fqdir / f"{srr}_1.fastq"
    r2 = fqdir / f"{srr}_2.fastq"
    if not r1.exists():
        r1 = fqdir / f"{srr}_1.fastq.gz"
        r2 = fqdir / f"{srr}_2.fastq.gz"

    outdir = fqdir / "split"
    outdir.mkdir(exist_ok=True)
    out_handles = {}
    for sex in ["male", "female"]:
        out_handles[sex] = (
            open(outdir / f"{sex}_1.fastq", "w"),
            open(outdir / f"{sex}_2.fastq", "w"),
        )

    # Process reads
    total, kept = 0, 0
    with open_fastq(r1) as f1, open_fastq(r2) as f2:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline(); p1 = f1.readline(); q1 = f1.readline()
            h2 = f2.readline(); s2 = f2.readline(); p2 = f2.readline(); q2 = f2.readline()

            total += 1
            bc = s1.strip()[:16]  # 16bp cell barcode
            sex = bc_to_sex.get(bc)
            if sex:
                o1, o2 = out_handles[sex]
                o1.write(h1); o1.write(s1); o1.write(p1); o1.write(q1)
                o2.write(h2); o2.write(s2); o2.write(p2); o2.write(q2)
                kept += 1

            if total % 1000000 == 0:
                print(f"[{lib}] processed {total:,} reads, kept {kept:,} ({kept/total*100:.1f}%)", flush=True)

    for o1, o2 in out_handles.values():
        o1.close(); o2.close()

    print(f"[{lib}] DONE: {total:,} total, {kept:,} kept ({kept/total*100:.1f}%)", flush=True)

if __name__ == "__main__":
    main()
