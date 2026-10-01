#!/usr/bin/env python3
"""
Align split FASTQs with minimap2 (mappy) and extract splice junctions.
Usage: align_and_count.py SRR LIB_LABEL SEX
Output: data/intermediate/junctions/<LIB>_<SEX>_junctions.tsv
        (chrom, donor, acceptor, strand, count)
"""
import sys
import mappy
from pathlib import Path
from collections import Counter

def main():
    srr, lib, sex = sys.argv[1], sys.argv[2], sys.argv[3]
    # Optional 4th arg: input FASTQ suffix (e.g., "_ds5m" for downsampled)
    suffix = sys.argv[4] if len(sys.argv) > 4 else ""
    base = Path.home() / "workspace" / "ML_Ozgun"
    fq = base / "data" / "raw" / "sra" / srr / "split" / f"{sex}_2{suffix}.fastq"
    outdir = base / "data" / "intermediate" / "junctions"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{lib}_{sex}_junctions.tsv"

    print(f"[{lib} {sex}] building minimap2 index...", flush=True)
    ref = str(base / "data" / "reference" / "dmel_BDGP6.54.dna.primary_assembly.fa")
    aligner = mappy.Aligner(ref, preset="splice")
    assert aligner, "failed to build index"

    junctions = Counter()
    total, mapped, spliced = 0, 0, 0

    print(f"[{lib} {sex}] aligning...", flush=True)
    with open(fq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip()
            f.readline(); f.readline()  # plus, qual
            total += 1

            hits = list(aligner.map(seq))
            if not hits:
                continue
            # Take primary alignment (highest mapq)
            aln = max(hits, key=lambda x: x.mapq)
            if aln.mapq < 10:  # low quality, skip
                continue
            mapped += 1

            # Parse CIGAR for N (spliced) operations
            # mappy cigar is list of [length, op]; op 3 = N (refskip/intron)
            ref_pos = aln.r_st
            for ln, op in aln.cigar:
                if op == 3:  # N: splice junction
                    # Filter: require intron >= 50bp (Drosophila introns are typically >50bp)
                    # This removes false-positive micro-indels mis-called as splicing
                    if ln >= 50:
                        donor = ref_pos
                        acceptor = ref_pos + ln
                        # Record junction (1-based, like GTF)
                        junctions[(aln.ctg, donor+1, acceptor+1, aln.strand)] += 1
                        spliced += 1
                if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                    ref_pos += ln

            if total % 1000000 == 0:
                print(f"[{lib} {sex}] {total:,} reads, {mapped:,} mapped, {spliced:,} junctions", flush=True)

    print(f"[{lib} {sex}] DONE: {total:,} reads, {mapped:,} mapped, {len(junctions):,} unique junctions", flush=True)
    with open(out, "w") as o:
        o.write("chrom\tdonor\tacceptor\tstrand\tcount\n")
        for (chrom, d, a, s), c in junctions.most_common():
            o.write(f"{chrom}\t{d}\t{a}\t{s}\t{c}\n")
    print(f"[{lib} {sex}] wrote {out}", flush=True)

if __name__ == "__main__":
    main()
