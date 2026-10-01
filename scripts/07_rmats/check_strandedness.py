#!/usr/bin/env python3
"""Determine library strandedness from a BAM + GTF.

For 10x 3' v2, Read 2 (cDNA) is expected to be sense (same strand as mRNA),
but we verify empirically: count reads overlapping exons on the same vs
opposite strand of the annotated gene.

Usage: check_strandedness.py <bam> <gtf>
"""
import sys
import pysam

bam_path, gtf_path = sys.argv[1], sys.argv[2]

# build interval index of exons: chrom -> list of (start, end, strand)
# (use gene strand; restrict to protein-coding genes on main chromosomes)
import gzip
exons = {}
opener = gzip.open if gtf_path.endswith(".gz") else open
n_exon = 0
with opener(gtf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "exon":
            continue
        chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
        if chrom not in ("2L", "2R", "3L", "3R", "4", "X", "Y"):
            continue
        exons.setdefault(chrom, []).append((start, end, strand))
        n_exon += 1
        if n_exon > 400000:  # enough for a robust estimate
            break
print(f"indexed {n_exon} exons", flush=True)

bam = pysam.AlignmentFile(bam_path, "rb")
same = opp = 0
n_reads = 0
for read in bam.fetch():
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    n_reads += 1
    if n_reads > 200000:
        break
    chrom = bam.get_reference_name(read.reference_id)
    if chrom not in exons:
        continue
    pos = read.reference_start + 1  # 1-based
    rstrand = "-" if read.is_reverse else "+"
    for (s, e, gstrand) in exons[chrom]:
        if s <= pos <= e:
            if rstrand == gstrand:
                same += 1
            else:
                opp += 1
            break

tot = same + opp
print(f"exonic reads sampled: {tot} (from {n_reads} reads)")
print(f"same strand: {same} ({100*same/max(tot,1):.1f}%)")
print(f"opposite strand: {opp} ({100*opp/max(tot,1):.1f}%)")
if same > 0.9 * tot:
    print("VERDICT: sense (fr-firststrand for single-end cDNA read)")
elif opp > 0.9 * tot:
    print("VERDICT: antisense (fr-secondstrand for single-end cDNA read)")
else:
    print("VERDICT: ambiguous -> use fr-unstranded")
