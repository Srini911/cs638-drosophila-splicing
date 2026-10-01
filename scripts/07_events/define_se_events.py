#!/usr/bin/env python3
"""
Define skipped-exon (SE) events from GTF.
Usage: define_se_events.py
Output: metadata/se_events.tsv (event_id, gene, chrom, strand, e1_s, e1_e, e2_s, e2_e, e3_s, e3_e)
Event ID: gene:chrom:e1s-e1e:e2s-e2e:e3s-e3e:strand
"""
import gzip
from pathlib import Path
from collections import defaultdict

def parse_gtf(gtf_path):
    """Parse GTF, return dict transcript_id -> list of exons."""
    transcripts = defaultdict(list)
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            attrs = fields[8]
            # Extract transcript_id and gene_id/gene_name
            tid = gid = gname = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tid = attr.split('"')[1]
                elif attr.startswith('gene_id'):
                    gid = attr.split('"')[1]
                elif attr.startswith('gene_name'):
                    gname = attr.split('"')[1]
            if tid:
                transcripts[tid].append((chrom, start, end, strand, gid, gname or gid))
    return transcripts

def main():
    base = Path.home() / "workspace" / "ML_Ozgun"
    gtf = base / "data" / "reference" / "dmel_BDGP6.54.63.gtf"
    out = base / "metadata" / "se_events.tsv"

    print("Parsing GTF...", flush=True)
    transcripts = parse_gtf(gtf)
    print(f"Found {len(transcripts)} transcripts", flush=True)

    events = {}  # event_id -> (gene, chrom, strand, e1, e2, e3)
    for tid, exons in transcripts.items():
        if len(exons) < 3:
            continue
        # Sort by genomic position (accounting for strand)
        strand = exons[0][3]
        if strand == "+":
            exons_sorted = sorted(exons, key=lambda x: x[1])
        else:
            exons_sorted = sorted(exons, key=lambda x: x[1], reverse=True)
        # For each consecutive triplet, define SE event
        for i in range(len(exons_sorted) - 2):
            e1, e2, e3 = exons_sorted[i], exons_sorted[i+1], exons_sorted[i+2]
            chrom = e1[0]
            gene = e1[5]
            # Use genomic coordinates (1-based, inclusive)
            # For + strand: e1 is upstream, e3 is downstream
            # For - strand: e1 is downstream in genomic coords, but upstream in transcript
            # To avoid confusion, store in transcript order
            eid = f"{gene}:{chrom}:{e1[1]}-{e1[2]}:{e2[1]}-{e2[2]}:{e3[1]}-{e3[2]}:{strand}"
            if eid not in events:
                events[eid] = (gene, chrom, strand, e1[1], e1[2], e2[1], e2[2], e3[1], e3[2])

    print(f"Defined {len(events)} SE events", flush=True)
    with open(out, "w") as o:
        o.write("event_id\tgene\tchrom\tstrand\te1_s\te1_e\te2_s\te2_e\te3_s\te3_e\n")
        for eid, (gene, chrom, strand, e1s, e1e, e2s, e2e, e3s, e3e) in sorted(events.items()):
            o.write(f"{eid}\t{gene}\t{chrom}\t{strand}\t{e1s}\t{e1e}\t{e2s}\t{e2e}\t{e3s}\t{e3e}\n")
    print(f"Wrote {out}", flush=True)

if __name__ == "__main__":
    main()
