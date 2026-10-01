#!/usr/bin/env python3
"""
Compute PSI for SE events from junction counts.
Usage: compute_psi.py
Input:  metadata/se_events.tsv, data/intermediate/junctions/*_junctions.tsv
Output: data/intermediate/psi_matrix.tsv (events x samples), data/intermediate/psi_qc.tsv
PSI = (I1 + I2) / (I1 + I2 + 2*S)  [average inclusion, see notes]
Requires min_total >= 10 junction reads per event per sample.
"""
import sys
from pathlib import Path
from collections import defaultdict

def load_junctions(junc_file):
    """Return dict (chrom, donor, acceptor) -> count (strand-agnostic for matching)."""
    d = {}
    with open(junc_file) as f:
        next(f)  # header
        for line in f:
            chrom, donor, acceptor, strand, count = line.strip().split("\t")
            # Use (chrom, donor, acceptor) as key; sum over strands (should be consistent)
            key = (chrom, int(donor), int(acceptor))
            d[key] = d.get(key, 0) + int(count)
    return d

def main():
    base = Path.home() / "workspace" / "ML_Ozgun"
    events_file = base / "metadata" / "se_events.tsv"
    junc_dir = base / "data" / "intermediate" / "junctions"
    out_matrix = base / "data" / "intermediate" / "psi_matrix.tsv"
    out_qc = base / "data" / "intermediate" / "psi_qc.tsv"

    # Load events
    events = []
    with open(events_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            eid, gene, chrom, strand = parts[0], parts[1], parts[2], parts[3]
            e1s, e1e, e2s, e2e, e3s, e3e = map(int, parts[4:10])
            # Junctions (1-based): donor = first base of intron = e1_e+1, acceptor = first base of next exon = e2_s
            j_inc1 = (chrom, e1e+1, e2s)
            j_inc2 = (chrom, e2e+1, e3s)
            j_skip = (chrom, e1e+1, e3s)
            events.append((eid, gene, j_inc1, j_inc2, j_skip))

    print(f"Loaded {len(events)} events", flush=True)

    # Load junction counts per sample
    junc_files = sorted(junc_dir.glob("*_junctions.tsv"))
    samples = [f.stem.replace("_junctions", "") for f in junc_files]
    print(f"Found {len(samples)} samples: {samples}", flush=True)

    # Compute PSI per event per sample
    # psi[event][sample] = PSI or None
    psi = defaultdict(dict)
    qc_stats = []
    for jf, sample in zip(junc_files, samples):
        jd = load_junctions(jf)
        n_computed = 0
        for eid, gene, j1, j2, j3 in events:
            # Try exact match first, then with +1 offset (to handle coordinate confusion)
            c1 = jd.get(j1, 0)
            c2 = jd.get(j2, 0)
            c3 = jd.get(j3, 0)
            total = c1 + c2 + c3
            if total >= 10:
                # PSI = inclusion / (inclusion + skipping)
                # Inclusion = (c1 + c2) / 2 (average of two junctions)
                inc = (c1 + c2) / 2.0
                psi_val = inc / (inc + c3) if (inc + c3) > 0 else None
                if psi_val is not None:
                    psi[eid][sample] = psi_val
                    n_computed += 1
        qc_stats.append((sample, n_computed))
        print(f"[{sample}] {n_computed} events with PSI", flush=True)

    # Write matrix (events with PSI in at least 50% of samples)
    min_samples = max(2, len(samples) // 2)
    kept_events = [eid for eid in psi if len(psi[eid]) >= min_samples]
    print(f"Keeping {len(kept_events)} events (PSI in >={min_samples} samples)", flush=True)

    with open(out_matrix, "w") as o:
        o.write("event_id\t" + "\t".join(samples) + "\n")
        for eid in sorted(kept_events):
            vals = [str(round(psi[eid].get(s, float("nan")), 4)) if s in psi[eid] else "NA" for s in samples]
            o.write(eid + "\t" + "\t".join(vals) + "\n")

    with open(out_qc, "w") as o:
        o.write("sample\tn_events_with_psi\n")
        for s, n in qc_stats:
            o.write(f"{s}\t{n}\n")

    print(f"Wrote {out_matrix} ({len(kept_events)} events x {len(samples)} samples)", flush=True)

if __name__ == "__main__":
    main()
