#!/usr/bin/env python3
"""
Pathway enrichment: test if variable SE events are enriched in neurodegeneration pathways.
Uses top variable events (by PSI std dev) vs background (all tested events).
Simple Fisher's exact test with a curated gene set.
Usage: pathway_enrichment.py
"""
import csv
from pathlib import Path
from scipy.stats import fisher_exact

# Curated Drosophila neurodegeneration-related genes (from literature)
# This is a small illustrative set; a full analysis would use GO/KEGG.
NEURODEG_GENES = {
    "tau", "Appl", "park", "PINK1", "DJ-1beta", "LRRK2", "SNCA",
    "HTT", "ATXN1", "ATXN3", "SOD1", "TDP-43", "FUS",
    "presenilin", "Psn", "nicastrin", "Aph-1",
    "parkin", "Pink1",
}

def gene_from_event(event_id):
    """Extract gene symbol from event_id (format: GENE:chrom:start-end:...)."""
    return event_id.split(":")[0]

def main():
    base = Path.home() / "workspace" / "ML_Ozgun"
    infile = base / "data" / "intermediate" / "psi_matrix.tsv"
    outdir = base / "results" / "pathways"
    outdir.mkdir(parents=True, exist_ok=True)

    # Read PSI matrix, compute std dev per event
    import numpy as np
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        samples = [c for c in reader.fieldnames if c != "event_id"]
        event_stds = []
        for row in reader:
            eid = row["event_id"]
            vals = [float(row[s]) for s in samples if row[s] != "NA"]
            if len(vals) >= 2:
                std = np.std(vals)
                event_stds.append((eid, std, gene_from_event(eid)))

    # Top 10% most variable
    event_stds.sort(key=lambda x: x[1], reverse=True)
    n_top = max(10, len(event_stds) // 10)
    top_events = event_stds[:n_top]
    all_genes = set(g for _, _, g in event_stds)
    top_genes = set(g for _, _, g in top_events)

    print(f"Total events: {len(event_stds)}, top variable: {n_top}", flush=True)
    print(f"Total genes: {len(all_genes)}, top genes: {len(top_genes)}", flush=True)

    # Fisher's exact test for neurodegeneration enrichment
    # Contingency: [in_top_and_neuro, in_top_not_neuro; not_top_and_neuro, not_top_not_neuro]
    top_neuro = len(top_genes & NEURODEG_GENES)
    top_not = len(top_genes - NEURODEG_GENES)
    bg_neuro = len((all_genes - top_genes) & NEURODEG_GENES)
    bg_not = len((all_genes - top_genes) - NEURODEG_GENES)

    table = [[top_neuro, top_not], [bg_neuro, bg_not]]
    odds, p = fisher_exact(table)

    print(f"\nNeurodegeneration gene enrichment:", flush=True)
    print(f"  Top variable genes in neuro set: {top_neuro}/{len(top_genes)}", flush=True)
    print(f"  Background genes in neuro set: {bg_neuro}/{len(all_genes - top_genes)}", flush=True)
    print(f"  Odds ratio: {odds:.2f}, p-value: {p:.3g}", flush=True)

    # Write results
    with open(outdir / "pathway_enrichment.txt", "w") as o:
        o.write(f"Top variable SE events: {n_top} / {len(event_stds)}\n")
        o.write(f"Neurodegeneration gene set size: {len(NEURODEG_GENES)}\n")
        o.write(f"Top genes overlapping: {top_neuro}\n")
        o.write(f"Odds ratio: {odds:.3f}\n")
        o.write(f"Fisher p-value: {p:.3g}\n")
        o.write(f"\nInterpretation: {'Enriched' if p < 0.05 else 'Not enriched'} ")
        o.write(f"({'p<0.05' if p < 0.05 else 'p>=0.05'})\n")

    # Also list top variable genes
    with open(outdir / "top_variable_genes.csv", "w", newline="") as o:
        w = csv.writer(o)
        w.writerow(["event_id", "gene", "psi_std"])
        for eid, std, gene in top_events[:50]:
            w.writerow([eid, gene, f"{std:.4f}"])

    print(f"\nWrote {outdir}/pathway_enrichment.txt", flush=True)

if __name__ == "__main__":
    main()
