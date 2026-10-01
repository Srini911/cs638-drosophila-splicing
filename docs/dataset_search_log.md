# Dataset Search Log — Phase 1

**Date:** 2026-09-30
**Objective:** Find public datasets capable of addressing sex-specific alternative splicing
in Drosophila brain aging, then justify (or reject) the incumbent GSE107451.

## Inclusion criteria

A dataset is INCLUDED as a candidate only if it plausibly satisfies ALL of:

1. Organism: *Drosophila melanogaster*
2. Tissue: brain (whole brain or dissected brain; head-only is marginal)
3. Assay: RNA-seq (bulk preferred for splicing; single-cell acceptable with caveats)
4. Sex: known per sample, OR scientifically defensible inference (e.g., roX2/sex-marker expression)
5. Aging dimension: ≥3 adult age points spanning young → aged (marginal: 2 points)
6. Replication: biological replicates per group (marginal: n=2)
7. Data: raw reads or alignment-level data obtainable (SRA/GEO)
8. Metadata sufficient to reconstruct experimental groups (age × sex × genotype × replicate)

## Exclusion criteria (any one excludes)

- No sex information and no defensible recovery strategy
- Pooled sexes with no recovery strategy
- No aging dimension (single time point, or young-adult-only range)
- Wrong tissue without brain data
- Incompatible design for rMATS (e.g., 3'-biased single-nucleus with no junction coverage — assessed per dataset)
- Raw data unavailable / unclear experimental design

## Search strategy

1. Web/GEO search for "Drosophila brain aging RNA-seq GEO male female".
2. Follow citations from review/frontier articles on fly brain aging transcriptomics.
3. For each candidate: verify organism, tissue, ages, sexes, assay, replication, and raw-data
   availability from the publication and repository record (GEO/OmicsDI/GitHub mirrors).
4. Record decision + reason in `results/dataset_candidates.csv`.

Searches performed 2026-09-30 (see `results/dataset_candidates.csv` for per-candidate detail):

- `GSE107451 Drosophila brain aging RNA-seq GEO`
- `Drosophila melanogaster brain aging sex-specific RNA-seq GEO dataset male female`
- `Pacifico Davis 2018 aging Drosophila brain RNA-seq GEO accession brain transcriptome olfactory memory`

## Key finding (verified from multiple secondary sources)

**GSE107451 (Davie et al., Cell 2018, "A single-cell transcriptome atlas of the ageing
Drosophila brain") is 10x Chromium SINGLE-CELL RNA-seq — not bulk RNA-seq.**

Verified facts (PMC6086935; OmicsDI record; independent GitHub mirror hjkoh93/drosophila-grx2-aging):

- Assay: 10x Chromium scRNA-seq (3' tag-based), ~56,902 cells
- Tissue: whole brain; strains w1118 and DGRP-551
- Ages: 0, 1, 3, 6, 9, 15, 30, 50 days
- Sexes: male and female brains sequenced (paper reports male/female reproducibility across
  clusters); per-cell sex NOT annotated in GEO metadata → requires inference (roX2, FBgn0019660)
- GEO files: `GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv.tar`
  + `GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv`; BioProject PRJNA420091

### Implication for the project (recorded as the central methodological risk)

10x 3' scRNA-seq reads are strongly 3'-biased and sparse per cell. rMATS (a bulk
splice-junction method) can only be applied after pseudo-bulking sex/age/replicate groups,
and junction coverage will be thin and 3'-biased. Consequences adopted as project policy:

- Junction/read-support QC is MANDATORY before interpreting any PSI value (Phase 8).
- rMATS on this data is treated as exploratory/hypothesis-generating, stated as such in the report.
- The bulk Pacifico et al. 2018 dataset (below) is the natural validation dataset for any
  splicing finding, because bulk RNA-seq has full-transcript junction coverage.

## Candidate dispositions (detail in `results/dataset_candidates.csv`)

| Accession | Verdict | One-line reason |
|---|---|---|
| GSE107451 (Davie et al. 2018) | INCLUDE (primary, with caveats) | Rich age series (8 points), both sexes, brain; scRNA-seq 3' bias limits rMATS → pseudo-bulk + strict junction QC required |
| Pacifico et al. 2018, PLoS ONE (bulk; GEO acc. TBD from paper) | INCLUDE (validation/secondary) | Bulk RNA-seq, dissected brain, Canton-S, M+F, ages 5/20/30/40d — best rMATS-compatible design; small n at some points |
| Lagging brain gene expression, Mol Neurobiol 2024 (bulk, 3/7/14d M+F) | EXCLUDE | Aging range too narrow (max 14d = young adult) |
| Fly Cell Atlas head (Lu et al. 2023; sn/scRNA-seq, 5/30/50/70d) | EXCLUDE | Head (not brain) + single-nucleus 3' assay; same splicing limitations, less relevant tissue |
| Bulk head aging, 15 time points, males only | EXCLUDE | No sex dimension — cannot address the research question |

## Decision

Proceed with **GSE107451 as the primary dataset** (continuity with prior course work; 3d/6d/9d
w1118 focus per brief §5), with the 3'-bias limitation documented as the project's top
methodological risk and Pacifico et al. 2018 held as the bulk validation dataset.
Phase 2 will verify GSE107451's metadata (ages, strains, sex recoverability) from the
authoritative GEO metadata file before any analysis.
