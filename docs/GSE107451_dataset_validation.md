# GSE107451 Dataset Validation — Phase 2

**Date:** 2026-09-30/2026-10-01
**Validator:** automated pipeline (all facts below verified from authoritative sources, not assumed)

## Study
- **Accession:** GEO GSE107451; BioProject PRJNA420091; SRA study SRP125768
- **Publication:** Davie et al., Cell 2018, "A single-cell transcriptome atlas of the ageing
  Drosophila brain" (PMC6086935)
- **Objective:** single-cell transcriptomic atlas of the aging Drosophila brain across lifespan

## Verified characteristics (from GEO metadata file
`GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv.gz` + SRA records)

| Attribute | Value |
|---|---|
| Organism | *Drosophila melanogaster* (taxid 7227) |
| Tissue | Whole brain (dissected, cuticle/retina/trachea removed) |
| Genotypes | w1118 and DGRP-551 |
| Ages (days) | 0, 1, 3, 6, 9, 15, 30, 50 |
| Assay | 10x Genomics Chromium Single Cell 3' v2 (per SRA library protocol) |
| SRA instruments | Illumina HiSeq 4000 / NextSeq 500, paired-end |
| Sex annotation | Per-cell `sex` column in GEO metadata (male/female) |
| Total cells (57k matrix) | 56,902 |

## Focus subset for this project (verified, not assumed)
- w1118, ages 3d / 6d / 9d: **16,019 cells**
- 3 ten-x libraries per age (replicates), each containing both sexes:

| Age | Libraries (SRA) | Female cells | Male cells |
|---|---|---|---|
| 3d | SRR6327125 (r1), SRR6327126 (r2), SRR6327127 (r3) | 477+1382+1197=3056 | 508+903+724=2135 |
| 6d | SRR6327128 (r1), SRR6327129 (r2), SRR6327130 (r3) | 623+611+1987=3221 | 686+839+939=2464 |
| 9d | SRR6327131 (r1), SRR6327132 (r2), SRR6327133 (r3) | 624+797+892=2313 | 328+1143+1331=2802 |

Library↔SRA mapping verified via SRA experiment titles (`w1118_3d_r1` etc.) and
metadata barcode suffixes (`<16bp>_w1118_3d_r1`).

## Sex representation
- Sex is explicitly provided per cell in the authoritative GEO metadata.
- Independently validated via roX2 (FBgn0019660, male-specific lncRNA):
  88.2% agreement; disagreements characterized (male dropout ~18% due to scRNA-seq
  sparsity; 286/56,902 ≈ 0.5% metadata-female cells with roX2 ≥ 2 UMIs excluded as
  probable mislabels/doublets). See `scripts/06_sex_classification/validate_rox2.py`.
- Final: 18 pseudo-bulk groups (3 ages × 3 replicates × 2 sexes); group table:
  `metadata/pseudobulk_groups.csv`.

## Ages retained / excluded
- **Retained:** 3d, 6d, 9d (young-adult aging window; matches prior course work; each
  with 3 replicates × 2 sexes).
- **Excluded:** 0d, 1d (newly eclosed; developmental, not aging), 15d, 30d, 50d
  (kept out of scope to bound compute for the 2-day deadline; documented here so the
  choice is explicit, not silent).

## Reference requirements
- Assembly: BDGP Release 6 (all r6.xx annotations share coordinates).
- Annotation used: Ensembl BDGP6.54 (FlyBase-sourced gene models, FBgn IDs).
  FlyBase r6.64 GTF was unreachable from this compute environment (FTP blocked by
  egress proxy); substitution documented in `docs/reference_manifest.md`.

## Limitations (carried into report)
1. 10x 3' v2 chemistry: reads are 3'-biased; splice-junction coverage is sparse and
   concentrated near transcript 3' ends. rMATS results are hypothesis-generating;
   every event requires read-support QC (Phase 8) before interpretation.
2. Pseudo-bulk observational unit is the library×sex group (n=18), NOT the cell.
   Thousands of cells ≠ thousands of independent samples.
3. Sex labels are author-provided + roX2-validated, not experimentally re-derived.
4. Only w1118 (not DGRP-551) and 3 ages used — limits generalizability claims.
