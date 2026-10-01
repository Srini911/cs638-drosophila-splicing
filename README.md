# Sex-Specific Splicing in *Drosophila* Brain Aging (GSE107451)

CS638 project: Can machine learning identify sex-specific alternative-splicing signatures of brain aging?

**Result:** No — 10x 3′ scRNA-seq cannot detect SE splicing due to 3′ bias (95% PSI=1.0). ML accuracy 50% (chance). Negative result documented honestly.

## Quick Start

### 1. Data (already downloaded for 2 libraries)
- Raw sralites: `data/raw/sralite/`
- Split FASTQs: `data/raw/sra/SRR*/split/`
- Junctions: `data/intermediate/junctions/`
- PSI matrix: `data/intermediate/psi_matrix.tsv`

### 2. Run Analysis
```bash
# PSI (already computed)
python3 scripts/07_events/compute_psi.py

# Prepare ML matrix
python3 scripts/11_ml/prepare_ml_matrix.py

# Train logistic regression (leave-one-library-out CV)
python3 scripts/12_ml/train_logreg.py

# Differential splicing (paired)
python3 scripts/10_diff/diff_splicing.py

# Pathways
python3 scripts/13_pathways/pathway_enrichment.py

# Figures
python3 scripts/14_figures/make_figures.py
```

### 3. Results
- `results/ml/`: OOF predictions, metrics, coefficients.
- `results/diff_splicing/`: Paired male-vs-female results.
- `results/pathways/`: Enrichment test.
- `figures/main/`: 3 main figures.
- `REPORT.md`: Full report.

## Methods (Summary)

1. **Split FASTQs** by barcode/sex: `scripts/05_alignment/split_fastq_by_barcode.py`
2. **Align** with minimap2 (mappy), count junctions from CIGAR: `scripts/06_align/align_and_count.py`
3. **Define SE events** from GTF: `scripts/07_events/define_se_events.py`
4. **Compute PSI**: `scripts/07_events/compute_psi.py`
5. **ML**: Logistic regression, leave-one-library-out CV: `scripts/12_ml/train_logreg.py`

**Deviations from plan:** minimap2 instead of STAR (index build failed); custom Python instead of rMATS (build failed); no UMI dedup. See `docs/methods_notes.md`.

## Key Files

- `metadata/`: Barcodes, pseudobulk groups, cell metadata.
- `docs/`: Methods notes, dataset validation.
- `.gitignore`: Excludes FASTQ, BAM, SRA, temp files.

## Limitations

- n=4 pseudo-bulks (2 libraries at 3d); 6d/9d not processed (slow downloads).
- 95% PSI=1.0 due to 10x 3′ bias — fundamental limitation.
- SE events only; no A3SS/A5SS/MXE/RI.

## Author

Srinivas — CS638, 2026-10-01.
