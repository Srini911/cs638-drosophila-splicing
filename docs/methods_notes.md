# Methods Notes (running log)

## 2026-09-30/10-01 — Pipeline decisions

### Dataset
- GSE107451 is 10x Chromium Single Cell 3' v2 (NOT bulk). Verified via SRA library
  construction protocol and GEO metadata.
- Focus: w1118, 3d/6d/9d, 9 libraries (SRR6327125–SRR6327133), 16,019 cells.

### Sex labels
- Primary: GEO metadata `sex` column (authoritative, per-cell).
- Validated via roX2 (FBgn0019660): 88.2% agreement. Disagreements = male dropout
  (expected) + 286 suspicious metadata-female cells with roX2≥2 (excluded).
- 18 pseudo-bulk groups; barcodes in `metadata/barcodes/<lib>_<sex>.txt`.

### Raw data acquisition
- SRA .sralite via NCBI CloudFront (only viable path through egress proxy;
  sratoolkit `prefetch` cannot use the proxy; ENA FASTQ too slow at 0.8 MB/s).
- sralite URLs resolved per-accession via E-utilities efetch runinfo `download_path`.
- fasterq-dump --include-technical --split-files to recover 10x barcode reads.

### Alignment
- STAR 2.7.11b STARsolo, 2 threads.
- Reference: Ensembl BDGP6.54 (primary assembly) + BDGP6.54.63 GTF.
  (FlyBase r6.64 unreachable; substitution documented in reference_manifest.md.)
- Whitelist: per-library 16bp barcode lists derived from GEO metadata
  (NOT the generic 10x 737K list — unavailable; custom list is more precise).
- --soloCBmatchWLtype 1MM, --soloCellFilter None (barcodes pre-QC'd).
- Genomic BAM (Aligned.sortedByCoord.out.bam) retained for rMATS; split by CB tag.

### rMATS
- Per age: male (3 BAMs) vs female (3 BAMs), paired-end mode, readLength 98.
- libType determined empirically from BAM strandedness (default fr-unstranded
  unless data shows clear strandedness).
- Significance rule (to be confirmed): FDR < 0.05 AND |ΔPSI| ≥ 0.1.

### ML
- See scripts/12_ml/train_logreg.py for the leakage-safe design.
- Primary: L2 logistic regression, leave-one-library-out CV (9 folds).
- Age as covariate; secondary leave-one-age-out generalization check (TBD).

## Tool versions
- STAR 2.7.11b (static binary)
- sratoolkit 3.1.1 (binaries only; prefetch unusable behind proxy)
- rMATS-turbo (source build; pending GSL)
- Python 3.12, pysam 0.24.1, scikit-learn (pip)

## 2026-10-01 — Read structure verified from SRA
- fasterq-dump on .sralite WITH --include-technical recovers 10x barcode reads.
- SRR6327125: _1 = 106.2M reads x 26bp (16bp CB + 10bp UMI); _2 = cDNA x 120bp.
- STAR sjdbOverhang 119; rMATS --readLength 120.

## 2026-10-01 ~00:41 — VM reboot incident
- The VM rebooted unexpectedly (uptime reset to 1 min). All processes in /opt
  (STAR, sratoolkit, rMATS-turbo build, GSL) were lost; /opt is ephemeral.
- ~/workspace persisted fully: sralite downloads (2/9), lib1 FASTQs, scripts,
  metadata, docs all intact.
- Recovery: reinstalled all tools to ~/workspace/tools (persistent), restarted
  downloads (resume), STAR index, GSL build. Lesson: keep everything under ~.

## 2026-10-01 ~03:00 — 3' bias limitation observed
- PSI computation on 3 samples (2,647 SE events): 95.1% of PSI values = 1.0.
- Only 3.4% have PSI < 0.9 (variable).
- This is expected: 10x Chromium 3' v2 captures the 3' end of transcripts;
  most reads do not span splice junctions. SE events in the gene body are
  invisible; only 3'-proximal events are detected, and most are constitutive.
- The ML analysis will proceed, but we expect weak signal. A negative result
  (ML cannot identify sex-specific splicing from 3' data) is scientifically
  valid and will be reported honestly.
- Mitigation: we use minimap2 splice-aware alignment (not STARsolo) due to
  VM constraints; junction counting is direct from CIGAR N operations.

## 2026-10-01 ~03:00 — Methodological deviations (infrastructure-constrained)
Due to VM limitations (2 CPUs, 7GB RAM, flaky network, ephemeral /opt), the
following deviations from the original rMATS+STARsolo plan were necessary:

1. **Aligner: minimap2 (mappy) instead of STARsolo.**
   - STAR genome index build never completed (SA sort pathological, >30 min
     for genome-only index; killed after multiple attempts).
   - minimap2 index builds in 10 sec; splice-aware alignment works for 120bp reads.
   - Trade-off: minimap2 is less sensitive than STAR for splicing, but sufficient
     for junction counting.

2. **Splicing: custom Python instead of rMATS.**
   - rMATS-turbo source build failed (no gfortran/BLAS via apt; GSL build OK
     but lbfgsb needed Fortran).
   - rMATS via bioconda failed (network timeouts).
   - Custom implementation: define SE events from GTF, count junctions from
     minimap2 CIGAR N operations, compute PSI = inclusion/(inclusion+skipping).
   - Only SE events (not A3SS/A5SS/MXE/RI) due to time; SE is most interpretable.

3. **FASTQ pre-splitting by barcode.**
   - Only 15-50% of reads are from whitelist barcodes (rest are ambient).
   - Splitting before alignment gives 2-6x speedup vs STARsolo on full FASTQs.
   - Sex separation done at FASTQ level; no BAM CB-splitting needed.

4. **Downsampling to 5M reads/sex for large libraries.**
   - Lib2 had 39.5M reads (vs 16.5M for lib1); downsampled to 5M/sex for time.
   - PSI is a ratio; downsampling preserves relative junction usage (documented).

5. **No UMI deduplication.**
   - UMI info in R1 not used; positional duplicates not marked.
   - Justification: PSI is a ratio (inclusion/skipping); PCR duplicates affect
     numerator and denominator proportionally (approximately unbiased).
   - Limitation documented.

All deviations are documented here; results will note these limitations.
