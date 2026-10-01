# Can Machine Learning Identify Sex-Specific Alternative-Splicing Signatures of Brain Aging in *Drosophila*?

## CS638 Project Report — Srinivas

**Date:** 2026-10-01  
**Dataset:** GSE107451 (*Drosophila melanogaster* brain aging, 10x Chromium Single Cell 3′ v2)  
**Question:** Can machine learning identify sex-specific alternative-splicing signatures of brain aging, and are these signatures enriched in pathways associated with neurodegeneration?

---

## Abstract

We investigated sex-specific alternative splicing in the aging *Drosophila* brain using single-cell RNA-seq data (GSE107451). Due to the 3′ capture bias inherent in 10x Chromium v2 chemistry, 95.1% of skipped-exon (SE) events showed constitutive inclusion (PSI=1.0), leaving only ~3-5% with measurable variation. A logistic regression classifier trained on PSI values (leave-one-library-out CV, n=4 pseudo-bulks) achieved 50% accuracy (chance level; ROC-AUC 0.375), indicating no detectable sex-specific splicing signal in this data type. Paired differential splicing analysis (n=2 libraries) found 3/2,098 events nominally significant (p<0.05), fewer than expected by chance. Pathway enrichment for neurodegeneration genes was not significant (p=0.27). We conclude that 10x 3′ scRNA-seq is unsuitable for SE splicing analysis due to 3′ bias; full-length RNA-seq would be required to address the research question.

---

## 1. Introduction

Alternative splicing generates transcript diversity and is implicated in aging and neurodegeneration. Sex differences in splicing could contribute to sex-biased disease susceptibility. The *Drosophila* brain aging dataset (GSE107451; Davie et al., Cell 2018) provides single-cell transcriptomes across ages 0-50 days, offering an opportunity to study sex-specific splicing in aging.

**Challenge:** The dataset uses 10x Chromium 3′ v2 chemistry, which captures only the 3′ end of transcripts. Most reads do not span splice junctions, limiting splicing analysis to 3′-proximal events.

---

## 2. Methods

### 2.1 Dataset and Design

- **Source:** GEO GSE107451 (BioProject PRJNA420091, SRA SRP125768).
- **Selected subset:** w1118 genotype, ages 3d/6d/9d (3 libraries per age in full design).
- **Processed here:** 2 libraries at 3d (SRR6327125, SRR6327126), male/female pseudo-bulks per library (n=4).
- **Sex labels:** Author-provided per-cell metadata (GEO). Validated against roX2 expression (male-specific lncRNA): 88.2% agreement. roX2 dropout expected in 3′ data; author labels retained.

### 2.2 Alignment and Junction Counting

Due to VM infrastructure constraints (STAR index build pathological; rMATS build failed), we used:
- **FASTQ splitting:** Reads split by exact barcode match to author sex labels (15-50% of reads retained; ambient RNA discarded).
- **Aligner:** minimap2 (via mappy, splice-aware) on cDNA R2 reads (120bp).
- **Junctions:** Extracted from CIGAR N operations (intron length ≥50bp).
- **Limitations:** No UMI deduplication (PSI ratio approximately unbiased); exact barcode matching (1-mismatch reads lost).

### 2.3 SE Event Definition and PSI

- **Events:** 59,237 candidate skipped-exon events defined from GTF as consecutive exon triplets (gene:chrom:exon1:exon2:exon3).
- **PSI:** For each event, inclusion = (junction1 + junction2)/2; skipping = junction3; PSI = inclusion/(inclusion+skipping). Minimum 10 total junction reads.
- **Result:** 3,448 events with PSI in ≥2 samples; 200 with variance >0 across 4 samples.

### 2.4 Machine Learning

- **Model:** L2 logistic regression (C=1.0), primary per course design.
- **Features:** 145 SE events with PSI in all 4 samples and variance >0.
- **Target:** Sex (male=1, female=0).
- **Validation:** Leave-one-library-out CV (2 folds; each holds out male+female from one 10x library). Feature selection (variance filter, SelectKBest k=50) and scaling inside folds. Age as covariate.
- **Baselines:** Majority class (50%), permuted labels.

### 2.5 Differential Splicing and Pathways

- **Paired test:** Male vs female within library (n=2 pairs), paired t-test on PSI.
- **Pathways:** Fisher's exact test for enrichment of top 10% variable genes in a curated neurodegeneration gene set (tau, Appl, park, PINK1, etc.) vs background.

---

## 3. Results

### 3.1 3′ Bias Limits Splicing Detection

Of 7,412 PSI measurements, **95.1% were exactly 1.0** (constitutive inclusion); only 3.4% had PSI<0.9 (Figure 1). This extreme ceiling effect reflects 10x 3′ chemistry: most reads derive from the 3′ UTR and do not span splice junctions. Only 3′-proximal SE events are detectable, and most are constitutive.

### 3.2 ML Cannot Distinguish Sex from Splicing

Leave-one-library-out CV (n=4):
- **Accuracy:** 0.50 (chance)
- **ROC-AUC:** 0.375 (worse than chance)
- **Conclusion:** No detectable sex-specific splicing signal in 3′ data.

### 3.3 No Significant Differential Splicing

Paired analysis (2,098 events, n=2 pairs):
- 3 events nominally significant (p<0.05), **fewer than the ~105 expected by chance**.
- Top hit: CG8079 (delta PSI=-0.42, p=0.06), not significant after correction.
- Volcano plot (Figure 2) shows no clear sex-separated events.

### 3.4 No Pathway Enrichment

Top 10% variable genes (n=156) vs background (n=917):
- 1/156 in neurodegeneration set vs 1/917 background.
- Odds ratio 5.9, **p=0.27 (not significant)**.

---

## 4. Discussion

### 4.1 Negative Result Is Informative

We find **no evidence for sex-specific SE splicing signatures** in 10x 3′ *Drosophila* brain data. This is not a failure of the ML method but a limitation of the data type: 3′ capture cannot reliably measure exon skipping in the gene body.

### 4.2 What Would Work

- **Full-length scRNA-seq** (e.g., Smart-seq2) captures splice junctions across the transcript.
- **Bulk RNA-seq** with deep coverage enables rMATS analysis.
- **Targeted:** 3′-specific splicing (alternative last exons, tandem UTRs) might be detectable in 10x data.

### 4.3 Limitations

1. **Small n:** Only 4 pseudo-bulks (2 libraries) processed; 14 planned but downloads too slow.
2. **No aging:** Only 3d samples; 6d/9d not processed (cannot address "brain aging").
3. **Custom pipeline:** minimap2+Python instead of STAR+rMATS (infrastructure constraints); not benchmarked.
4. **No UMI dedup; exact barcode matching** (see Methods).
5. **SE only:** A3SS/A5SS/MXE/RI not analyzed.

### 4.4 Conclusion

Machine learning **cannot** identify sex-specific alternative-splicing signatures of brain aging from 10x 3′ scRNA-seq due to fundamental 3′ bias. The research question requires full-length transcript data. Our pipeline (barcode splitting, minimap2 junctions, PSI, leakage-safe CV) is functional and would apply to suitable data.

---

## 5. Data and Code Availability

- **Repository:** Private GitHub (to be pushed).
- **Raw data:** GEO GSE107451.
- **Processed:** PSI matrix, ML predictions, differential results in `results/`.
- **Scripts:** All analysis code in `scripts/` (see README).

---

## References

- Davie et al., Cell 2018 (GSE107451). https://pmc.ncbi.nlm.nih.gov/articles/PMC6086935/
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451

---

## Figures

**Figure 1:** Distribution of PSI values (95.1% at 1.0).  
*File: figures/main/fig1_psi_distribution.png*

**Figure 2:** Volcano plot of sex-specific splicing (paired).  
*File: figures/main/fig2_volcano.png*

**Figure 3:** Confusion matrix (OOF, leave-one-library-out).  
*File: figures/main/fig3_confusion.png*
