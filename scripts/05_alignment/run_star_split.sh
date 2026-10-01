#!/bin/bash
# Phase 6b — Regular STAR alignment of barcode-split FASTQs.
# Usage: run_star_split.sh SRR LIB_LABEL SEX
#   SRR       : SRA run accession
#   LIB_LABEL : e.g. w1118_3d_r1
#   SEX       : male or female
# Input:  data/raw/sra/<SRR>/split/<SEX>_2.fastq (cDNA; _1 is barcodes, not aligned)
# Output: data/intermediate/star/<LIB_LABEL>_<SEX>/Aligned.sortedByCoord.out.bam
# Note: FASTQs are pre-split by cell barcode+sex, so no CB/UMI handling needed.
#       UMI deduplication is skipped (documented limitation); PSI is a ratio
#       and approximately robust to PCR duplicates.
set -euo pipefail
SRR="$1"; LIB="$2"; SEX="$3"
BASE="$HOME/workspace/ML_Ozgun"
FQDIR="$BASE/data/raw/sra/$SRR/split"
OUT="$BASE/data/intermediate/star/${LIB}_${SEX}"
INDEX="$BASE/data/reference/star_index"
GTF="$BASE/data/reference/dmel_BDGP6.54.63.gtf"

mkdir -p "$OUT"
R2="$FQDIR/${SEX}_2.fastq"
if [ ! -f "$R2" ]; then echo "missing $R2"; exit 1; fi

echo "[$LIB $SEX] aligning $(wc -l < "$R2" | awk '{print $1/4}') reads"

~/workspace/tools/STAR --runThreadN 2 \
  --genomeDir "$INDEX" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 119 \
  --readFilesIn "$R2" \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 2000000000 \
  --outFileNamePrefix "$OUT/" \
  > "$OUT/star.log" 2>&1

echo "[$LIB $SEX] done"
grep -E "Number of input reads|Uniquely mapped reads %" "$OUT/Log.final.out" || true
