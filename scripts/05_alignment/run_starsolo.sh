#!/bin/bash
# Phase 6 — STARsolo alignment of one 10x v2 library.
# Usage: run_starsolo.sh SRRxxxxxxx LIB_LABEL
#   SRRxxxxxxx : SRA run accession (FASTQs already in data/raw/sra/<SRR>/)
#   LIB_LABEL  : e.g. w1118_3d_r1 (used for barcode whitelist + output naming)
# Produces: data/intermediate/starsolo/<LIB_LABEL>/Aligned.sortedByCoord.out.bam
set -euo pipefail
SRR="$1"; LIB="$2"
BASE="$HOME/workspace/ML_Ozgun"
FQDIR="$BASE/data/raw/sra/$SRR"
OUT="$BASE/data/intermediate/starsolo/$LIB"
WL="$BASE/metadata/barcodes/${LIB}.txt"   # 16bp barcodes for this library (from GEO metadata)
INDEX="$BASE/data/reference/star_index"

mkdir -p "$OUT"
# fasterq-dump writes uncompressed .fastq; fall back to .fastq.gz if present
if [ -f "$FQDIR/${SRR}_2.fastq.gz" ]; then
  R2="$FQDIR/${SRR}_2.fastq.gz"; R1="$FQDIR/${SRR}_1.fastq.gz"; ZCAT="--readFilesCommand zcat"
else
  R2="$FQDIR/${SRR}_2.fastq"; R1="$FQDIR/${SRR}_1.fastq"; ZCAT=""
fi
# read length of cDNA read
RLEN=$(head -2 "$R2" | tail -1 | wc -c); RLEN=$((RLEN-1))
echo "[$LIB] cDNA read length: $RLEN"

STAR --runThreadN 2 \
  --genomeDir "$INDEX" \
  --sjdbGTFfile "$BASE/data/reference/dmel_BDGP6.54.63.gtf" \
  --sjdbOverhang 119 \
  --readFilesIn "$R2" "$R1" \
  $ZCAT \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WL" \
  --soloCBmatchWLtype 1MM \
  --soloUMIlen 10 \
  --soloCellFilter None \
  --soloFeatures GeneFull \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes CB UB \
  --limitBAMsortRAM 2000000000 \
  --outFileNamePrefix "$OUT/" \
  > "$OUT/starsolo.log" 2>&1

echo "[$LIB] done. log: $OUT/starsolo.log"
grep -E "Number of input reads|Uniquely mapped reads %" "$OUT/Log.final.out" || true
