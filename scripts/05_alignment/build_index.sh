#!/bin/bash
# Phase 5 — build STAR genome index for Drosophila BDGP6.54 (Release 6 assembly).
# sjdbOverhang 119 (10x cDNA read length 120, verified from FASTQ).
# genomeSAindexNbases 12 (required for small 137MB genome; 14 causes segfault).
set -euo pipefail
BASE="$HOME/workspace/ML_Ozgun"
REF="$BASE/data/reference"
INDEX="$REF/star_index"
mkdir -p "$INDEX"

STAR --runMode genomeGenerate --runThreadN 2 \
  --genomeDir "$INDEX" \
  --genomeFastaFiles "$REF/dmel_BDGP6.54.dna.primary_assembly.fa" \
  --sjdbGTFfile "$REF/dmel_BDGP6.54.63.gtf" \
  --sjdbOverhang 119 \
  --genomeSAindexNbases 12 \
  > "$REF/star_index_build.log" 2>&1
echo "index built: $INDEX"
