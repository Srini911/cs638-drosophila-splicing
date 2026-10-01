#!/bin/bash
# Phase 9 — rMATS alternative splicing: male vs female for one age.
# Usage: run_rmats.sh AGE LIBTYPE
#   AGE: 3, 6, or 9 (days)
#   LIBTYPE: fr-unstranded | fr-firststrand | fr-secondstrand
#            (determined empirically from BAM strandedness; see logs/strandedness.txt)
# Inputs: data/intermediate/pseudobulk/w1118_<age>d_r{1,2,3}_{male,female}.bam
# Outputs: results/rmats/<age>d/
set -euo pipefail
export PATH="$HOME/workspace/tools/gsl/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/workspace/tools/gsl/lib:${LD_LIBRARY_PATH:-}"
AGE="$1"; LIBTYPE="${2:-fr-unstranded}"
BASE="$HOME/workspace/ML_Ozgun"
PB="$BASE/data/intermediate/pseudobulk"
OUT="$BASE/results/rmats/${AGE}d"
TMP="$BASE/data/intermediate/rmats_tmp/${AGE}d"
RMATS_PY="$HOME/workspace/tools/rmats-turbo/rmats.py"
GTF="$BASE/data/reference/dmel_BDGP6.54.63.gtf.gz"
READLEN=120

mkdir -p "$OUT" "$TMP"
# zcat GTF to a plain file for rMATS (it needs an uncompressed GTF path)
if [ ! -f "$BASE/data/reference/dmel_BDGP6.54.63.gtf" ]; then
  zcat "$GTF" > "$BASE/data/reference/dmel_BDGP6.54.63.gtf"
fi

B1=$(ls "$PB"/w1118_${AGE}d_r{1,2,3}_male.bam | paste -sd,)
B2=$(ls "$PB"/w1118_${AGE}d_r{1,2,3}_female.bam | paste -sd,)
echo "$B1" > "$TMP/b1.txt"; echo "$B2" > "$TMP/b2.txt"
echo "male BAMs: $B1"
echo "female BAMs: $B2"

python3 "$RMATS_PY" \
  --b1 "$TMP/b1.txt" --b2 "$TMP/b2.txt" \
  --gtf "$BASE/data/reference/dmel_BDGP6.54.63.gtf" \
  --od "$OUT" --tmp "$TMP" \
  -t paired --readLength "$READLEN" \
  --libType "$LIBTYPE" \
  --nthread 2 --tstat 2 \
  > "$BASE/logs/rmats_${AGE}d.log" 2>&1
echo "rMATS ${AGE}d complete -> $OUT"
