#!/bin/bash
# Master driver — per 10x library: download sralite -> fasterq-dump -> STARsolo
# -> split BAM by sex -> cleanup. Processes libraries SEQUENTIALLY to bound disk.
# Usage: pipeline_align_all.sh
# Logs: logs/pipeline_<lib>.log ; status: logs/pipeline_status.txt
set -uo pipefail
BASE="$HOME/workspace/ML_Ozgun"
SRA="$BASE/data/raw/sralite"
FQROOT="$BASE/data/raw/sra"
PB="$BASE/data/intermediate/pseudobulk"
LOGDIR="$BASE/logs"
mkdir -p "$PB" "$LOGDIR"
export PATH="$HOME/workspace/tools/sratoolkit.3.1.1-ubuntu64/bin:$PATH"

# lib -> SRR mapping (w1118 3d/6d/9d)
declare -A LIBS=(
  [w1118_3d_r1]=SRR6327125 [w1118_3d_r2]=SRR6327126 [w1118_3d_r3]=SRR6327127
  [w1118_6d_r1]=SRR6327128 [w1118_6d_r2]=SRR6327129 [w1118_6d_r3]=SRR6327130
  [w1118_9d_r1]=SRR6327131 [w1118_9d_r2]=SRR6327132 [w1118_9d_r3]=SRR6327133
)
ORDER=(w1118_3d_r1 w1118_3d_r2 w1118_3d_r3 w1118_6d_r1 w1118_6d_r2 w1118_6d_r3 w1118_9d_r1 w1118_9d_r2 w1118_9d_r3)

get_url() { grep "^$1 " /tmp/sralite_urls.txt | cut -d' ' -f2; }

process_lib() {
  local lib="$1" srr="$2"
  local log="$LOGDIR/pipeline_${lib}.log"
  echo "===== $lib ($srr) $(date) =====" | tee -a "$log"
  local fqdir="$FQROOT/$srr"; mkdir -p "$fqdir"

  # 1. download sralite (resume-capable)
  if [ ! -f "$SRA/${srr}.sralite.done" ]; then
    echo "[$lib] downloading sralite..." | tee -a "$log"
    curl -sL --retry 5 --retry-delay 15 -C - --max-time 10800 \
      -o "$SRA/${srr}.sralite" "$(get_url $srr)" >>"$log" 2>&1 \
      || { echo "[$lib] DOWNLOAD FAILED" | tee -a "$log"; return 1; }
    echo done > "$SRA/${srr}.sralite.done"
  fi

  # 2. fasterq-dump with technical reads (barcodes)
  if [ ! -f "$fqdir/fqdump.done" ]; then
    echo "[$lib] fasterq-dump..." | tee -a "$log"
    ln -sf "$SRA/${srr}.sralite" "$SRA/${srr}.sra"
    fasterq-dump "$SRA/${srr}.sra" --include-technical --split-files \
      -O "$fqdir" -e 2 -t "$fqdir/tmp" >>"$log" 2>&1 \
      || { echo "[$lib] FASTERQ-DUMP FAILED" | tee -a "$log"; return 1; }
    echo done > "$fqdir/fqdump.done"
    ls -la "$fqdir"/*.fastq >>"$log" 2>&1
  fi

  # 3. STARsolo
  if [ ! -f "$BASE/data/intermediate/starsolo/$lib/Aligned.sortedByCoord.out.bam" ]; then
    echo "[$lib] STARsolo..." | tee -a "$log"
    bash "$BASE/scripts/05_alignment/run_starsolo.sh" "$srr" "$lib" >>"$log" 2>&1 \
      || { echo "[$lib] STARSOLO FAILED" | tee -a "$log"; return 1; }
  fi

  # 4. split by sex
  if [ ! -f "$PB/${lib}_male.bam" ]; then
    echo "[$lib] splitting by sex..." | tee -a "$log"
    python3 "$BASE/scripts/07_rmats/split_bam_by_sex.py" \
      "$BASE/data/intermediate/starsolo/$lib/Aligned.sortedByCoord.out.bam" \
      "$lib" "$PB" >>"$log" 2>&1 \
      || { echo "[$lib] SPLIT FAILED" | tee -a "$log"; return 1; }
  fi

  # 5. cleanup heavy intermediates (keep split BAMs)
  echo "[$lib] cleaning up..." | tee -a "$log"
  rm -f "$fqdir"/*.fastq "$SRA/${srr}.sralite" "$SRA/${srr}.sra"
  rm -rf "$fqdir/tmp" "$BASE/data/intermediate/starsolo/$lib/Aligned.sortedByCoord.out.bam"
  echo "[$lib] COMPLETE $(date)" | tee -a "$log"
  echo "$lib COMPLETE" >> "$LOGDIR/pipeline_status.txt"
  return 0
}

for lib in "${ORDER[@]}"; do
  if grep -q "^$lib COMPLETE" "$LOGDIR/pipeline_status.txt" 2>/dev/null; then
    echo "[$lib] already complete, skipping"
    continue
  fi
  process_lib "$lib" "${LIBS[$lib]}" || { echo "PIPELINE STOPPED at $lib"; exit 1; }
done
echo "ALL LIBRARIES COMPLETE $(date)" | tee -a "$LOGDIR/pipeline_status.txt"
