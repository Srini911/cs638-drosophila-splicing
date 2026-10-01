#!/bin/bash
# Download all 9 sralite files sequentially (bandwidth-optimal: single stream).
# Usage: download_all_sralite.sh
# Reads URLs from /tmp/sralite_urls.txt ; writes *.done markers.
set -u
SRA="$HOME/workspace/ML_Ozgun/data/raw/sralite"
mkdir -p "$SRA"
while read -r srr url; do
  out="$SRA/${srr}.sralite"
  if [ -f "$SRA/${srr}.done" ]; then echo "[$srr] already done, skipping"; continue; fi
  echo "[$srr] downloading $(date)"
  for i in $(seq 1 10); do
    if curl -sL --retry 3 -C - --max-time 10800 -o "$out" "$url"; then
      echo done > "$SRA/${srr}.done"
      echo "[$srr] complete $(date)"
      break
    fi
    echo "[$srr] retry $i"; sleep 20
  done
done < /tmp/sralite_urls.txt
echo "ALL SRALITE DOWNLOADS COMPLETE $(date)"
