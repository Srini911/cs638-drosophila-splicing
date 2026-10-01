#!/bin/bash
# Robust ENA FASTQ downloader for one SRR accession.
# Resumes partial downloads, retries transient failures, verifies byte counts
# against the ENA portal API before declaring success.
# Usage: fetch_ena.sh SRRxxxxxxx OUTDIR
set -u
SRR="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

# expected sizes from ENA API: "bytes1;bytes2"
META=$(curl -s --max-time 60 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp,fastq_bytes&format=json")
B1=$(echo "$META" | python3 -c "import json,sys; d=json.load(sys.stdin)[0]; print(d['fastq_bytes'].split(';')[0])")
B2=$(echo "$META" | python3 -c "import json,sys; d=json.load(sys.stdin)[0]; print(d['fastq_bytes'].split(';')[1])")
P1="https://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/00${SRR: -1}/${SRR}/${SRR}_1.fastq.gz"
# NOTE: ENA path layout is vol1/fastq/SRR632/005/SRR6327125/ ; derive generically:
PFX1=$(echo "$META" | python3 -c "import json,sys; d=json.load(sys.stdin)[0]; print(d['fastq_ftp'].split(';')[0])")
P1="https://${PFX1}"
PFX2=$(echo "$META" | python3 -c "import json,sys; d=json.load(sys.stdin)[0]; print(d['fastq_ftp'].split(';')[1])")
P2="https://${PFX2}"

dl_one() { # url outfile expected_bytes
  local url="$1" out="$2" exp="$3" tries=0
  while [ $tries -lt 12 ]; do
    curl -sL --retry 5 --retry-delay 10 -C - --max-time 3600 -o "$out" "$url"
    local got=$(stat -c%s "$out" 2>/dev/null || echo 0)
    if [ "$got" -eq "$exp" ]; then echo "OK $out ($got bytes)"; return 0; fi
    echo "partial $out: $got/$exp bytes, retrying... (attempt $((tries+1)))"
    tries=$((tries+1)); sleep 15
  done
  echo "FAILED $out after $tries attempts" 1>&2; return 1
}

dl_one "$P1" "${OUTDIR}/${SRR}_1.fastq.gz" "$B1" || exit 1
dl_one "$P2" "${OUTDIR}/${SRR}_2.fastq.gz" "$B2" || exit 1
echo "DONE ${SRR}" > "${OUTDIR}/done.txt"
