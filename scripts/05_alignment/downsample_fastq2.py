#!/usr/bin/env python3
"""
Memory-efficient downsample: two-pass, keeps random subset.
Usage: downsample_fastq2.py input.fastq output.fastq N [seed]
"""
import sys
import random

def main():
    infile, outfile, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42

    # Pass 1: count reads
    print("Counting reads...", flush=True)
    total = 0
    with open(infile) as f:
        for _ in f:
            total += 1
    total //= 4
    print(f"Total: {total:,} reads", flush=True)

    if total <= n:
        print("No downsampling needed", flush=True)
        import shutil
        shutil.copy(infile, outfile)
        return

    # Generate sorted random indices
    random.seed(seed)
    keep = set(random.sample(range(total), n))
    print(f"Selected {len(keep):,} indices", flush=True)

    # Pass 2: stream and keep
    kept = 0
    with open(infile) as fin, open(outfile, "w") as fout:
        idx = 0
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline(); p = fin.readline(); q = fin.readline()
            if idx in keep:
                fout.write(h); fout.write(s); fout.write(p); fout.write(q)
                kept += 1
            idx += 1
            if idx % 1000000 == 0:
                print(f"  {idx:,} / {total:,}...", flush=True)

    print(f"Done: kept {kept:,} reads", flush=True)

if __name__ == "__main__":
    main()
