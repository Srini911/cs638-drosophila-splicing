#!/usr/bin/env python3
"""
Downsample FASTQ to N reads (reservoir sampling).
Usage: downsample_fastq.py input.fastq output.fastq N [seed]
"""
import sys
import random

def main():
    infile, outfile, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42
    random.seed(seed)

    reservoir = []
    total = 0
    with open(infile) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline(); p = f.readline(); q = f.readline()
            total += 1
            if len(reservoir) < n:
                reservoir.append((h, s, p, q))
            else:
                j = random.randint(0, total - 1)
                if j < n:
                    reservoir[j] = (h, s, p, q)

    with open(outfile, "w") as o:
        for h, s, p, q in reservoir:
            o.write(h); o.write(s); o.write(p); o.write(q)

    print(f"Downsampled {total:,} -> {len(reservoir):,} reads (seed={seed})", flush=True)

if __name__ == "__main__":
    main()
