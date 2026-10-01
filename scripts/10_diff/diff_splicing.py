#!/usr/bin/env python3
"""
Differential splicing: male vs female, paired by library.
Uses PSI matrix; computes delta PSI and paired t-test.
Usage: diff_splicing.py
"""
import csv
import numpy as np
from pathlib import Path
from scipy import stats

def main():
    base = Path.home() / "workspace" / "ML_Ozgun"
    infile = base / "data" / "intermediate" / "psi_matrix.tsv"
    outdir = base / "results" / "diff_splicing"
    outdir.mkdir(parents=True, exist_ok=True)

    # Read PSI matrix
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        samples = [c for c in reader.fieldnames if c != "event_id"]
        events = []
        psi = {}
        for row in reader:
            eid = row["event_id"]
            events.append(eid)
            psi[eid] = {s: (float(row[s]) if row[s] != "NA" else np.nan) for s in samples}

    print(f"Samples: {samples}", flush=True)

    # Pair by library: w1118_3d_r1_female vs w1118_3d_r1_male, etc.
    # Extract library ID: w1118_3d_r1
    libs = sorted(set("_".join(s.split("_")[:3]) for s in samples))
    print(f"Libraries: {libs}", flush=True)

    results = []
    for eid in events:
        # Get paired values
        male_vals = []
        female_vals = []
        for lib in libs:
            m_sample = f"{lib}_male"
            f_sample = f"{lib}_female"
            if m_sample in samples and f_sample in samples:
                m = psi[eid][m_sample]
                f = psi[eid][f_sample]
                if not (np.isnan(m) or np.isnan(f)):
                    male_vals.append(m)
                    female_vals.append(f)

        # Require at least 2 pairs
        if len(male_vals) >= 2:
            male_vals = np.array(male_vals)
            female_vals = np.array(female_vals)
            delta = np.mean(male_vals - female_vals)

            # Paired t-test (if n>=2)
            if len(male_vals) >= 2:
                t, p = stats.ttest_rel(male_vals, female_vals)
            else:
                t, p = np.nan, np.nan

            # Mean PSI per sex
            mean_male = np.mean(male_vals)
            mean_female = np.mean(female_vals)

            results.append({
                "event_id": eid,
                "n_pairs": len(male_vals),
                "mean_psi_male": mean_male,
                "mean_psi_female": mean_female,
                "delta_psi": delta,
                "t_stat": t,
                "p_value": p,
            })

    # Sort by |delta_psi|
    results.sort(key=lambda x: abs(x["delta_psi"]), reverse=True)

    # Write
    outfile = outdir / "diff_splicing_paired.csv"
    with open(outfile, "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=["event_id", "n_pairs", "mean_psi_male",
                                          "mean_psi_female", "delta_psi", "t_stat", "p_value"])
        w.writeheader()
        w.writerows(results)

    print(f"Wrote {outfile}: {len(results)} events tested", flush=True)

    # Summary
    sig = [r for r in results if not np.isnan(r["p_value"]) and r["p_value"] < 0.05]
    print(f"Nominally significant (p<0.05): {len(sig)} / {len(results)}", flush=True)

    # Top 10 by |delta|
    print("\nTop 10 by |delta PSI|:", flush=True)
    for r in results[:10]:
        print(f"  {r['event_id'][:60]}: delta={r['delta_psi']:.3f}, p={r['p_value']:.3g}", flush=True)

if __name__ == "__main__":
    main()
