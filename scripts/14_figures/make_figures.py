#!/usr/bin/env python3
"""
Generate main figures for the report.
Usage: make_figures.py
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    base = Path.home() / "workspace" / "ML_Ozgun"
    figdir = base / "figures" / "main"
    figdir.mkdir(parents=True, exist_ok=True)

    # Fig 1: PSI distribution
    psis = []
    with open(base / "data" / "intermediate" / "psi_matrix.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        samples = [c for c in reader.fieldnames if c != "event_id"]
        for row in reader:
            for s in samples:
                if row[s] != "NA":
                    psis.append(float(row[s]))

    plt.figure(figsize=(8, 5))
    plt.hist(psis, bins=50, edgecolor="black", alpha=0.7)
    plt.xlabel("PSI (Percent Spliced In)")
    plt.ylabel("Count")
    plt.title(f"Distribution of PSI values (n={len(psis)}; {len(samples)} samples)")
    plt.axvline(1.0, color="red", linestyle="--", label="PSI=1.0 (constitutive)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figdir / "fig1_psi_distribution.png", dpi=150)
    plt.close()
    print("Fig1 done", flush=True)

    # Fig 2: Delta PSI volcano (paired)
    deltas = []
    pvals = []
    with open(base / "results" / "diff_splicing" / "diff_splicing_paired.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                d = float(row["delta_psi"])
                p = float(row["p_value"])
                if not np.isnan(p) and p > 0:
                    deltas.append(d)
                    pvals.append(p)
            except:
                pass

    deltas = np.array(deltas)
    pvals = np.array(pvals)
    neglogp = -np.log10(pvals)

    plt.figure(figsize=(8, 5))
    plt.scatter(deltas, neglogp, alpha=0.5, s=10)
    plt.xlabel("Delta PSI (male - female)")
    plt.ylabel("-log10(p-value)")
    plt.title("Volcano plot: sex-specific splicing (paired by library)")
    plt.axhline(-np.log10(0.05), color="red", linestyle="--", label="p=0.05")
    plt.axvline(0, color="gray", linestyle="-", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figdir / "fig2_volcano.png", dpi=150)
    plt.close()
    print("Fig2 done", flush=True)

    # Fig 3: ML OOF predictions (if available)
    ml_file = base / "results" / "ml" / "oof_predictions.csv"
    if ml_file.exists():
        import pandas as pd
        df = pd.read_csv(ml_file)
        plt.figure(figsize=(6, 6))
        # Confusion matrix
        from sklearn.metrics import confusion_matrix
        cm = confusion_matrix(df["true_sex"], df["pred_sex"])
        plt.imshow(cm, cmap="Blues")
        plt.colorbar()
        plt.xticks([0, 1], ["Female", "Male"])
        plt.yticks([0, 1], ["Female", "Male"])
        plt.xlabel("Predicted")
        plt.ylabel("True")
        plt.title("Confusion Matrix (OOF, leave-one-library-out)")
        for i in range(2):
            for j in range(2):
                plt.text(j, i, cm[i, j], ha="center", va="center", fontsize=14)
        plt.tight_layout()
        plt.savefig(figdir / "fig3_confusion.png", dpi=150)
        plt.close()
        print("Fig3 done", flush=True)

    print(f"Figures saved to {figdir}", flush=True)

if __name__ == "__main__":
    main()
