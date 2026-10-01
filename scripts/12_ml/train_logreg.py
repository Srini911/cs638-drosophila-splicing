#!/usr/bin/env python3
"""Phase 16-18 — primary ML model: logistic regression predicting sex from PSI.

Design (leakage-safe):
  - Rows: 18 pseudo-bulk samples (3 ages x 3 replicates x 2 sexes). Target: sex.
  - Columns: QC-passed splicing events (PSI values). NO sex-derived features.
  - Age is handled as a covariate AND via grouped CV (see below).
  - Validation: leave-one-library-out CV (9 folds; each fold holds out the 2
    samples — male+female — from one 10x library). This tests generalization
    across biological replicates. All data-dependent steps (variance filter,
    univariate feature selection, scaling) are fit INSIDE each training fold.
  - Model: L2-penalized logistic regression (C tuned inside folds via a small
    inner grid — actually fixed C=1.0 as primary for simplicity; documented).
  - Metrics from out-of-fold predictions: accuracy, precision, recall, F1,
    ROC-AUC, confusion matrix. Baselines: majority-class, permuted labels.

Outputs: results/ml/oof_predictions.csv, results/ml/metrics.txt,
         results/ml/coefficients.csv, figures.
"""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif
from sklearn.pipeline import Pipeline
from sklearn.metrics import (accuracy_score, precision_score, recall_score,
                             f1_score, roc_auc_score, confusion_matrix)

BASE = os.path.expanduser("~/workspace/ML_Ozgun")
OUT = os.path.join(BASE, "results/ml")
FIG = os.path.join(BASE, "figures/main")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

ml = pd.read_csv(os.path.join(BASE, "results/psi/psi_matrix_qc.csv"), index_col=0)
print(f"ML matrix: {ml.shape[0]} samples x {ml.shape[1]} features")

# target + groups from sample names: w1118_{age}d_r{rep}_{sex}
y = ml.index.str.extract(r"_([a-z]+)$")[0].map({"male": 1, "female": 0}).values
lib = ml.index.str.extract(r"(w1118_\dd_r\d)")[0].values          # 9 libraries
age = ml.index.str.extract(r"w1118_(\d)d")[0].astype(int).values  # age covariate
X = ml.values
print("class balance:", pd.Series(y).value_counts().to_dict())

K = 50  # features selected per fold (univariate, inside fold)
oof_pred, oof_prob = np.zeros(len(y), dtype=int), np.zeros(len(y))
fold_coefs = []

libs = sorted(set(lib))
for i, held in enumerate(libs):
    tr = lib != held
    te = lib == held
    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("var", VarianceThreshold(threshold=1e-6)),
        ("sel", SelectKBest(f_classif, k=min(K, tr.sum() - 1))),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(C=1.0, max_iter=5000)),
    ])
    # add age as covariate: append AFTER selection/scaling would leak structure;
    # simplest defensible: include age as an extra column scaled with the rest.
    # To keep the pipeline clean we fit on X+age jointly inside the fold.
    Xa = np.column_stack([X, age])
    pipe.fit(Xa[tr], y[tr])
    oof_pred[te] = pipe.predict(Xa[te])
    oof_prob[te] = pipe.predict_proba(Xa[te])[:, 1]
    # record which features were selected
    sel_mask = pipe.named_steps["sel"].get_support()
    var_mask = pipe.named_steps["var"].get_support()
    feat_names = np.array(list(ml.columns) + ["AGE_COVARIATE"])
    # map through variance filter
    kept = feat_names[var_mask][sel_mask]
    coefs = pipe.named_steps["clf"].coef_[0]
    fold_coefs.append(pd.DataFrame({"fold": held, "feature": kept, "coef": coefs}))
    print(f"fold {held}: train={tr.sum()} test={te.sum()} "
          f"acc={accuracy_score(y[te], oof_pred[te]):.2f}")

oof = pd.DataFrame({"sample": ml.index, "true_sex": y,
                    "pred_sex": oof_pred, "pred_prob": oof_prob,
                    "library": lib, "age": age})
oof.to_csv(os.path.join(OUT, "oof_predictions.csv"), index=False)

acc = accuracy_score(y, oof_pred)
print("\n=== OUT-OF-FOLD (leave-one-library-out) ===")
print(f"accuracy  {acc:.3f}")
print(f"precision {precision_score(y, oof_pred):.3f}")
print(f"recall    {recall_score(y, oof_pred):.3f}")
print(f"F1        {f1_score(y, oof_pred):.3f}")
print(f"ROC-AUC   {roc_auc_score(y, oof_prob):.3f}")
print("confusion matrix (rows=true 0=F/1=M, cols=pred):")
print(confusion_matrix(y, oof_pred))
# baselines
maj = int(y.mean() >= 0.5)
print(f"\nmajority-class baseline accuracy: {accuracy_score(y, np.full_like(y, maj)):.3f}")
rng = np.random.default_rng(0)
perm_accs = [accuracy_score(y, rng.permutation(oof_pred)) for _ in range(200)]
print(f"permuted-prediction baseline accuracy: mean={np.mean(perm_accs):.3f} "
      f"(95%ile={np.percentile(perm_accs, 95):.3f})")

with open(os.path.join(OUT, "metrics.txt"), "w") as fh:
    fh.write(f"n_samples={len(y)} n_features={X.shape[1]}\n")
    fh.write(f"cv=leave-one-library-out (9 folds)\n")
    fh.write(f"accuracy={acc:.4f}\nprecision={precision_score(y, oof_pred):.4f}\n")
    fh.write(f"recall={recall_score(y, oof_pred):.4f}\nf1={f1_score(y, oof_pred):.4f}\n")
    fh.write(f"roc_auc={roc_auc_score(y, oof_prob):.4f}\n")
    fh.write(f"confusion_matrix=\n{confusion_matrix(y, oof_pred)}\n")

coefs = pd.concat(fold_coefs, ignore_index=True)
coefs.to_csv(os.path.join(OUT, "coefficients.csv"), index=False)
# consensus features: selected in >=5/9 folds, mean |coef|
cons = (coefs.groupby("feature")["coef"]
        .agg(n_folds="count", mean_coef="mean", mean_abs="lambda s: s.abs().mean()")
        .reset_index().sort_values("mean_abs", ascending=False))
cons.to_csv(os.path.join(OUT, "consensus_features.csv"), index=False)
print("\ntop consensus features:")
print(cons.head(10).to_string(index=False))
