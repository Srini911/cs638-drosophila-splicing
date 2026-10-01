# Project Requirements Audit — CS 638 Applied Machine Learning Semester Project

**Date:** 2026-09-30
**Status:** DRAFT — professor-provided materials were NOT found locally; requirements below are
reconstructed from the user's MASTER PROJECT INSTRUCTION (the authoritative spec available).
Items marked NEEDS USER INPUT must be confirmed by the student.

## Audit method

Searched the local VM for professor materials: project instructions, grading rubric, lecture
notes, syllabus, and deadline notices.

- Searched `~` and `/tmp` (depth 5) for `*638*`, `*ozgun*`, `*syllabus*`, `*rubric*`, `*.pdf`,
  `*.docx`, `*.pptx`.
- Result: **no CS 638 professor files found.** Only unrelated files (resume docs) exist.
- The Unity HPC paths referenced in the project brief
  (`/project/pi_mukulika_ray_umb_edu/ML_Ozgun`,
  `/project/pi_mukulika_ray_umb_edu/Aging_Dataset_GSE107451/Wild_Datasets_GSE107451`)
  **do not exist on this VM** — prior project outputs are not accessible here.

## Requirement table

| Requirement | Source | How we will satisfy it | Project artifact |
|---|---|---|---|
| Answer: "Can ML identify sex-specific alternative-splicing signatures of brain aging, enriched in neurodegeneration pathways?" (Drosophila brain) | Master instruction §1 | Full pipeline: rMATS → PSI matrix → leakage-safe ML → pathway enrichment | `report/final_report.*`, `results/` |
| Start from dataset discovery; independently justify GSE107451 | Master instruction §1, §4–5 | Phase 1 search log + candidate table; Phase 2 validation doc | `docs/dataset_search_log.md`, `results/dataset_candidates.csv`, `docs/GSE107451_dataset_validation.md` |
| Logistic regression as primary classifier unless course materials justify otherwise | Master instruction §2 | Logistic regression first; comparisons only if course-appropriate | `scripts/12_ml/`, `results/ml/` |
| No fabrication of data/results/citations/conclusions at any step | Master instruction §2 | Every claim traced to a file; claim-evidence table | `docs/claim_evidence_table.csv` |
| Leakage-safe validation (feature selection/scaling inside folds; no test influence) | Master instruction §2, §14–15, §18 | Group-aware CV (e.g., leave-one-group-out over replicates); pipelines fit on train folds only | `scripts/12_ml/`, `results/ml/` |
| Report training vs held-out performance separately; baselines + sanity checks | Master instruction §18–19 | Out-of-fold prediction table; majority-class/permutation baselines | `results/ml/` |
| Conservative biological interpretation; no causal/disease claims beyond evidence | Master instruction §22, §26 | Interpretation section distinguishes result vs literature vs hypothesis | `report/final_report.*` |
| Reproducibility: environment, versions, seeds, commands recorded | Master instruction §26, §31 | `environment.yml`, run logs, relative paths | `environment.yml`, `logs/` |
| Clean project directory; no raw/large data in git | Master instruction §27, §29 | Scaffold per §27; strong `.gitignore` | `ML_Ozgun/`, `.gitignore` |
| Professional README | Master instruction §28 | README with all listed sections | `README.md` |
| PRIVATE GitHub repo; no secrets/tokens; show upload plan before pushing | Master instruction §29 | Pre-push audit + file list review with user | remote repo (private) |
| Final report with honest methods/results/interpretation/limitations split | Master instruction §30, §36 | Report written only from actual outputs | `report/final_report.*` |
| Figure/table plan telling a coherent story | Master instruction §31, §37 | 8-figure sequence adapted to actual results | `figures/main/`, `figures/supplementary/` |
| Four audits before completion (scientific, ML, repo, report) | Master instruction §32 | Audit checklists executed and recorded | `docs/` audit notes |
| FINAL_DELIVERABLES.md + PROJECT_REPRODUCTION.md | Master instruction §33 | Generated at completion | `FINAL_DELIVERABLES.md`, `PROJECT_REPRODUCTION.md` |
| Work phase by phase; stop and diagnose on failure; never fake outputs | Master instruction §40 | Phase gates; validation before destructive/expensive steps | this workflow |
| Deadline for final submission | NEEDS USER INPUT | User to provide date | — |
| Grading rubric weights (e.g., model novelty vs rigor vs report) | NEEDS USER INPUT | User to provide rubric or course page | — |
| Required report sections/length, presentation requirement | NEEDS USER INPUT | User to provide; defaults to §30 structure | `report/` |
| Allowed/prohibited ML methods per lectures | NEEDS USER INPUT | User to provide; default: logistic regression + minimal course-typical comparisons | `scripts/12_ml/` |
| Whether Unity HPC access is available for STAR/rMATS compute | NEEDS USER INPUT | Determines heavy-compute path (see 2-day plan) | `scripts/05_alignment/`, `scripts/07_rmats/` (SLURM) |

## Decisions / defaults adopted in absence of professor materials

1. Primary model: logistic regression (per §2).
2. Validation: group-aware cross-validation with all data-dependent steps inside folds (per §14).
3. rMATS significance rule: FDR < 0.05 and |ΔPSI| ≥ 0.1 as an initial, documented rule —
   **to be re-confirmed against any professor threshold once materials are provided.**
4. No model zoo; no deep learning; no promoter analysis in the core pipeline (optional §29 only).
