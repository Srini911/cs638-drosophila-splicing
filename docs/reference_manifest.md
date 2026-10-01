# Reference Manifest — Phase 5

**Date:** 2026-10-01

## Reference used

| Item | Value |
|---|---|
| Assembly | BDGP Release 6 (Drosophila melanogaster) |
| Annotation | Ensembl BDGP6.54, genebuild 2023-09 (release 116) |
| FASTA | `data/reference/dmel_BDGP6.54.dna.primary_assembly.fa.gz` (8 chromosomes: 2L, 2R, 3L, 3R, 4, X, Y, mitochondrion) |
| GTF | `data/reference/dmel_BDGP6.54.63.gtf.gz` |
| FASTA md5 | 897605ff8cb67a23424dc26ae1aad46f |
| GTF md5 | 4aaa5a1bb9d44253d93d23398c2b8afb |
| Downloaded | 2026-10-01 from https://ftp.ensembl.org/pub/release-116/ |
| Gene IDs | FBgn (FlyBase-sourced gene models; roX2 = FBgn0019660 confirmed present) |

## Documented substitution

Prior course work used the **FlyBase r6.64 GTF**. That file could not be obtained in
this compute environment: `ftp.flybase.net` (FTP and HTTPS) is unreachable through
the environment's egress proxy (verified 2026-09-30: FTP returns empty reply,
HTTPS connection fails).

Substitution rationale (explicit, not silent):
- All FlyBase r6.xx annotations share the **BDGP Release 6 assembly** — genomic
  coordinates are identical between r6.64 and BDGP6.54. rMATS event coordinates are
  therefore directly comparable.
- Ensembl's *D. melanogaster* annotation is built from FlyBase gene models, so
  FBgn identifiers (used throughout this project, e.g. roX2) are preserved.
- FASTA and GTF come from the same Ensembl release → guaranteed compatibility
  (no mixed releases).

Consequence: gene-model differences between r6.64 and BDGP6.54 may cause small
differences in event counts vs. prior work. Any such disagreement will be
investigated per the brief (§10), not hidden.

## STAR index
- Built with STAR 2.7.11b: `scripts/05_alignment/build_index.sh`
- sjdbOverhang 119 (10x cDNA read length 120bp, verified from SRR6327125 FASTQ)
- Location: `data/reference/star_index/` (not committed to git; reproducible via script)
