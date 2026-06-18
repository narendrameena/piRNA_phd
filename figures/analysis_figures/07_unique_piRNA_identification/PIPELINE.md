# 07 — Strain- & timepoint-specific (UNIQUE) piRNA identification

**What these figures are.** Identification/classification of piRNAs private to one strain: length QC, the
strain-specific caller, the genuinely-unique-vs-expressed-elsewhere split, the SNP-variant route, and PCA.
Pilot C57BL/6NJ, CAST/EiJ, SPRET/EiJ × {E16.5,P12.5,P20.5}; scaled to 16 strains where noted (`*16`).
Re-traced from `unique_pirna/METHODS_review_notes.md`; the **four-route framework is BioMNI triple-verified**.

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1 · piRNA sequence layer.** cutadapt 5.0 → STAR 2.7.11b (unmasked strain genome, piRNA params) → collapsed
distinct sequences per strain × tp × rep.

**S2 · length window (data-driven).** distinct-sequence length distribution; window = FWHM of the bulk.
**Result #:** **24–32 nt** (FWHM 26–30) captures **96.7 %** of the unique set.

**S3 · strain-specific caller.** present **≥2/3 replicates** in strain AND absent (<2/3) in every other strain,
**∩ edgeR (R 4.2.3) quasi-likelihood DA** (`filterByExpr`, BH-FDR<0.05, logFC>0). **Result #:** strain-specific
candidates E16.5 **305 k** / P12.5 256 k / P20.5 274 k (≈**835 k** total).

**S4 · SEQUENCE level — genuinely-unique vs expressed-elsewhere (STAR genome-anchored).** re-test each
candidate **sequence** for expression in the OTHER 15 strains' expressed pools, anchored by **STAR 2.7.11b**
mapping to genome with the pipeline's piRNA params, **only** relaxing `--outFilterMismatchNmax 0→3` (catches
SNP-variants). expressed-exact / SNP-variant = **NOT unique** (same aligner as the rest of the pipeline =
tool-consistent). **Result #:** genuinely unique C57BL/6NJ **109 810** (58.5 %), CAST **142 751** (61.8 %),
SPRET **202 177** (61.9 %).

**S5 · mismatch cutoff (SNP-variant).** data-driven from the mismatch distribution (no magic number).

**S6 · LOCUS level — orthology by PANGENOME (cactus halLiftover), NOT STAR.** for the genuinely-unique
candidates, project the locus GRCm39→each strain through the cactus HAL (`lift_presence16.sh` →
`classify_unique16_locus.py`): orthologous locus **exists elsewhere** → *conserved-but-silent* (expression
change); locus in **no other strain** → *strain-private locus* (locus gain). **This resolves the
divergence-vs-absence confound that sequence mapping alone cannot** (a present-but-diverged locus is NOT
mis-called a new locus) — the same lesson as theme 03. STAR answers the **sequence** question; the pangenome
answers the **locus** question.

**S7 · PCA.** **DESeq2 (R 4.2.3)** size-factor norm → PCA (top-500 variable). **Result #:** PC1 variance rises
**63 → 75 → 87 %** (E16.5→P12.5→P20.5); SPRET most distinct.

**S8 · figures.** matplotlib (Python 3.11.15). → length, counts, step-4, SNP-variant, PCA, concept (18 figs).

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| cutadapt | 5.0 | trim sRNA | 20–36 nt, `--discard-untrimmed` |
| STAR | 2.7.11b | **sequence-level** genome-anchored uniqueness test | piRNA params; Step-4 mismatch 0→3 |
| cactus / halLiftover | v2.9.3 | **locus-level** cross-strain orthology (conserved-but-silent vs strain-private) | `halLiftover <HAL> GRCm39 <bed> <strain>` |
| edgeR | R 4.2.3 | quasi-likelihood DA (data-driven floor) | `filterByExpr`, QL, BH-FDR<0.05, logFC>0 |
| DESeq2 | R 4.2.3 | PCA normalisation | size factors; top-500 variable |
| Python | 3.11.15 | figures | matplotlib |

## INPUTS  `edger16/{tp}.*`, `unique16/*.step4_classified.csv.gz`, `pca16/*.pca.csv` → `data/SourceData_*.csv`.
## OUTPUTS (`figures/`, 18)  length/window, pilot/strain-specific counts, step4 classification, SNP-variant, PCA, four-routes concept.

## DOUBLE-VERIFICATION
- Thresholds data-driven (FWHM 24–32; edgeR FDR/logFC; mismatch cutoff) — no magic numbers.
- Four-route classification **BioMNI triple-verified** (METHODS §Biological grounding).
- Uniqueness rechecked (pool-level 100 % true-unique; raw-per-rep caveat in `project_locus_figure_redesign`).
