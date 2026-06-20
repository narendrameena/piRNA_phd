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

**S3 · strain-specific caller.** present **≥2/3 replicates** in strain AND absent in every other strain,
**∩ edgeR (R 4.2.3) quasi-likelihood DA** (`filterByExpr`, BH-FDR<0.05, logFC>0). **ABSENCE RULE — ADOPTED 2026-06-18: ≥2-read**
(an other strain counts as "having it" only at **≥2 reads total**; a single read = index-hopping/contamination noise — was the
looser "<2/3 replicates"; BioMNI 3/3 signed off). **Result # (16-strain `edger16` ∩ presence/absence):** 451,972 loose →
**404,769 ≥2-read** (`strain_specific_DA_2read.csv.gz`). The absence ladder (loose→≥2-read→strict) + read-mass coverage are in
`Fig_strain_specific_DA16_decomposition` / `Fig_strain_specific_coverage16`; ≥1-rep vs ≥2-rep presence is moot (filterByExpr forces all 3 reps for a private call).

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

**S8 · figures.** matplotlib (Python 3.11.15). → length, counts, step-4, SNP-variant, PCA, concept, + the 16-strain
strain-specific-DA / decomposition / coverage / class-by-strain×timepoint set (22 figs).

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
## OUTPUTS (`figures/`, 22)  length/window, pilot/strain-specific counts, **`Fig_strain_specific_DA16`** (+ replicate error bars), **`_decomposition`** (DA-only vs presence/absence-only vs intersection), **`Fig_strain_specific_coverage16`** (% of total piRNA read-mass), step4 classification (+ `class_by_strain_timepoint16`, `class_strain_timepoint16`), SNP-variant, PCA, four-routes concept.

## DOUBLE-VERIFICATION
- Thresholds data-driven (FWHM 24–32; edgeR FDR/logFC; mismatch cutoff) — no magic numbers.
- Four-route classification **BioMNI triple-verified** (METHODS §Biological grounding).
- Uniqueness rechecked (pool-level 100 % true-unique; raw-per-rep caveat in `project_locus_figure_redesign`).
- **≥2-read absence adoption (2026-06-18) is robustness-validated:** it removes ONLY "expressed-elsewhere (exact)" candidates (87,441→40,238); every genuinely-unique class is 100 % preserved (strain-private 20,846, conserved-but-silent 86,115, SNP-variant 217,559) → TE-origin / depth-confound / PCA figures are byte-identical and unchanged. BioMNI 3/3 signed off; full per-figure audit in the `Fig_strain_specific_DA16` note.

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# S3 edgeR QL-DA (>=2-read presence/absence intersection):
bash analysis/claude_biomni_analysis/unique_pirna/run_edger16.sh
# S4 sequence-level uniqueness (STAR mismatch 0->3) + step-4 classification:
bash analysis/claude_biomni_analysis/unique_pirna/cand_self16.sh
bash analysis/claude_biomni_analysis/unique_pirna/run_classify_step4.sh
# S6 locus orthology by pangenome (halLiftover) -> conserved-but-silent vs strain-private:
bash analysis/claude_biomni_analysis/unique_pirna/lift_presence16.sh
python analysis/claude_biomni_analysis/unique_pirna/classify_unique16_locus.py
# S7 PCA (DESeq2):
bash analysis/claude_biomni_analysis/unique_pirna/run_pca.sh
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_concept_four_routes.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_classes16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_numbers_summary.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_unique.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_unique16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pilot_unique_pirna_counts.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pirna_length_window.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_snp_variant_coord.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_snp_variant_nucleotide.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_spret_split.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_class_by_strain_timepoint16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_class_strain_timepoint16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_classification.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_classification16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_DA.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_DA16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_DA16_decomposition.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_coverage16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_timepoint_combos16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique16_class_breakdown.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_length.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_length16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_length_byclass16.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_timepoint.py
$PY figures/analysis_figures/07_unique_piRNA_identification/code/build_replicate_DA_detection.py
Rscript figures/analysis_figures/07_unique_piRNA_identification/code/combine_pca16_byclass.R
Rscript figures/analysis_figures/07_unique_piRNA_identification/code/combine_pca16_classes.R
Rscript figures/analysis_figures/07_unique_piRNA_identification/code/coverage_probe.R
$PY figures/analysis_figures/07_unique_piRNA_identification/code/depth_confound_check.py
Rscript figures/analysis_figures/07_unique_piRNA_identification/code/edger16_2read.R
Rscript figures/analysis_figures/07_unique_piRNA_identification/code/edger16_coverage.R
Rscript figures/analysis_figures/07_unique_piRNA_identification/code/presence_only_counts.R
Rscript figures/analysis_figures/07_unique_piRNA_identification/code/threshold_probe.R
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_concept_four_routes.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_classes16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_numbers_summary.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_unique.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pca_unique16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pilot_unique_pirna_counts.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_pirna_length_window.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_snp_variant_coord.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_snp_variant_nucleotide.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_spret_split.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_class_by_strain_timepoint16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_class_strain_timepoint16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_classification.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_step4_classification16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_DA.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_DA16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_DA16_decomposition.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_strain_specific_coverage16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_timepoint_combos16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique16_class_breakdown.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_length.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_length16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_length_byclass16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/Fig_unique_pirna_timepoint.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/build_replicate_DA_detection.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/combine_pca16_byclass.R`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/combine_pca16_classes.R`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/coverage_probe.R`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/depth_confound_check.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/edger16_2read.R`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/edger16_coverage.R`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/presence_only_counts.R`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/code/threshold_probe.R`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/run_edger16.sh` — S3 edgeR quasi-likelihood DA (16 strains)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cand_self16.sh` — S4 STAR genome-anchored uniqueness (mismatch 0->3)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/run_classify_step4.sh` — S4 step-4 expressed/SNP-variant/unique split driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/run_step4_map16.sh` — S4 step-4 mapping (16 strains)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/lift_presence16.sh` — S6 locus halLiftover GRCm39->each strain
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/classify_unique16_locus.py` — S6 conserved-but-silent vs strain-private locus classifier
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/run_pca.sh` — S7 DESeq2 PCA driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/run_classify_unique16.sh` — 16-strain unique classification driver
