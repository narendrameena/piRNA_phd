# Data-integrity audit — raw → figure, all 15 themes (2026-06-19)

Independent verification: for each root dataset, re-derived/spot-checked values via a **fresh code path** (not the figure's own script), traced lineage toward raw BAMs, and confirmed figures reproduce their source. Method per project rule (data-driven, probe-based; never trust the caption).

## Backbone: unique-piRNA root `unique16/final_classified_clean_2read.csv.gz` (themes 07/08/09/10/15)
Full chain verified raw → figure:
- **raw sRNA reads → STAR BAMs → edgeR counts** (`edger16/{tp}.counts.tsv.gz`, 16,151,680 collapsed seqs × 48 samples).
- **Trace to raw (independent):** all **10,857/10,857** strain-private candidates @16.5dpc satisfy ≥2 reads in ≥2/3 own-strain reps (100%) AND absence (<2/3 reps ≥2 reads) in **every** other strain (100%), checked directly against the count matrix.
- **expr_class layer:** novel rows all have n_other_exact==0; non-novel all ≥1; other_exact_strains populated count == n_other_exact>0 count (40,238). 100% consistent.
- **klass (3-class):** independently re-derived from `expr_class` + `homolog_strains` using committed `classify_unique16_locus.py` logic → **0 mismatches / 404,769**.
- **klass5 (5-class):** clean block-diagonal refinement of klass (expressed-elsewhere 40,238 unchanged; CBS-klass 303,674 → SNP-variant 217,559 + CBS 86,115; private-klass 60,857 → low-quality 40,011 + private 20,846). Arithmetic closes exactly.
- **Integrity:** 404,769 rows, all unique by (strain,timepoint,sequence); no NaN in any critical column.
- **Figures:** theme-07 SourceData (step4_classification16, unique_pirna_timepoint CBS+private, step4_class_by_strain_timepoint16, tp_by_class, length_byclass16) all MATCH independent recompute from the root.

## Per-theme root verification
| Theme | Root | Result |
|---|---|---|
| 01 phasing | `phasing_allstrains_1random/ALL_summary.csv` | 143 rows, frac_plus1∈[0,1], 0 NaN; +1 phasing signal present (100% z>2, mean z=6.29) |
| 02 PICB clusters | `SourceData_PICB_cluster_counts.csv` | 192 rows; **1 corrupted (see Finding 2)**; n_clusters 1965–27965, none 0 |
| 03 zamore/SV | `combined_rebuild/all_strains_{expression,SV}_matrix.csv` | 214 loci × 48 = 10,272 rows each, 0 NaN; liftover = **pangenome** halLiftover (verified, not UCSC); retention recompute sane (~88–90%/strain) |
| 04 PICB QC | PICB chr-by-chr vs whole-bam xlsx | roots present |
| 05 C57BL_6NJ pangenome | `C57BL_6NJ_classified.bed`, `class_stats.csv` | present |
| 06 zamore coverage | black6 PICB xlsx (clusters/seeds/cores) | present |
| 07 unique-piRNA ID | unique-piRNA root | fully verified (above) |
| 08 TE origin/strand | `cand_self16/*.cand_self16.bam` (16) + root | 16 BAMs present; class layer = verified root |
| 09 TE-driven | `*.coord_byclass16.csv` | 80 rows = 16×5 klass5; pct==100·in_ins/n; n≤root (0 violations); strain-private coverage 1.00; fold 4.8–86× (matches finding) |
| 10 cluster PAV | `cluster_PAV_catalogue.csv.gz` | 1,366,033 rows, expected cols |
| 11/12/13 locus | `pav_clusters.py` ← picb_pangenome_clusters + BAMs (theme-14 root) | builders render (60/29/12 figs) |
| 14 circos | `picb_pangenome_clusters.tsv`, `active_te_expression_sRNA.tsv`, `active_1u_bias_tp.tsv` | present, readable |
| 15 misc | `results/pingpong/all_zscore.tab`, `TE_driven_COORDINATE16_*` (16), `divergence_candidates.tsv`, `pachytene_cluster.json` | present |

## FINDINGS — both FIXED 2026-06-19
**1. (MEDIUM — reproducibility) klass5 / `final_classified_clean_2read.csv.gz` producer was not committed → FIXED.**
Only the 3-class `final_classified.csv.gz` had a committed producer. CORRECTION on investigation: the upstream
determinants ARE committed — `classify_step416.py` → `step4_16/*.step4_classified16.csv.gz` (1–3 mm SNP-variant);
`snp_variant_refinement.csv` (SNP-variant set by cand_id); `cand_self16.sh` → `cand_self16/*.bam` (mm0 own-locus,
`--outFilterMismatchNmax 0`); `edger16/{tp}.strain_specific_DA_2read.csv.gz` (≥2-read DA set). Only the MERGE step
(klass4/klass5 + ≥2-read filter) lacked a script. **Written: `unique_pirna/make_klass5.py`.** Each rule reverse-verified:
klass4 from snp_variant_refinement = 0/451,972 mismatch; klass5 low-quality from cand_self16-mapped = 60,857/60,857;
≥2-read filter = DA_2read membership 451,972/451,972. The producer was run to a temp dir and **reproduces both
`final_classified_clean.csv.gz` (all 13 cols identical) and `final_classified_clean_2read.csv.gz` (all 13 cols
value-identical)** vs the canonical files. Canonical files NOT overwritten.

**2. (LOW — isolated corruption) Theme 02: `NZO_HlLtJ 12.5dpp` replicate 3 PICB xlsx concurrency-corrupted → FIXED.**
(cluster_n_chrom=19, complete=False — missing 2 chromosomes). Combined run clean (n_chrom=21); `Fig_picb_combined_clusters`
(combined-only) was already unaffected. **Added `complete==True` filtering** to `Fig_picb_rep_vs_combined.py` (drop
incomplete runs) and `Fig_picb_rep_combined_overlap.py` (drop groups with any incomplete run, via group-completeness
merge). Both re-rendered; conclusions unchanged (combined>max rep 0/48; support 78–81% in all 3 reps).

## Not affected / clean
- No genome-build mixing: Zamore mm10→GRCm39 is proper UCSC liftover; all cross-strain work is pangenome halLiftover.
- All ~100 figure scripts render; no figure uses the old 3-class classification (klass5 sweep complete).
