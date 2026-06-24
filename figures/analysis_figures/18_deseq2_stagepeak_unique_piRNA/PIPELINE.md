# 18 — DESeq2 stage-peak (27/30 nt) WITHIN-TIMEPOINT unique piRNAs

**What this theme is.** A self-contained, **DESeq2-based** re-derivation of strain-specific ("unique") piRNAs,
restricted to the stage-characteristic piRNA length classes and judged **strictly within each developmental
stage**. It is **additive** — theme 07 (the edgeR-based unique-piRNA pipeline) and every other theme are left
untouched; this theme keeps the DESeq2 + stage-peak version alongside them. Built 2026-06-23 on user direction.

Three method decisions, each **data-driven and verified**:
1. **DESeq2, not edgeR**, for the differential-abundance (DA) engine — proven on our data + BioMNI 3/3.
2. **27/30 nt stage-peak length filter applied BEFORE DESeq2** — tested vs filter-after.
3. **Uniqueness judged WITHIN timepoint** (not pooled across stages) — BioMNI 3/3, biologically grounded.

---

## STEP-BY-STEP: raw data → figure

**S0 · inputs (reused, native).** Per-tp count matrices `unique_pirna/edger16/{tp}.{counts.tsv.gz,seqs.txt.gz,
samples.tsv}` (25–32 nt sRNA sequences × 48 libraries = 16 strains × 3 reps; counts header = sample order;
`libsize_window` = total 25–32 nt reads per library). All counts are integer, multimapper-weighted.

**S1 · benchmark edgeR vs DESeq2 (which DA method?).** `benchmark_da_methods.R` (SLURM `run_bench_da.sh`):
one-vs-rest per strain on a 300k-feature representative sample, 3 tp × 5 strain-label permutations. Metrics
(design verified by BioMNI 3/3; refs below): real-data concordance, **permutation-null false-positive control**,
**p-value calibration** (uniformity). **Result: DESeq2 is correctly calibrated; edgeR is ~2× anti-conservative
(null p<0.05 ≈ 0.09–0.10 vs target 0.05) and gives ~10–21× more null false-positives.** edgeR's ~20–25% extra
real calls are mostly false positives. → **DESeq2 adopted.** (`Fig_benchmark_da`)

**S2 · stage-peak length filter BEFORE DESeq2.** Per tp keep only the stage-characteristic peak length(s):
**E16.5 → 27 nt; P12.5 → 27 & 30 nt; P20.5 → 30 nt** (pre-pachytene ~27 / pachytene ~30; user-specified exact
modal lengths). `deseq_filter_order_test.R`: raw → length-filter → `filterByExpr` → DESeq2 one-vs-rest
(padj<0.05 & log2FC>0) → + ≥2-read presence/absence rule (present in focal X = ≥2 of 3 reps; absent elsewhere =
every other strain <2 reads total) → **stage-peak strain-specific candidates** (Order A = production).

**S3 · filter-order test (before vs after).** Same script also runs Order B (DESeq2 on the full 25–32 nt set,
then subset to 27/30). They differ only in size factors, dispersion trend and the BH denominator. **Result:
filtering BEFORE recovers ~6–7 % more stage-peak candidates** (less multiple-testing burden, length-focused
normalization), Jaccard 0.85–0.93. → **filter-before adopted.** (`Fig_filter_order`)

**S4 · candidates.** `deseq16_lenfilt/all_orderA_stagepeak_candidates.csv.gz` = **15,967** (E16.5 3,777 / P12.5
9,717 / P20.5 2,473). These are a **100 % subset of the edgeR candidate set at strain|tp|seq**, so every
classification determinant (cand_self mm0; step4 cand_to_Y SNP; Cactus/halLiftover locus-presence) is **reused —
no STAR/halLiftover recompute** (`classify_deseq_stagepeak.py`).

**S5 · WITHIN-TIMEPOINT classification (6 classes).** A piRNA (strain X, stage T, seq) is non-unique ONLY for
**same-stage** evidence elsewhere (BioMNI 3/3: stage-specific MILI/MIWI2 vs MIWI machinery → compare within
stage). Determinants: `ee_withintp_diag.csv.gz` (per-tp pools: same-/other-/no-stage expression elsewhere),
`loci/present_in_{Y}.bed` (halLiftover), `snp_variant_refinement_withintp.csv` (per-tp SNP), `cand_self16/*.bam`.
  - **expressed elsewhere (same stage)** — shared at stage → NOT unique (= 0 here: no candidate's exact seq is in another strain's same-tp pool)
  - **SNP-variant (1–3 mm, same stage)** — 1–3 mm allele of a same-stage-expressed piRNA → NOT unique
  - **low-quality** — no mm0 own-genome locus → NOT unique
  - **UNIQUE: strain-private locus** — no ortholog elsewhere (insertion-driven)
  - **UNIQUE: conserved-but-silent** — ortholog present, silent at T elsewhere (regulatory divergence)
  - **UNIQUE: stage-shifted (heterochronic)** — exact seq expressed elsewhere ONLY at a different stage
    (developmental-timing divergence; within its own stage it is strain-specific)

**Result:** **genuinely-unique within-tp = 10,724 (67 %)** = conserved-but-silent 6,560 + stage-shifted 2,725 +
strain-private 1,439. Not unique: SNP-variant 4,394 + low-quality 849. Strict-sequence subset (excl. heterochronic)
= 7,999. (`Fig_deseq_stagepeak_unique`)

**S6 · SNP-variant effect.** `Fig_snp_variant_effect.py`: SNP-refinement removes **4,394 strain alleles
(29 % of the naive unique set; −24/−30/−35 % per stage)** — these only look strain-specific but are 1–3 mm
alleles of same-stage piRNAs. 88 % are single-SNP; transition/transversion ≈ 1.37.

**S7 · PCA.** `Fig_pca_stagepeak.py`: per-tp PCA of the within-tp unique-piRNA expression across 48 libraries
(CPM by `libsize_window`, log2). **Wild-derived strains separate by their divergent unique repertoire**
(SPRET dominates PC1, up to 83 % at P20.5); the 12 classical strains cluster tightly.

**S8 · TE family of origin (feed the stage-shifted class in).** `Fig_te_origin_stageshift.py` (method = theme-08
`Fig_TE_private_families16`, extended to 3 classes): per strain, intersect each mechanism's loci (cand_self16
mm0, PanSN `#`-stripped → chrN) with the per-strain RepeatMasker BED → largest-overlap TE family. **Result (TE-
derived, lower bound):** strain-private **50 %** / **stage-shifted 45 %** / conserved-but-silent **31 %**. Strain-
private and stage-shifted are both dominated by young, active TEs — **LTR/ERVK 45 %/39 % + LINE/L1 17 %/18 %** —
whereas conserved-but-silent is more diverse/older (SINE/Alu, ERVL-MaLR). → **heterochronic (stage-shifted)
uniqueness is TE-driven by the same young ERVK/IAP + L1 families as insertion-driven uniqueness** (NOT the older
families of the regulatory class). **BioMNI 3/3** confirm (young ERVK/L1 = main source of recent strain-variable
piRNAs; timing shift plausibly via LTR-promoter activity / DNA-methylation timing / A-MYB pachytene activation).
(`Fig_te_origin_stageshift`)

**S9 · sense/antisense to TE (silencing-competence).** `Fig_sense_antisense_stageshift.py` (theme-08 convention:
orientation relative to the TE strand from RM `.out`; antisense = base-pairs the TE transcript = silencing-
competent; reuses the classification-independent per-locus orientation cache). **Result:** **stage-shifted 56 %
antisense** and **conserved-but-silent 57 %** are antisense-biased (silencing-competent), whereas **strain-private
(new insertions) 47 % ≈ 50/50 (dual-strand)**. So the heterochronic class is a **HYBRID** — young-TE ORIGIN
(ERVK/L1) like insertion-driven (S8), but the antisense SILENCING signature of conserved loci: shared young-TE
silencing loci expressed at different developmental stages across strains. **BioMNI 3/3** confirm (new insertions
dual-strand de-novo; conserved loci antisense-biased; stage-shifted antisense ⇒ silencing-competent).
(`Fig_sense_antisense_stageshift`)

## TOOLS
| Tool | What/why | Key params |
|---|---|---|
| DESeq2 1.50.2 (R, env biomni_e1) | DA engine (adopted) | `~0+strain`, one-vs-rest numeric contrast, padj<0.05 & log2FC>0 |
| edgeR 4.8.2 | benchmark comparator + `filterByExpr` | TMM on libsize_window; QL F-test |
| pysam / pandas / numpy | classification + figures | determinant reuse; Liberation Sans, vector |
| BioMNI (3 agents) | method + biology grounding | triple-verified (DESeq2 choice; within-tp; heterochrony) |

## INPUTS
`unique_pirna/edger16/{tp}.*` (counts/seqs/samples) · `unique16/loci/present_in_{Y}.bed` (halLiftover) ·
`unique16/snp_variant_refinement_withintp.csv` · `cand_self16/*.bam` · `unique16/ee_withintp_diag.csv.gz`.

## OUTPUTS (`figures/`, PDF+SVG+PNG + `.note.md`; data in `data/`)
`Fig_benchmark_da` (DESeq2 vs edgeR) · `Fig_filter_order` (27/30 before vs after) ·
`Fig_deseq_stagepeak_unique` (the unique set, 3 mechanisms) · `Fig_snp_variant_effect` ·
`Fig_pca_stagepeak` · `Fig_te_origin_stageshift` (TE family of origin) · `Fig_sense_antisense_stageshift` (TE strand / silencing-competence) · `Fig_unique_expression_heatmap` (expression of all unique piRNAs across strain × tp) · `Fig_unique_expression_heatmap_classical` / `_wild` (clade-split heatmaps; 554 classical vs 10,170 wild = 95% wild). Source data: `data/SourceData_*.csv`.

## VERDICT
For stage-characteristic strain-unique piRNAs: use **DESeq2** (calibrated; edgeR over-calls), filter to the
**27/30 nt stage peaks BEFORE DESeq2** (more sensitive), and judge uniqueness **within developmental stage**
(stage-specific biology). This yields **10,724 within-tp genuinely-unique piRNAs (67 %)** across three biological
mechanisms — insertion (strain-private), regulatory divergence (conserved-but-silent), and **heterochronic
timing divergence (stage-shifted)** — with SNP-refinement removing 29 % of naive-unique strain alleles. Wild-
derived strains (esp. SPRET) carry the bulk of unique piRNAs (PCA).

## REFERENCES (PMID + DOI verified against PubMed/publishers 2026-06-23)
- Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR. *Bioinformatics* 26(1):139–140. PMID **19910308** · DOI **10.1093/bioinformatics/btp616**
- Love MI, Huber W, Anders S (2014). DESeq2. *Genome Biology* 15:550. PMID **25516281** · DOI **10.1186/s13059-014-0550-8**
- Soneson C, Delorenzi M (2013). Comparison of DE methods for RNA-seq. *BMC Bioinformatics* 14:91. PMID **23497356** · DOI **10.1186/1471-2105-14-91**
- Schurch NJ et al. (2016). How many biological replicates / which DE tool. *RNA* 22(6):839–851. PMID **27022035** · DOI **10.1261/rna.053959.115**
- Aravin A et al. (2006). MILI-bound small RNAs in mouse testes. *Nature* 442:203–207. PMID **16751777** · DOI **10.1038/nature04916**
- Girard A, Sachidanandam R, Hannon GJ, Carmell MA (2006). piRNAs bind mammalian Piwi. *Nature* 442:199–202. PMID **16751776** · DOI **10.1038/nature04917**
- Li XZ et al. (2013). A-MYB initiates pachytene piRNA production. *Molecular Cell* 50:67–81. PMID **23523368** · DOI **10.1016/j.molcel.2013.02.016**

## SCRIPTS (full paths)
Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`PY=…/envs/biomni_e1/bin/python`, `R=…/envs/biomni_e1/bin/Rscript`).

*Compute (`analysis/claude_biomni_analysis/unique_pirna/`):*
- `benchmark_da_methods.R` + `run_bench_da.sh` — edgeR-vs-DESeq2 benchmark (SLURM array, 3 tp)
- `export_nullpv_hist.R` — null p-value histograms for the calibration panel
- `deseq_filter_order_test.R` + `run_filter_order.sh` — filter-order test; produces Order-A (production) candidates
- `classify_deseq_stagepeak.py` — within-tp 6-class classification (reuses determinants)
- `ee_withintp_diag.py` — per-tp same/other/no-stage expression-elsewhere diagnostic
- `build_pools16_pertp.py`, `classify_step416_pertp.py`, `aggregate_within_tp_snp.py` — per-tp pools + within-tp SNP set

*Figures (`figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA/code/`):*
- `Fig_benchmark_da.py` · `Fig_filter_order.py` · `Fig_deseq_stagepeak_unique.py` · `Fig_snp_variant_effect.py` · `Fig_pca_stagepeak.py` · `Fig_te_origin_stageshift.py` · `Fig_sense_antisense_stageshift.py` · `Fig_unique_expression_heatmap.py` · `Fig_unique_expression_heatmap_split.py` (classical + wild)
