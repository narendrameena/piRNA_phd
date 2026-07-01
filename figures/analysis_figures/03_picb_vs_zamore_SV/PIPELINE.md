# 03 — Zamore pachytene piRNA loci across 16 strains (PANGENOME, reframed 2026-06-17)

**What these figures are.** Cross-strain analysis of the **214 Zamore conserved pachytene piRNA-precursor
loci**: structural **retention**, **expression** (developmental + per-strain), expression **concentration**,
and **sequence-divergence ↔ expression-variability** — all on the **verified combined-run PICB clusters**,
with cross-strain projection through the **cactus pangenome graph (halLiftover)**, not UCSC chain liftOver.

> **Why reframed.** The 2026-05-21 version used UCSC `liftOver` and concluded "~20 % of conserved piRNA loci
> are structurally disrupted (SV-driven)." That was an **artifact of UCSC's ≥95 %-identity requirement**
> (it scores sequence-*diverged* loci as *absent*): the pangenome shows only **0.6 %** are truly
> not-projected and **98.8 % are retained**. Pachytene loci = conserved position, divergent sequence
> (Yu 2021 *Nat Commun* PMID 33397987; Li 2013 *Mol Cell* PMID 23523368). Old figures kept in
> `superseded_UCSC_liftover/`.

---

## STEP-BY-STEP: raw data → figure (tool · version · parameters · result)

> Every number below was **recomputed from the saved data files** to double-verify (see end).

**S0 · sRNA adapter trimming.** Tool **cutadapt 5.0**. Per strain × timepoint × replicate FASTQ.
Params: 3′ adapter `TGGAATTCTCGGGTGCCAAGG`, length filter `-m 20 -M 36`, `--discard-untrimmed`.
→ piRNA-length (20–36 nt) trimmed reads.

**S1 · sRNA alignment to each strain's UNMASKED genome.** Tool **STAR 2.7.11b**.
Params (piRNA-critical): `--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600
--alignIntronMax 1 --alignEndsType EndToEnd --scoreDelOpen -10000 --scoreInsOpen -10000` (zero mismatch,
highly permissive multimappers, no splicing/indels). → per-strain×tp BAM on REL-2205 unmasked assembly.

**S2 · PICB cluster calling (COMBINED run).** Tool **PICB (R 4.2.3 / Bioconductor)** on the
**replicate-pooled** BAM; FPM uses genome-wide `LIBRARY.SIZE`. → **48** verified
`*-{tp}.combined.xlsx` (16 strains × 3 tp), each spanning chr1–19+X; `clusters` sheet.
**Result #:** C57BL/6NJ P20.5 = **2 660** clusters.

**S3 · combined cluster BEDs.** Tool **Python 3.11.15 / openpyxl** (`clusters` sheet → BED, start−1,
bare-numeric chr, std chroms). → `combined_beds/{strain}_{tp}.bed`. **Result #: 48/48** BEDs written.

**S4 · Zamore reference loci (the ONLY UCSC liftOver).** Tool **UCSC `liftOver`** (2023 binary).
Input = published mm10 (GRCm38) pachytene annotation; chain `mm10ToMm39.over.chain.gz` (default
`-minMatch 0.95`). → **467** raw mm10 entries → grouped/stage-filtered to **214** GRCm39 loci
(`_loci_mm39_noprefix.bed`).

**S5 · cross-strain projection = PANGENOME halLiftover.** Tool **cactus v2.9.3 `halLiftover`**
(`singularity exec --bind /mnt <SIF> halLiftover <HAL> GRCm39 <bed> <strain> <out>`); each locus tiled into
**100-bp windows**; `fraction_lifted` = % of windows that project = graded retention.
**Result #:** **98.8 %** of locus×strain combinations retained (≥50 % span); most-divergent strain
SPRET/EiJ mean span retained **94.3 %**. **QC vs UCSC:** pangenome lifts **99.4 %** of loci vs UCSC
**79.8 %** — UCSC drops to **52/214** for SPRET (the artifact).

**S6 · expression call.** Tool **bedtools 2.31.1 `intersect -wa -u`** (projected locus ∩ that strain's
S3 combined clusters); quantitative FPM = Σ overlapping cluster `all_primary_FPM` from the pangenome
projection table. → status expressed / not_expressed / not_lifted.
**Result #:** expressed rises **E16.5 77.6 % → P12.5 85.8 % → P20.5 97.6 %** (3 342 loci-strain at P20.5);
true not_lifted = **0.6 %**.

**S7 · structural variants.** Tool **bcftools 1.21 `view -H -R`** on the 17-strain pangenome VCF
(tabix-indexed; GRCm39 frame); SV = INS/DEL **≥300 bp**, region **±50 kb** of each locus; intersect
**bedtools `window`/`intersect`**. → **Result #:** SV matrix **10 272** rows (214×16×3 windows);
CAST/EiJ direct-overlap SV loci = **119**.

**S8 · TE + genomic-region annotation.** TE: locus projected C57BL_6NJ via **halLiftover** ∩
**RepeatMasker** (C57BL_6NJ BED, col4 `name|class/family`). Region: locus ∩ **Ensembl GRCm39.115 gff3**
`biotype=`. → **Result #:** dominant TE class **SINE 108 / LTR 72 / LINE 28 / none 6** loci; mean TE
fraction **0.278**.

**S9 · statistics.** Tool **scipy** (`spearmanr`, `mannwhitneyu`, `kruskal`) + Lorenz/Gini.
→ **Result #:** Gini **0.903**; **25** loci → 90 % of output; cross-strain expression-rank Spearman
**ρ̄ = 0.88**; divergence↔CV Spearman classical **ρ=0.40** (p=1e-9), wild ρ=0.15; wild vs classical Gini
Mann–Whitney **p=0.042**. TE-class→divergence **Kruskal–Wallis p=0.38 → REJECTED, not claimed.**

**S10 · figures.** Tool **matplotlib (Python 3.11.15, Liberation Sans, vector PDF/SVG + 300 dpi PNG)**.
→ **5** figures (below).

---

## TOOLS (tool · version · what & why · key parameters)

| Tool | Version | What / why | Key parameters |
|---|---|---|---|
| cutadapt | 5.0 | trim sRNA 3′ adapter | `TGGAATTCTCGGGTGCCAAGG`; `-m 20 -M 36`; `--discard-untrimmed` |
| STAR | 2.7.11b | align sRNA to **unmasked** strain genome | `--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --alignIntronMax 1 --alignEndsType EndToEnd --scoreDelOpen/InsOpen -10000` |
| PICB | R 4.2.3 (Bioconductor) | piRNA clusters on **pooled** BAM; genome-wide FPM | `clusters` sheet; combined run |
| UCSC liftOver | 2023 binary | **only** mm10→mm39 of Zamore loci | chain `mm10ToMm39.over.chain.gz`; `-minMatch 0.95` |
| cactus / halLiftover | **v2.9.3** (SIF) | **pangenome** cross-strain projection | `singularity exec --bind /mnt … halLiftover <HAL> <src> <bed> <tgt> <out>`; 100-bp tiling |
| bcftools | 1.21 | SVs from 17-strain pangenome VCF | `view -H -R`; SV ≥300 bp; ±50 kb |
| bedtools | 2.31.1 | expression / TE / region / SV overlaps | `intersect -wa -u`, `-wa -wb -f 0.10`; `window -w` |
| RepeatMasker BED | per-strain REL-2205 | TE class/family | col4 `name|class/family` |
| Ensembl gff3 | GRCm39.115 | genomic-region class | `gene`/`ncRNA_gene` `biotype=` |
| Python | 3.11.15 (`biomni_e1`) | analysis + figures | pandas, numpy, scipy, matplotlib, openpyxl |

## INPUTS
`results/picb_result_combined/*/*.combined.xlsx` (48) · `…/C57BL_6NJ_pangenome/zamore_mm10.bed` (mm10) +
`_loci_mm39_noprefix.bed` (214) · `results/pangenome/output/mouse_17strain_pangenome.full.hal` ·
`…/mouse_17strain_pangenome.raw.vcf.gz`(+tbi) · `unique_pirna/cluster_pav/picb_pangenome_clusters.tsv` ·
`resources/repeatMasker/C57BL_6NJ_repeatmasker.bed` · `resources/black6/genome/Mus_musculus.GRCm39.115.chr.gff3` ·
`results/liftOverChainFiles/mm10ToMm39.over.chain.gz`

## OUTPUTS (5 figures, in `figures/`; PDF+SVG+PNG + `.note.md`)
| Figure | Headline result |
|---|---|
| `Fig_zamore_expression_pangenome` | expression 78→86→**98 %** into pachytene |
| `Fig_zamore_retention_pangenome` | **98.8 %** retained; QC pangenome 99.4 % vs UCSC 79.8 % |
| `Fig_zamore_expr_divergence_tests` | Gini **0.90**; 25 loci→90 %; ρ̄=0.88; div↔CV classical ρ=0.40 |
| `Fig_zamore_stage_time_heatmap` | wave: 99 pachytene loci fire at P20.5 |
| `Fig_zamore_retention_heatmap` | divergent loci concentrate in wild strains (SPRET) |

Data tables: `data/` (matrices `*_PANGENOME.csv`, `zamore_fraction_lifted.csv`, `zamore_locus_annotation.csv`,
`zamore_locus_expression_P20.5.csv`, `SourceData_*.csv`); stats log `data/zamore_expression_tests.txt`.

## DOUBLE-VERIFICATION
- **Every result number above was recomputed independently from the saved data files** (2026-06-17) and
  matches the figures: 214 loci, 98.8 % retained, 77.6/85.8/97.6 % expressed, not_lifted 0.6 %, pangenome
  99.4 % vs UCSC 79.8 % (SPRET 52/214), Gini 0.903, 25 loci→90 %, ρ̄ 0.88, dom-TE SINE 108/LTR 72/LINE 28.
- **Method steered by user:** UCSC liftOver used *only* mm10→mm39; all cross-strain projection = pangenome.
- **Biology triple-checked:** files + **BioMNI genomics** (mm10 build; ~100 pachytene loci; staging;
  few-loci-dominate) + the **papers** (`papers/Zamore_lab/`; Li 2013 PMID 23523368, Yu 2021 PMID 33397987).
- **Tested, not eyeballed:** TE-class→divergence rejected (Kruskal–Wallis p=0.38) and **not claimed**.
- **Caveat:** locus FPM = Σ overlapping clusters (aggregate; rank/scale conclusions are robust).
- BioMNI literature agent was degenerate (no PMIDs) → citations confirmed from the papers folder
  (`VERIFICATION_QUEUE.md`).

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# S2 PICB combined clusters:
bash analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh
# S4 Zamore loci mm10->mm39 (ONLY UCSC use):
liftOver zamore_mm10.bed mm10ToMm39.over.chain.gz loci_mm39.bed unmapped -minMatch=0.95
# S5 cross-strain projection (pangenome, per locus x strain):
singularity exec --bind /mnt <cactus.sif> halLiftover <HAL> GRCm39 <loci.bed> <strain> <out.bed>
# S7 SVs from 17-strain pangenome VCF:
bcftools view -H -R loci.bed results/pangenome/output/mouse_17strain_pangenome.raw.vcf.gz
# S3/S6/S8-S10 BEDs, expression, TE/region, stats + figures = the code/ scripts below.
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_PICB_vs_Zamore.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_SV_TE.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_SV_mechanism.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_expr_divergence_tests.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_expression_pangenome.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_retention_heatmap.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_retention_pangenome.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_stage_time_heatmap.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/all_strains_SV_figure.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/all_strains_all_timepoints_figure.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/annotate_zamore_loci.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/build_replicate_pct_expressed.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/pangenome_fraction_lifted.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/rebuild_picb_COMBINED.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/rebuild_zamore_COMBINED.py
$PY figures/analysis_figures/03_picb_vs_zamore_SV/code/zamore_expression_tests.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_PICB_vs_Zamore.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_SV_TE.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_SV_mechanism.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_expr_divergence_tests.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_expression_pangenome.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_retention_heatmap.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_retention_pangenome.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/Fig_zamore_stage_time_heatmap.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/all_strains_SV_figure.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/all_strains_all_timepoints_figure.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/annotate_zamore_loci.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/build_replicate_pct_expressed.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/pangenome_fraction_lifted.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/rebuild_picb_COMBINED.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/rebuild_zamore_COMBINED.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/strain_order.py`  _(imported helper)_
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/03_picb_vs_zamore_SV/code/zamore_expression_tests.py`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/run_picb_analysis_chunked.sh` — per-replicate PICB driver (cutadapt->STAR->PICB)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — PICB cluster calling (chunked, genome-wide LIBRARY.SIZE)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh` — combined (replicate-pooled) PICB driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_combine_script.R` — PICB on pooled BAM
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/build_picb_pangenome_fpm.py` — cluster -> GRCm39 pangenome projection + FPM
