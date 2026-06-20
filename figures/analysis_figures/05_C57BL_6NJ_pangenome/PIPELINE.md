# 05 — C57BL/6NJ pangenome starter analysis (P12.5 / P20.5)

**What these figures are.** The single-strain **pilot** that established the pangenome workflow later scaled
to 16 strains: developmental classification of C57BL/6NJ piRNA clusters + their expression, SV, TE and
Zamore-conservation content (Wong-2011 colourblind palette).

---

## STEP-BY-STEP (tool · version · parameters · result)

**S0–S2 · clusters.** cutadapt 5.0 → STAR 2.7.11b (unmasked C57BL/6NJ, piRNA params) → **PICB (R 4.2.3)**
at P12.5 and P20.5. **Result #:** P12.5 merged = **9 728** clusters; P20.5 merged = **1 968** clusters.

**S3 · developmental classification.** **Python 3.11.15 / bedtools 2.31.1** — classify each cluster by
timepoint expression. **Result #:** **3** dev classes (P12.5_only / shared_postnatal / pachytene).

**S4 · TE + SV content.** intersect clusters with **RepeatMasker** (C57BL/6NJ BED) and pangenome SVs.
**Result #:** **2 727** TE-sized SVs at clusters (GRCm39); per-class TE fraction in `C57BL_6NJ_TE_annotation.csv`.

**S5 · Zamore conservation.** Zamore loci (**mm10→mm39 UCSC liftOver**) ∩ C57BL/6NJ clusters.

**S6 · figures.** matplotlib (Python 3.11.15). → **7** figures (Fig1–7).

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| STAR | 2.7.11b | unmasked C57BL/6NJ alignment | piRNA params |
| PICB | R 4.2.3 | cluster calling (P12.5/P20.5) | `clusters` sheet |
| bedtools | 2.31.1 | dev-class / TE / SV / Zamore overlaps | `intersect` |
| RepeatMasker BED | per-strain | TE class/family | col4 `name|class/family` |
| UCSC liftOver | 2023 | Zamore mm10→mm39 | `mm10ToMm39.over.chain.gz` |
| Python | 3.11.15 | figures | matplotlib (Wong palette) |

## INPUTS  `results/picb_result*/C57BL_6NJ-*`; RepeatMasker BED; Zamore loci; pangenome SVs → `data/C57BL_6NJ_pangenome/*.csv`.
## OUTPUTS (`figures/`)  Fig1 developmental_classification · Fig2 FPM_expression · Fig3 pangenome_SV_content · Fig4 TE_content · Fig5 timepoint_TE_analysis · Fig6 Zamore_liftover · Fig7 SV_expression.

## DOUBLE-VERIFICATION
- Numbers recomputed from `C57BL_6NJ_pangenome/*.csv` (P12.5 9728 / P20.5 1968 clusters; 3 dev classes; 2727 SVs).
- mm10→mm39 liftover enforced for Zamore; TE class from RepeatMasker (not inferred).
- This pilot's conclusions were re-tested at 16-strain scale (themes 03, 10) — concordant direction.

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# PICB clusters for C57BL/6NJ P12.5 & P20.5:
bash workflow/scripts/run_picb_analysis_chunked.sh
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/05_C57BL_6NJ_pangenome/code/C57BL_6NJ_pangenome_figures.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/05_C57BL_6NJ_pangenome/code/C57BL_6NJ_pangenome_figures.py`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/run_picb_analysis_chunked.sh` — per-replicate PICB driver (cutadapt->STAR->PICB)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — PICB cluster calling (chunked, genome-wide LIBRARY.SIZE)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh` — combined (replicate-pooled) PICB driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_combine_script.R` — PICB on pooled BAM
