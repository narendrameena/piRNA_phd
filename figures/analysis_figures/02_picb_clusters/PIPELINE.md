# 02 — PICB piRNA clusters per strain (counts & replicate reproducibility)

**What these figures are.** How many piRNA clusters PICB calls per strain × timepoint, and whether pooling
replicates before PICB changes **which** clusters are found (reproducibility QC). PICB = piRNA Cluster Builder.

---

## STEP-BY-STEP: raw data → figure (tool · version · parameters · result)

**S0 · trim.** **cutadapt 5.0** — `TGGAATTCTCGGGTGCCAAGG`, 20–36 nt, `--discard-untrimmed`.

**S1 · align.** **STAR 2.7.11b** on **unmasked** strain genome (piRNA params: `--outFilterMismatchNmax 0
--outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --alignIntronMax 1 --alignEndsType EndToEnd`).
→ per-replicate + replicate-pooled BAMs.

**S2 · PICB cluster calling.** **PICB (R 4.2.3 / Bioconductor)**, genome-wide `LIBRARY.SIZE`; `clusters`
sheet. Two modes: per-replicate (rep 1/2/3) and **combined** (reps pooled before PICB).
**Result #:** **48** combined runs (16 strains × 3 tp); median **12 713** clusters (range 2 174–25 958);
combined ≈ mean of the 3 replicates.

**S3 · reproducibility overlap.** **bedtools 2.31.1** — same-strand interval overlap of replicate vs combined
cluster BEDs; support composition (in 3/2/1/0 reps). **Result #:** a single replicate already recovers
**92.1 %** (median) of combined clusters; **<1 %** of combined clusters are unique to combined.

**S4 · figures.** **matplotlib (Python 3.11.15)** — grouped bars (**error bar = ±SD across 3 replicates**),
single-rep-vs-combined, and the support-composition / recovery overlap. → **3** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| cutadapt | 5.0 | trim adapter | `TGGAATTCTCGGGTGCCAAGG`, 20–36 nt |
| STAR | 2.7.11b | unmasked alignment | piRNA params (0 mismatch, 800 multimap, EndToEnd) |
| PICB | R 4.2.3 | cluster calling | `clusters` sheet; combined = reps pooled |
| bedtools | 2.31.1 | rep↔combined overlap | same-strand interval intersect |
| Python | 3.11.15 | figures | matplotlib, error bars = ±SD |

## INPUTS
`results/picb_result/{strain}/{strain}-{tp}.{rep}.xlsx` (per-rep) + `results/picb_result_combined/…combined.xlsx`
→ `source_data/SourceData_PICB_cluster_counts.csv` (replicate 1/2/3/combined),
`SourceData_PICB_rep_combined_overlap.csv`.

## OUTPUTS (`figures/`, PDF+SVG+PNG + `.note.md`)
`Fig_picb_combined_clusters` (counts/strain, ±SD error bars) · `Fig_picb_rep_vs_combined` ·
`Fig_picb_rep_combined_overlap` (support composition + recovery).

## DOUBLE-VERIFICATION
- Numbers recomputed from `SourceData_PICB_cluster_counts.csv`: 48 combined, median 12 713 (2 174–25 958);
  per-rep example C57BL/6NJ E16.5 = 19933/18810/18786, combined 18855 (combined ≈ rep mean → ±SD error bars valid).
- The "combining replicates doesn't change which clusters are found" conclusion is read directly from the
  >99 % overlap — a **method-QC** result, not a biological claim. Cross-ref theme 04 (PICB library-size QC).

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# PICB clusters per replicate (cutadapt -> STAR -> PICB; calls picb_script_chunked.R):
bash workflow/scripts/run_picb_analysis_chunked.sh
# PICB on replicate-pooled BAM (combined run; calls picb_combine_script.R):
bash analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh
# S3 replicate<->combined overlap is computed inside the figure scripts (interval overlap).
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/02_picb_clusters/code/Fig_picb_combined_clusters.py
$PY figures/analysis_figures/02_picb_clusters/code/Fig_picb_rep_combined_overlap.py
$PY figures/analysis_figures/02_picb_clusters/code/Fig_picb_rep_vs_combined.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/02_picb_clusters/code/Fig_picb_combined_clusters.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/02_picb_clusters/code/Fig_picb_rep_combined_overlap.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/02_picb_clusters/code/Fig_picb_rep_vs_combined.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/02_picb_clusters/code/strain_order.py`  _(imported helper)_

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/run_picb_analysis_chunked.sh` — per-replicate PICB driver (cutadapt->STAR->PICB)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — PICB cluster calling (chunked, genome-wide LIBRARY.SIZE)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh` — combined (replicate-pooled) PICB driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_combine_script.R` — PICB on pooled BAM
