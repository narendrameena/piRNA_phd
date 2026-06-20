# 06 — Zamore pachytene-locus coverage (C57BL/6 P12.5 / P20.5)

**What these figures are.** Direct sRNA-seq **coverage** of the Zamore conserved pachytene piRNA loci +
PICB cluster architecture — a recovery/QC check that the pipeline sees the canonical pachytene piRNA genes.

> **Note:** C57BL/6 (not 6NJ) = **external public data**; per memory `project_black6_qc` the P12.5/P20.5
> "replicates" are byte-identical duplicate files — depth/replication treated with caution (QC only).

---

## STEP-BY-STEP (tool · version · parameters · result)

**S0–S1 · reads → genome.** cutadapt 5.0 (C57BL/6 public adapters) → STAR 2.7.11b (unmasked GRCm39, piRNA
params).

**S2 · Zamore loci.** published pachytene annotation, **mm10 → mm39 UCSC liftOver**, grouped by stage.
**Result #:** **214** Zamore genes — **Pachytene 99 / Prepachytene 83 / Hybrid 32**.

**S3 · coverage.** **bedtools 2.31.1 / samtools 1.21** — sRNA BAM coverage per Zamore gene + over PICB
clusters; detection rate per stage; per-stage CDFs; pachytene heatmap. **Result #:** detection rates for
**4** stage groups; top **15** pachytene loci by FPM (`Fig3_top_pachytene_loci_FPM.csv`).

**S4 · figures.** matplotlib (Python 3.11.15). → **3** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| cutadapt | 5.0 | trim C57BL/6 public sRNA | TruSeq adapters (see `project_black6_qc`) |
| STAR | 2.7.11b | align to GRCm39 | piRNA params |
| UCSC liftOver | 2023 | Zamore mm10→mm39 | `mm10ToMm39.over.chain.gz` |
| bedtools / samtools | 2.31.1 / 1.21 | per-gene coverage | `coverage`/`intersect` |
| Python | 3.11.15 | figures | matplotlib |

## INPUTS  C57BL/6 sRNA BAMs; Zamore loci (mm39) → `data/Fig{1,2,3}_*.csv`, `P12_5_P20_5_zamore_coverage_per_gene.csv`.
## OUTPUTS (`figures/`)  Fig1 PICB_cluster_architecture · Fig2 Zamore_gene_coverage · Fig3 coverage_detail.

## DOUBLE-VERIFICATION
- Stage counts recomputed (Pachytene 99 / Prepachytene 83 / Hybrid 32) — identical to themes 03 & 15.
- mm10→mm39 liftover enforced; coverage from BAMs (not inferred).
- **Caveat:** external C57BL/6 data with duplicate-file "replicates" — recovery QC only, excluded from thesis.

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# cutadapt -> STAR (GRCm39 unmasked, piRNA params); Zamore loci mm10->mm39:
liftOver zamore_mm10.bed mm10ToMm39.over.chain.gz zamore_mm39.bed unmapped -minMatch=0.95
# per-gene + per-cluster coverage (bedtools/samtools) is done inside the figure scripts below.
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/06_zamore_coverage/code/Fig_black6_pirna_qc_persample.py
$PY figures/analysis_figures/06_zamore_coverage/code/generate_P12_5_P20_5_figures.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/06_zamore_coverage/code/Fig_black6_pirna_qc_persample.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/06_zamore_coverage/code/generate_P12_5_P20_5_figures.py`
