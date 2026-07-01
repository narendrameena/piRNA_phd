# 14 — Circos: TE expression vs piRNA silencing across the 16-strain pangenome

**What these figures are.** 10 genome-wide circos plots; rings = the 16 strains, axis = GRCm39 in 2 Mb bins.
Per TE family and bin they contrast **TE expression** (sense-to-TE small RNA / RNA) against **piRNA silencing**
(antisense-to-TE piRNA), genome-wide.

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1 · reads/clusters.** cutadapt 5.0 → STAR 2.7.11b (unmasked) → PICB (R 4.2.3); RNA-seq for TE expression
(featureCounts on `results/TE_expression_rna/`).

**S2 · pangenome frame.** **cactus v2.9.3 halLiftover** → GRCm39 2 Mb bins.

**S3 · TE family signal per bin.** **bedtools 2.31.1** ∩ **RepeatMasker**; per bin × TE family compute
**sense-to-TE** (expression) and **antisense-to-TE** (silencing piRNA) — strand is **TE-relative**, never
genomic ±. **Result #:** GRCm39 binned at **2 Mb**; 16 strain rings.

**S4 · figures.** matplotlib (Python 3.11.15) circos suite. **Result #:** **10** circos figures (per-strain +
combined TE×ping-pong `Fig_circos_combined16`).

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| STAR / PICB | 2.7.11b / R 4.2.3 | sRNA reads + clusters | piRNA params |
| featureCounts | (subread) | TE-family RNA expression | per-family Geneid |
| cactus / halLiftover | v2.9.3 | project to GRCm39 bins | strain→GRCm39 |
| bedtools | 2.31.1 | bin × TE-family overlap | `intersect`; 2 Mb bins |
| RepeatMasker BED | per-strain | TE family + strand | col4 + strand |
| Python | 3.11.15 | circos figures | matplotlib |

## INPUTS  per-strain sRNA BAMs; `results/TE_expression_rna/*.TE.featureCounts`; RepeatMasker; HAL → `data/` (binned TE matrices).
## OUTPUTS (`figures/`, 10)  `Fig_circos_*16` (per-strain coverage + combined TE × ping-pong).

## DOUBLE-VERIFICATION
- sense = TE expression, antisense = piRNA silencing — strand **relative to the TE** (RepeatMasker), the
  corrected definition (`session_resume` lesson); genomic strand = cluster architecture, not sense/antisense.
- TE expression (RNA) and piRNA silencing (sRNA) are independent measurements (cross-checked).

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# S1 TE RNA expression (featureCounts) + sRNA PICB; S2 project clusters/bins -> GRCm39:
python analysis/claude_biomni_analysis/unique_pirna/build_picb_pangenome_fpm.py
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_1u16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_clusters16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_combined16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_coverage16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_integrated16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_picb16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_pingpong16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_private16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_te16.py
$PY figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_te_pirna16.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_1u16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_clusters16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_combined16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_coverage16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_integrated16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_picb16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_pingpong16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_private16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_te16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/code/Fig_circos_te_pirna16.py`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/build_picb_pangenome_fpm.py` — cluster -> GRCm39 pangenome projection + FPM
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/run_picb_analysis_chunked.sh` — per-replicate PICB driver (cutadapt->STAR->PICB)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — PICB cluster calling (chunked, genome-wide LIBRARY.SIZE)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh` — combined (replicate-pooled) PICB driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_combine_script.R` — PICB on pooled BAM
