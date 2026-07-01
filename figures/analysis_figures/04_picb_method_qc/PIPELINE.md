# 04 — PICB method QC (library-size mode comparison + mapping complexity)

**What these figures are.** Methodological QC underpinning every PICB result (themes 02, 03, 05, 10):
(a) PICB FPM is only correct with a **genome-wide LIBRARY.SIZE**, and (b) the sRNA BAMs have the
mapping-complexity PICB assumes. Exists because of a real bug (memory `project_picb_concurrency_bug`).

---

## STEP-BY-STEP (tool · version · parameters · result)

**S0–S1 · reads → BAM.** cutadapt 5.0 + STAR 2.7.11b (piRNA params; unmasked genome) — as theme 02.

**S2 · three PICB processing modes.** **PICB (R 4.2.3)** on the same BAM, varying `LIBRARY.SIZE`:
Mode 1 = correct chr-by-chr (genome-wide size passed to each per-chr call); Mode 2 = naive chr-by-chr
(per-chromosome size → inflates FPM); Mode 3 = whole-BAM (genome-wide from `PICBload`).
**Result #:** genome-wide library size example LS = **11 277 390** reads; **Mode 1 ≈ Mode 3** (FPM agree),
Mode 2 diverges.

**S3 · mapping-complexity diagnostics.** **samtools 1.21 / Python 3.11.15** on STAR BAMs — unique-mapping
fraction, NH-tag distribution, multimapping breakdown, inflation ratio, piRNA-size enrichment.
**Result #:** **12** diagnostic panels (`mapping_complexity/01–12`) confirm the high-multimapper regime
piRNAs require (not an artefact).

**S4 · figures.** matplotlib (Python 3.11.15). → 2 mode-comparison figures + 12-panel complexity report.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| STAR | 2.7.11b | sRNA alignment | piRNA params (multimap 800) |
| PICB | R 4.2.3 | 3 library-size modes | genome-wide vs per-chr `LIBRARY.SIZE` |
| samtools | 1.21 | NH/multimap stats | flagstat / NH tag |
| Python | 3.11.15 | figures | matplotlib |

## INPUTS  `results/picb_result`, `results/picb_result_combined` (xlsx); STAR BAM NH stats.
## OUTPUTS (`figures/`)  `PICB_chr_vs_wholebam_comparison`, `PICB_chrbychr_vs_wholebam` + `mapping_complexity/01–12`.

## DOUBLE-VERIFICATION
- Mode 1 vs Mode 3 are two independent correct routes that must agree → the cross-check itself.
- xlsx completeness verified by chr1–19+X coverage (memory `project_picb_concurrency_bug`).
- This is **methods QC**, not a biological claim.

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# 3 PICB library-size modes (Mode1 chr-by-chr genome-wide, Mode2 naive, Mode3 whole-BAM):
bash analysis/claude_biomni_analysis/run_picb_comparison_slurm.sh
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/04_picb_method_qc/code/picb_chr_vs_whole_bam_comparison.py
$PY figures/analysis_figures/04_picb_method_qc/code/picb_comparison_figure.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/04_picb_method_qc/code/picb_chr_vs_whole_bam_comparison.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/04_picb_method_qc/code/picb_comparison_figure.py`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/run_picb_comparison_slurm.sh` — 3-mode PICB library-size comparison driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — chunked PICB (Mode 1/2)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_wholebam.R` — whole-BAM PICB (Mode 3)
