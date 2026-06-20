# 12 — TE-creation source loci (strain-private TE insertion → new piRNA cluster)

**What these figures are.** One figure per **strain-private TE-insertion locus** where a TE inserted in a
single (carrier) strain and a piRNA cluster formed there — shown across all 16 strains (genome PAV ●/○ +
pangenome FPM). Same 3-panel anatomy as theme 11, windowed on the insertion.

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1–S5 · reads → clusters → pangenome → TE/PAV → read layer.** identical to theme 11 §S1–S5
(cutadapt 5.0 / STAR 2.7.11b / PICB R 4.2.3 / cactus v2.9.3 halLiftover / bedtools 2.31.1 + RepeatMasker / pysam).

**S2b · carrier-cluster window SNAP fix.** the strain-private TE-insertion master table gives a one-directional
breakpoint window that can miss the adjacent cluster; **Python 3.11.15** gap-merges (gap+8 kb) the carrier's
pangenome clusters on the locus chromosome and snaps the window to the best-overlapping cluster group
(±3 kb; 45 kb cap → fallback breakpoint ±12 kb). **Result #:** fixed **9/27** loci that previously KeyError'd;
all **27** verified read-based ≥99 % carrier coverage.

**S6 · figures.** matplotlib (Python 3.11.15) — `make_source_pav_multi.py` (+ `make_source_pav.py` single).
EXPR_MIN = 50. **Result #:** **27** source loci → **29** figures (incl. CAST ERVL chr10 single+multi).

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| STAR / PICB | 2.7.11b / R 4.2.3 | reads + clusters | piRNA params; combined |
| cactus / halLiftover | v2.9.3 | pangenome projection + genome PAV | strain↔GRCm39 |
| bedtools | 2.31.1 | gap-merge + insertion ∩ cluster | `merge`(gap 8 kb)/`intersect` |
| pysam | (Python 3.11.15) | carrier read coverage / 1U / antisense | 24–32 nt; TE-relative |
| Python | 3.11.15 | window snap + figures | gap+8 kb, ±3 kb, 45 kb cap; matplotlib |

## INPUTS  `pangenome_te/source_loci_master_creation_private.tsv` (carrier,chrom,start,end,te,strand,FPM,class,gbreadth); `picb_pangenome_clusters.tsv`; BAMs; RepeatMasker; HAL.
## OUTPUTS (`figures/`, 29)  `Fig_source_pavmulti_{carrier}_{TE}_{chr}_{pos}` (+ single).

## DOUBLE-VERIFICATION
- The window-snap fix was verified **read-based** (all 27 carriers ≥99 % coverage) after the KeyError diagnosis.
- Simple-union snap over-broadening (LP_J 223 kb) caught + capped (45 kb → breakpoint ±12 kb fallback).
- TE-creation causal claim = **[finding]** (shares theme 09's status); 1U/antisense = *established*.

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# project each strain's combined PICB clusters -> GRCm39 via cactus halLiftover:
python analysis/claude_biomni_analysis/unique_pirna/build_picb_pangenome_fpm.py   # -> picb_pangenome_clusters.tsv
# source-loci master table -> windowed 3-panel views (make_source_pav*.py below).
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/12_creation_source_loci/code/make_source_pav.py
$PY figures/analysis_figures/12_creation_source_loci/code/make_source_pav_multi.py
bash figures/analysis_figures/12_creation_source_loci/code/render_srccreation.sh
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/12_creation_source_loci/code/make_source_pav.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/12_creation_source_loci/code/make_source_pav_multi.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/12_creation_source_loci/code/pav_clusters.py`  _(imported helper)_
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/12_creation_source_loci/code/render_srccreation.sh`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/build_picb_pangenome_fpm.py` — cluster -> GRCm39 pangenome projection + FPM
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/run_picb_analysis_chunked.sh` — per-replicate PICB driver (cutadapt->STAR->PICB)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — PICB cluster calling (chunked, genome-wide LIBRARY.SIZE)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh` — combined (replicate-pooled) PICB driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_combine_script.R` — PICB on pooled BAM
