# 13 — Divergence loci (genome-conserved cluster, strain-restricted expression)

**What these figures are.** Loci where the piRNA-cluster LOCUS is present in (nearly) all 16 strain genomes,
yet the cluster is **expressed in only a few strains** — cluster gain/loss by **regulatory/sequence
divergence**, not by locus gain. The contrasting mechanism to theme 09/12 (TE-insertion creation).

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1–S5 · reads → clusters → pangenome → TE/PAV → read layer.** identical to theme 11 §S1–S5
(cutadapt 5.0 / STAR 2.7.11b / PICB R 4.2.3 / cactus v2.9.3 halLiftover / bedtools 2.31.1 + RepeatMasker / pysam).

**S6 · divergence-locus selection.** **Python 3.11.15** — keep loci whose **genome PAV = present in ≥ most
strains** (halLiftover ●) but whose **pangenome FPM = expressed in only a few** strains (present-but-silent
elsewhere). **Result #:** **12** divergence loci (`divergence_loci_final.tsv`).

**S7 · figures.** matplotlib (Python 3.11.15) — `make_pav_locus_multi.py` (genome ●/○ vs expression FPM make
the divergence visible). **Result #:** **12** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| STAR / PICB | 2.7.11b / R 4.2.3 | reads + clusters | piRNA params; combined |
| cactus / halLiftover | v2.9.3 | genome PAV (locus present) vs pangenome FPM (expressed) | strain↔GRCm39 |
| bedtools | 2.31.1 | TE ∩ cluster | `intersect` |
| pysam | (Python 3.11.15) | per-strain coverage / 1U / antisense | 24–32 nt; TE-relative |
| Python | 3.11.15 | divergence selection + figures | present-genome ∩ silent-expression; matplotlib |

## INPUTS  `pangenome_te/divergence_loci_final.tsv` (12); `picb_pangenome_clusters.tsv`; BAMs; RepeatMasker; HAL.
## OUTPUTS (`figures/`, 12)  `Fig_divergence_pav_{locus}` (PDF+SVG+PNG + `.note.md`).

## DOUBLE-VERIFICATION
- Genome presence (halLiftover ●) and expression (FPM) are **separate** axes → divergence = present-but-silent,
  distinct from absence (the theme-03 lesson applied at the cluster level).
- Mechanism (divergence-driven, broad across wild+classical, ~83 % of strain-restricted loci) is
  **BioMNI-verified** (theme 10 / `project_te_driven_finding`); cross-ref theme 09 (insertion-driven).

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# project each strain's combined PICB clusters -> GRCm39 via cactus halLiftover:
python analysis/claude_biomni_analysis/unique_pirna/build_picb_pangenome_fpm.py   # -> picb_pangenome_clusters.tsv
# present-genome ∩ silent-expression selection -> make_pav_locus_multi.py below.
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/13_divergence_loci/code/make_pav_locus_multi.py
bash figures/analysis_figures/13_divergence_loci/code/render_divergence.sh
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/13_divergence_loci/code/make_pav_locus_multi.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/13_divergence_loci/code/pav_clusters.py`  _(imported helper)_
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/13_divergence_loci/code/render_divergence.sh`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/build_picb_pangenome_fpm.py` — cluster -> GRCm39 pangenome projection + FPM
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/run_picb_analysis_chunked.sh` — per-replicate PICB driver (cutadapt->STAR->PICB)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — PICB cluster calling (chunked, genome-wide LIBRARY.SIZE)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh` — combined (replicate-pooled) PICB driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_combine_script.R` — PICB on pooled BAM
