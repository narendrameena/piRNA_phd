# 17 — PICB clusters vs Trinity precursors: which better identifies piRNA source loci?

**What this theme is.** A head-to-head, data-driven comparison of the two precursor-finding strategies used in
this project: **PICB** (piRNA Cluster Builder — a *reference-genome* cluster caller on mapped sRNA; theme 02) vs
**Trinity** (*de-novo* RNA-seq assembly → contigs ≥500 bp covered by sRNA; theme 16). "Better" is judged on
multiple axes, not one number: **(i)** how many loci each finds, **(ii)** reciprocal genomic concordance
(specificity + sensitivity), **(iii)** fragmentation / locus length, and **(iv)** the decisive test — **how much
of the actual piRNA output (read mass) each locus set captures**.

> **Build compatibility — VERIFIED before any intersection** (mandatory, per project rule). Both sets live in the
> SAME strain REL-2205 assembly: chr1 length **194,686,469** is identical between the PanSN minimap2 target
> (`resources/REL-2205-Assembly/{sid}_chromosomes_MT.fasta`, headers `{sid}#1#chrN`) and the PICB refFasta
> (`resources/PICB/refFasta/{sid}_chromosomes_MT.fasta`, headers `chrN`). Only the chromosome-name convention
> differs → strip `{sid}#1#chr` ↔ `N` to intersect. No liftover, no build mix.

---

## STEP-BY-STEP: raw data → figure

**S1 · PICB cluster intervals.** from `figures/analysis_figures/_shared_data/picb_pangenome_clusters.tsv`
(combined-replicate PICB; theme 02/03). Use the strain-**own** coordinates (`own_chrom/own_start/own_end`),
de-duplicate, `bedtools merge` per strain × tp → PICB clusters.

**S2 · Trinity precursor intervals.** 100/100 contigs (**RPM>100 & RPKM>100**, theme 16) per sample → their
genomic loci from `results/filter_precursors_bed/{sid}/{sid}-{tp}.{rep}.bed` (minimap2-mapped; 99.5 % of 100/100
precursors map). Strip PanSN prefix, union the 3 replicates, `bedtools merge` → Trinity precursor loci.

**S3 · reciprocal overlap (`bedtools intersect`).** per strain × tp: % Trinity loci overlapping any PICB cluster;
% PICB clusters overlapping any Trinity locus; fragmentation = distinct Trinity contigs per recovered cluster;
bp coverage; median locus length. → `picb_vs_trinity/overlap_per_strain_tp.csv` (`compute_overlap.py`).
**Results (16 strains, systematic):** PICB **~2–4× more loci** (E16.5 ~15,770 / P12.5 ~12,275 / P20.5 ~3,006 vs
Trinity ~3,816 / 7,147 / 1,490); **~59–69 % of Trinity precursors overlap a PICB cluster** (~1/3 off-cluster);
Trinity **recovers ~10 % / 18 % / 9 % of PICB clusters**; fragmentation **2.6× → 5.7×** (worst at pachytene);
PICB median length ~2,600–3,100 bp vs Trinity ~1,300–1,900 bp.

**S4 · DECISIVE — piRNA read capture.** of **total piRNA** (25–32 nt reads, multimapper-weighted; `samtools view -L`
+ length filter), fraction mapping inside PICB clusters vs Trinity precursor **EXON blocks** (bed12tobed6 — NOT the
intron-spanning genomic span, which inflates Trinity 5–7×; see `Fig_capture_methodology_test`). Subset: pachytene
P20.5 ×3 (SPRET/CAST/C57BL_6NJ) + P12.5 ×2. → `picb_vs_trinity/read_capture_pirna.csv`.
**Result (thesis 100/100 filter; 25–32 nt piRNA; all 16 strains, 46/48 — C3H_HeJ/NOD_ShiLtJ E16.5 totals pending a
65 GB recompute):** **PICB ≥ Trinity, developmental:** **E16.5 PICB ~4× Trinity** (19.9 vs 4.8 %), **P12.5 ~2×**
(8.9 vs 4.4 %), **P20.5 ~tied** (26.3 vs 26.1 %; at 100/100 Trinity catches up even more — CAST P20.5 Trinity 20.3 %
slightly exceeds PICB). 5-sample piRNA-specific (25–32 nt): P20.5 11.2/11.0, 19.9/20.3, 26.0/25.6 %; P12.5 7.4/4.8,
8.8/3.8 %. NB sets overlap → not additive. Caveat: 24 nt is 1U-impure (46 % vs core 80 %) → 25–32 is the cleaner window.

**S5 · figures.** `Fig_picb_vs_trinity_concordance` (2×2: counts, overlap, fragmentation+length, bp) and
`Fig_picb_vs_trinity_readcapture` (read-mass capture).

## TOOLS
| Tool | What/why | Key params |
|---|---|---|
| bedtools | merge + reciprocal intersect | `merge`; `intersect -u` / `-wa -wb` |
| samtools | read-mass capture | `view -c -L <bed>`; `idxstats` (total) |
| Python (pandas/matplotlib) | aggregate + figures | Liberation Sans, vector |

## INPUTS
`_shared_data/picb_pangenome_clusters.tsv` (PICB) · `results/filter_precursors_bed/` + `all_trinity_filtred_precursors.csv.gz` (Trinity 200/200) · `results/STAR_srna_strain_wise/…/Aligned.sortedByCoord.out.bam` (read capture, PanSN).

## OUTPUTS (`figures/`, PDF+SVG+PNG + `.note.md`)
`Fig_picb_vs_trinity_concordance` · `Fig_picb_vs_trinity_readcapture` (decisive) · `Fig_pirna_length_window_test`
(TEST: 25–32 nt window) · `Fig_capture_methodology_test` (TEST: exon-blocks vs intron-span). Data:
`data/overlap_per_strain_tp.csv`, `data/read_capture_pirna.csv`, `data/SourceData_*.csv`.

## VERDICT (data-driven; triple-verified via BioMNI 2026-06-21 — genomics + general CONFIRM, literature non-committal, 0 contradict)
For defining **piRNA source loci genome-wide, PICB is the better method**: reference-anchored, calls **~2–4× more
loci**, recovers the long pachytene clusters, and captures **≥** the piRNA output (S4: ~tied at pachytene, ~1.7–2.5×
ahead at P12.5). **Trinity is a complement, not a primary caller**: annotation-free, yields precursor *transcripts*,
and at pachytene captures comparable piRNA mass with far fewer loci (precursors are few + dominant) — but it
over-segments clusters (2.6–5.8×), under-captures at P12.5, and most of its precursors recover only a small fraction
(<12%) of the full PICB cluster set. Consistent with thesis Ch.6 (Figs 6.4/6.5: Trinity over-estimates/fragments vs
the curated Zamore reference). **Use PICB for the source-locus catalogue; use Trinity to recover expressed precursor
transcripts (esp. at pachytene) where annotation is missing.** Caveats (BioMNI): both depend on sequencing depth;
both can yield false positives; validation needed.

script: /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/17_picb_vs_trinity/code/Fig_picb_vs_trinity_concordance.py
---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps:**
```bash
# reciprocal overlap (all 16 strains x 3 tp) -> overlap_per_strain_tp.csv
$PY analysis/claude_biomni_analysis/picb_vs_trinity/compute_overlap.py
# piRNA read capture (representative subset) -> read_capture.csv
$PY analysis/claude_biomni_analysis/picb_vs_trinity/read_capture.py
```

**Figure step:**
```bash
$PY figures/analysis_figures/17_picb_vs_trinity/code/Fig_picb_vs_trinity_concordance.py
$PY figures/analysis_figures/17_picb_vs_trinity/code/Fig_picb_vs_trinity_readcapture.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/17_picb_vs_trinity/code/Fig_picb_vs_trinity_concordance.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/17_picb_vs_trinity/code/Fig_picb_vs_trinity_readcapture.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/17_picb_vs_trinity/code/Fig_pirna_length_window_test.py`  _(TEST: 25–32 nt window)_
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/17_picb_vs_trinity/code/Fig_capture_methodology_test.py`  _(TEST: exon-blocks vs intron-span)_
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/17_picb_vs_trinity/code/strain_order.py`  _(imported helper)_

*Compute (`analysis/claude_biomni_analysis/picb_vs_trinity/`):*
- `…/picb_vs_trinity/compute_overlap.py` — reciprocal overlap / fragmentation / bp (Trinity EXON blocks), 16 strains × 3 tp
- `…/picb_vs_trinity/read_capture.py` — all-sRNA read-mass capture (samtools `view -c -L`; exon blocks)
- `…/picb_vs_trinity/read_capture_pirna.py` — **piRNA-specific** capture (25–32 nt; total piRNA denominator)

*Upstream pipeline (produce the inputs):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/rules/picb_cluster.smk` — PICB clusters
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/rules/trinity_assemblies.smk` + `filter_precursors.smk` + `workflow/scripts/python/filter_precursors.py` — Trinity precursors
