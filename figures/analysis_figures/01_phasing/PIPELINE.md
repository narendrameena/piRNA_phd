# 01 — piRNA phasing across spermatogenesis (16 strains)

**What these figures are.** piRNA **phasing** (head-to-tail +1-nt adjacency — the signature of
Zucchini/MitoPLD-dependent phased biogenesis) per strain × timepoint, and in detail for C57BL/6NJ.
Method follows **Almeida et al. 2025 *Genome Biol* (PMID 39844208)** (cited in every script header).

---

## STEP-BY-STEP: raw data → figure (tool · version · parameters · result)

**S0 · adapter trim.** **cutadapt 5.0** — adapter `TGGAATTCTCGGGTGCCAAGG`, `-m 20 -M 36`,
`--discard-untrimmed`. → piRNA-length reads (16.5dpc.2 regenerated from raw).

**S1 · 1-random-coordinate alignment.** **STAR 2.7.11b** — `--outSAMmultNmax 1 --outMultimapperOrder Random
--runRNGseed 777` on each strain's **unmasked** genome (piRNA params as theme 11 §2). One coordinate per read
→ no multi-locus double-counting. → per-sample BAM.

**S2 · length window.** filter to **24–32 nt** (piRNA window).

**S3 · phasing score.** **R 4.2.3** (`Rsamtools` + `GenomicRanges::follow`): per read find the nearest
downstream same-strand read (3′→5′ adjacency); a **+1-nt** gap = a phased pair; `frac_plus1` = phased / all
adjacent pairs; `zscore_plus1` = z of the +1 bin. → `ALL_summary.csv` (1 row/sample).
**Result #:** n = **143** samples; mean +1-phasing rises **E16.5 19.5 % → P12.5 27.9 % → P20.5 50.5 %**.

**S4 · figures.** **matplotlib (Python 3.11.15)** — small-multiples (bars = mean, **error bar = ±SD**, dots =
replicates), overlaid lines, per-replicate, and C57BL/6NJ %+z-score. → **5** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| cutadapt | 5.0 | trim sRNA adapter | `TGGAATTCTCGGGTGCCAAGG`, 20–36 nt, `--discard-untrimmed` |
| STAR | 2.7.11b | **1-random** alignment (avoids multimapper inflation) | `--outSAMmultNmax 1 --outMultimapperOrder Random --runRNGseed 777` |
| R | 4.2.3 | phasing score | `GenomicRanges::follow` 3′→5′; 24–32 nt; +1-nt = phased |
| Python | 3.11.15 | figures | matplotlib, Liberation Sans, vector |

## INPUTS
`results/STAR_srna_*_1random/…` (1-coord BAMs) → `phasing_allstrains_1random/ALL_summary.csv`,
`phasing_C57BL_6NJ_1random/ALL_summary.csv`. Per-sample distance distributions: `figures/per_sample_qc/` (143).

## OUTPUTS (`figures/`, PDF+SVG+PNG + `.note.md`)
`Fig_phasing_allstrains` (16-strain small-multiples) · `Fig_phasing_allstrains_lines` · `Fig_phasing_perReplicate`
· `Fig_phasing_C57BL_6NJ_timepoints` (+z-score) · `Fig_phasing_C57BL_6NJ_P20.5_pachytene` (exploratory).
Source data: `data/SourceData_Fig_phasing_*.csv`.

## DOUBLE-VERIFICATION
- Result numbers recomputed from `phasing_allstrains_ALL_summary.csv` (143 samples; 19.5/27.9/50.5 %).
- 1-random controls multimapper inflation; method matches the cited Almeida protocol (*established* biology —
  phasing rises to the pachytene stage; no novel claim).
- Excluded earlier non-1-random attempts (`phasing_test/`, `phasing_C57BL_6NJ_exact/`).
