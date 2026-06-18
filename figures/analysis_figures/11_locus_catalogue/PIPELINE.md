# 11 — per-locus piRNA-cluster catalogue (pangenome)

**What these figures are.** A 3-panel "anatomy" view per piRNA locus across 16 strains: (A) pangenome PICB-FPM
per strain × timepoint, (B) the carrier's own-genome coverage + TE track + antisense-to-TE silencing bar +
per-tp FPM, (C) base-resolution nucleotide ruler (1U; red 5′ arrow = antisense-to-TE silencing).

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1 · trim/align.** cutadapt 5.0 (`TGGAATTCTCGGGTGCCAAGG`, 20–36 nt) → STAR 2.7.11b (unmasked strain genome;
`--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --alignIntronMax 1
--alignEndsType EndToEnd`).

**S2 · clusters.** PICB (R 4.2.3) combined run → `clusters` sheet.

**S3 · pangenome projection.** **cactus v2.9.3 halLiftover** strain→GRCm39 → `picb_pangenome_clusters.tsv`
(every strain's clusters in a common frame + FPM). **No pairwise liftover.**

**S4 · TE + genome PAV.** **bedtools 2.31.1** ∩ **RepeatMasker** (TE family/strand) + **halLiftover** genome
presence/absence (●/○) per strain.

**S5 · read layer.** **pysam (Python 3.11.15)** `fetch_primary` → per-position 5′-end density (1U), per-tp
FPM, antisense-to-TE fraction (TE-relative strand).

**S6 · figures.** matplotlib (Python 3.11.15) — `make_pav_locus{,_multi,_single}.py` + `pav_clusters.py`
(shared helpers; `te_strand_bar` = TE-family×strand top bar + antisense silencing second bar, width = %on-TE).
**Result #:** **20** loci × (main + multi + single) = **60** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| cutadapt | 5.0 | trim | 20–36 nt, `--discard-untrimmed` |
| STAR | 2.7.11b | unmasked alignment | piRNA params |
| PICB | R 4.2.3 | clusters | combined run |
| cactus / halLiftover | v2.9.3 | pangenome projection + genome PAV | strain↔GRCm39 |
| bedtools | 2.31.1 | TE ∩ cluster | `intersect` |
| RepeatMasker BED | per-strain | TE family/strand | col4 `name|class/family`, strand |
| pysam | (Python 3.11.15) | 1U / antisense / FPM | 24–32 nt; TE-relative strand |
| Python | 3.11.15 | figures | matplotlib, Liberation Sans |

## INPUTS  `picb_pangenome_clusters.tsv` (→`_shared_data/`), per-strain BAMs, RepeatMasker, HAL; `catalogue_loci.tsv` (20 loci).
## OUTPUTS (`figures/`, 60)  `Fig_pav_locus_{locus}[_multi|_single]` (PDF+SVG+PNG + `.note.md`).

## DOUBLE-VERIFICATION
- sense/antisense is **TE-relative** (RepeatMasker strand), not genomic ± (memory `session_resume` lesson).
- Coverage flips by TE strand; pangenome FPM is the ground-truth presence (verified before any PAV claim).
- Second-bar geometry locked after read-count verification (`project_locus_figure_redesign`).
