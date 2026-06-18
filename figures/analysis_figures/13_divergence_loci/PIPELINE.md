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
