# 09 — TE-driven piRNA-cluster evolution

**What this theme is.** Evidence that some strain-specific piRNA clusters are **created by a polymorphic TE
insertion**. Flagship: the CAST/EiJ ERVK provirus on chr5. Upstream reads/clusters = themes 07/11.

> **Status [finding · ⏳ pending BioMNI].** "A polymorphic TE insertion nucleates a new antisense-1U silencing
> cluster" is a project conclusion (`VERIFICATION_QUEUE.md`). The read-level 1U + antisense-to-TE measurements
> rest on *established* piRNA biology; the causal "TE-driven creation" claim is not yet triple-confirmed.

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1 · reads + clusters.** cutadapt 5.0 → STAR 2.7.11b (unmasked) → PICB (R 4.2.3) — themes 07/11.

**S2 · strain-private TE insertions.** **Python 3.11.15** parse the 16-strain pangenome SV VCF
(`parse_insertions16.py`; singleton carrier, ≥40 bp) → lift to strain coords.

**S3 · coordinate TE-driven test.** **bedtools 2.31.1** intersect each insertion with the carrier's PICB
clusters **in one coordinate frame** (`te_driven_coord_test.py`); a candidate must satisfy: cluster ≈
insertion, dominant family = **real TE** (LINE/SINE/LTR, not Simple_repeat), **1U ≥ 80 %**, **high
antisense-to-TE** (sense + low-1U ⇒ TE transcript, rejected). **Result #:** coordinate candidate rows CAST/EiJ
**7 487**, SPRET/EiJ **7 538** (then family/1U/antisense-filtered to high-confidence loci). This rigour caught
false positives (a (CCGT)n microsatellite; a sense/low-1U TE transcript).

**S4 · flagship validation.** **pysam (Python 3.11.15)** read-level at the CAST chr5 ERVK locus.
**Result #:** **14 312** P20.5 primary piRNAs, **~98 % 1U**, **~100 % antisense-to-ERVK**.

**S5 · figures.** matplotlib — flagship Panel A (pangenome FPM) / B (coverage + TE×strand + antisense
silencing bar) / C (base-resolution 1U); plus pangenome TE-driven test + locus gallery. → **33** figs
(incl. **19** in `figures/locus_gallery_TE/`).

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| STAR / PICB | 2.7.11b / R 4.2.3 | reads + clusters | piRNA params |
| bcftools / Python | 1.21 / 3.11.15 | parse pangenome SV insertions | singleton ≥40 bp |
| cactus / halLiftover | v2.9.3 | lift insertions/loci across strains | `halLiftover` |
| bedtools | 2.31.1 | insertion ∩ cluster (1 frame) | `intersect` |
| pysam | (Python 3.11.15) | read-level 1U / antisense-to-TE | 24–32 nt; 5′ base; TE-relative |
| Python | 3.11.15 | figures | matplotlib |

## INPUTS  pangenome SV VCF; PICB clusters; RepeatMasker → `data/TE_driven_COORDINATE_{CAST_EiJ,SPRET_EiJ}.csv`.
## OUTPUTS (`figures/`, 33)  `Fig_TEdriven_CAST_chr5_ERVK[_multi]` flagship + `Fig_te_driven_*16` + `figures/locus_gallery_TE/` (19 TE-locus views).

## DOUBLE-VERIFICATION
- Candidate filter (real-TE + 1U≥80 % + antisense) documented + shown to reject artefacts.
- Coordinate-based (not sequence-containment) test; cross-ref theme 13 (divergence, the contrasting mechanism).
- Causal "TE-driven creation" = **[finding]** (`VERIFICATION_QUEUE.md`); 1U/antisense = *established*.
