# 08 — TE/ncRNA origin and strand of strain-private piRNAs

**What these figures are.** Where strain-private piRNAs come from (TE families; lncRNA-derived pachytene
piRNAs) and their strand relative to the target (**antisense-to-TE = silencing**). Upstream = theme 07.

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1 · strain-private piRNAs.** exactly theme 07 (cutadapt 5.0 → STAR 2.7.11b → edgeR → STAR/pangenome split).

**S2 · TE origin.** **bedtools 2.31.1** intersect each strain-private piRNA locus with **RepeatMasker**
(per-strain BED, col4 `name|class/family`); tally TE class/family. **Result #:** dominant TE classes
**LTR/ERVK (IAP), LINE/L1, ERVL-MaLR, SINE**; TE-derived fraction ≈ **31/30/16 %** (C57/CAST/SPRET; lower
bound — index = main chr+MT). → `Fig_TE_private_families[16]`.

**S3 · strand / silencing.** sense vs antisense **relative to the TE feature**
(`antisense = (TE_strand=='-') != read.is_reverse`, **never** genomic ±). **Result #:** strain-private
piRNAs **51–57 % antisense-to-TE** vs common piRNAs 46–49 % (enriched for silencing-competent antisense).
→ `Fig_sense_antisense`.

**S4 · lncRNA (pachytene) route.** strain-private piRNAs from lncRNA loci; a **read-level confounding audit
(exclude protein_coding)** is applied before calling a cluster lncRNA-derived; cross-strain locus presence by
**cactus halLiftover** (pangenome). **[finding]** strain-private lncRNA piRNAs come mostly from CONSERVED
lncRNAs by divergence (not locus gain) — mirrors the TE finding. → `Fig_ncrna_driven_test[16]`,
`Fig_concept_four_routes_lncRNA`.

**S5 · figures.** matplotlib (Python 3.11.15). → **8** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| STAR/edgeR/DESeq2 | 2.7.11b / R 4.2.3 | strain-private calling (theme 07) | see theme 07 |
| bedtools | 2.31.1 | TE-family + gene overlap | `intersect -wa -wb` |
| RepeatMasker BED | per-strain | TE class/family | col4 `name|class/family` |
| cactus / halLiftover | v2.9.3 | lncRNA-locus cross-strain presence | `halLiftover <HAL> GRCm39 <bed> <strain>` |
| Python | 3.11.15 | strand calc + figures | TE-relative antisense; matplotlib |

## INPUTS  theme-07 unique sets; `resources/repeatMasker/*`; `gene_beds/` (protein_coding + lncRNA from v3.2 GFF3) → `data/SourceData_TE_private_families16*.csv`, `SourceData_ncrna_driven_test16.csv`.
## OUTPUTS (`figures/`, 8)  `Fig_TE_private_families[16]`, `Fig_sense_antisense`, `Fig_TE_timepoint_strain`, `Fig_unique_pirna_drivers`, `Fig_ncrna_driven_test[16]`, concept figs.

## DOUBLE-VERIFICATION
- Strand is **TE-relative** (verified in the sense/antisense builder), never genomic ±.
- lncRNA calls gated by the protein-coding read-level confounding audit.
- Divergence-vs-gain origin = **[finding]** (`VERIFICATION_QUEUE.md`); 1U + antisense rest on *established* biology.
