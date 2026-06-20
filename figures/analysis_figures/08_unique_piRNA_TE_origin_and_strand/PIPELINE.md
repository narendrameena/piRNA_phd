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

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# S2 TE origin (bedtools ∩ RepeatMasker):
bash analysis/claude_biomni_analysis/unique_pirna/run_coord_te16.sh
# S3 strand = sense/antisense relative to the TE (never genomic +/-):
bash analysis/claude_biomni_analysis/unique_pirna/run_sense_antisense.sh
bash analysis/claude_biomni_analysis/unique_pirna/te_sa_array.sh
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_TE_private_families.py
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_TE_private_families16.py
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_TE_timepoint_strain.py
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_concept_four_routes_lncRNA.py
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_ncrna_driven_test.py
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_ncrna_driven_test16.py
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_sense_antisense.py
$PY figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_unique_pirna_drivers.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_TE_private_families.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_TE_private_families16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_TE_timepoint_strain.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_concept_four_routes_lncRNA.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_ncrna_driven_test.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_ncrna_driven_test16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_sense_antisense.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/code/Fig_unique_pirna_drivers.py`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/run_coord_te16.sh` — S2 strain-private piRNA x TE-family overlap
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/run_sense_antisense.sh` — S3 TE-relative sense/antisense driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/te_sa_array.sh` — S3 sense/antisense SLURM array
