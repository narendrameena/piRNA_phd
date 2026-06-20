# 15 — piRNA biology & genomic context (16 strains)

**What these figures are.** Descriptive piRNA-biology characterisation supporting themes 07–14: biogenesis
signatures (1U + ping-pong), the developmental program, genomic-region (genic) context, pachytene-cluster
architecture, and lncRNA-derived example loci (incl. Gm10505). Not itself a strain-private/TE-driven claim.

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1 · reads.** cutadapt 5.0 → STAR 2.7.11b (unmasked strain genome, piRNA params) → per-strain piRNA reads.

**S2 · biogenesis signatures.** **pysam (Python 3.11.15)** — **1U** (5′-uridine bias) and **ping-pong**
(10-nt 5′–5′ overlap) per strain × timepoint. *Established* primary-biogenesis + slicer-amplification signatures
(used as cross-strain QC). → `Fig_biogenesis16`. **Result #:** 16 strains × 3 tp.

**S3 · developmental program.** population shift E16.5→P12.5→P20.5 (prepachytene→pachytene). → `Fig_developmental_program16`, `Fig_pachytene_cluster`.

**S4 · genomic-region (genic) context.** **bedtools 2.31.1** ∩ **Ensembl GRCm39 gene annotation** → fraction
of piRNAs from exon / intron / intergenic / lncRNA / TE, for **all** vs **unique** piRNAs and curated gene
lists. → `Fig_pirna_genic_features_*`, `Fig_pirna_genic_list3/4_*`, `Fig_pirna_genic_regions`. **Result #:**
6 genic-feature figures (list1–4).

**S5 · lncRNA-derived examples.** coverage views of lncRNA (pachytene) piRNA loci; **Gm10505** across all 16
strains; a **read-level confounding audit (exclude protein_coding)** gates any lncRNA-derived call. →
`Fig_gm10505_16strains`, `Fig_ncrna_examples`, `Fig_locus_example_ncRNA`.

**S6 · figures.** matplotlib (Python 3.11.15). **Result #:** **12** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| cutadapt / STAR | 5.0 / 2.7.11b | trim + unmasked alignment | piRNA params |
| pysam | (Python 3.11.15) | 1U / ping-pong (10-nt overlap) | 24–32 nt; 5′ base; 5′–5′ overlap |
| bedtools | 2.31.1 | genic-region overlap | `intersect` |
| Ensembl gff3/GTF | GRCm39 | gene/biotype regions | exon/intron/lncRNA/TE |
| Python | 3.11.15 | figures | matplotlib |

## INPUTS  per-strain piRNA reads; GRCm39 gene annotation; `results/.../genric_regions/*_count.csv` (featureCounts) → `data/` (genic tables).
## OUTPUTS (`figures/`, 12)  biogenesis16, developmental_program16, pachytene_cluster, gm10505_16strains, pirna_genic_* (6), ncrna_examples, locus_example_ncRNA.

## DOUBLE-VERIFICATION
- 1U / ping-pong / prepachytene→pachytene / lncRNA-derived pachytene piRNAs are all *established* — these
  figures **measure** them (characterisation/QC), not novel claims.
- lncRNA calls gated by the protein-coding confounding audit (`project_ncrna_driven_finding`).
- The divergence-driven origin of strain-private lncRNA piRNAs is the **[finding]** in theme 08, not here.

---

## SCRIPTS & COMMANDS (full paths)

Run from repo root `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA` (`export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"`; `PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python`).

**Compute steps — (re)generate the data the figures read:**
```bash
# cutadapt -> STAR (unmasked) -> PICB; 1U/ping-pong (pysam) + genic overlap (bedtools)
# are computed inside the figure scripts below.
bash workflow/scripts/run_picb_analysis_chunked.sh
```

**Figure step — render (`$PY` for .py, `Rscript` for .R, `bash` for .sh; `strain_order.py`/`pav_clusters.py` are imported helpers, not run):**
```bash
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_biogenesis16.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_developmental_program16.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_gm10505_16strains.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_locus_example_ncRNA.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_pachytene_cluster.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_pirna_genic_features_all.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_pirna_genic_regions.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/Fig_pirna_regions_alllists.py
$PY figures/analysis_figures/15_pirna_biology_misc/code/make_ncrna_examples.py
```

**All scripts (full paths):**

*Figure / analysis (`code/`):*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_biogenesis16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_developmental_program16.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_gm10505_16strains.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_locus_example_ncRNA.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_pachytene_cluster.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_pirna_genic_features_all.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_pirna_genic_regions.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/Fig_pirna_regions_alllists.py`
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/code/make_ncrna_examples.py`

*Upstream / compute:*
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/run_picb_analysis_chunked.sh` — per-replicate PICB driver (cutadapt->STAR->PICB)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_script_chunked.R` — PICB cluster calling (chunked, genome-wide LIBRARY.SIZE)
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/run_combined.sh` — combined (replicate-pooled) PICB driver
- `/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/R/picb_combine_script.R` — PICB on pooled BAM
