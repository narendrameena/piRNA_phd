# Analysis figures — piRNA cluster, phasing & SV analysis

Organized copies of the analysis figures. **Originals** stay next to their generating
scripts/data under `analysis/claude_biomni_analysis/` (and `results/`); these are curated
copies for easy access. Each figure is provided as **PDF + SVG + PNG** (vector, Liberation
Sans, ≥300 dpi). Strain order in every multi-strain figure = canonical thesis Fig 4.4 order
(median P20.5 PC1); timepoints in developmental order E16.5 → P12.5 → P20.5.

| Folder | Figure | Source script | Source data / input |
|---|---|---|---|
| **01_phasing** | `Fig_phasing_allstrains` (small-multiples, 16 strains × 3 tp) | `Fig_phasing_allstrains.py` | `phasing_allstrains_1random/ALL_summary.csv` ; source_data/SourceData_Fig_phasing_allstrains_* |
| | `Fig_phasing_allstrains_lines` (all strains overlaid) | `Fig_phasing_allstrains_lines.py` | same |
| | `Fig_phasing_perReplicate` (every replicate, horizontal) | `Fig_phasing_perReplicate.py` | source_data/SourceData_Fig_phasing_perReplicate.csv |
| | `Fig_phasing_C57BL_6NJ_timepoints` (+1% & z-score) | `Fig_phasing_timepoints.py` | phasing_C57BL_6NJ_1random/ALL_summary.csv |
| | `per_sample_qc/*.phasing.pdf` | per-sample distance distributions (143 samples, final all-strains run) | |
| **02_picb_clusters** | `Fig_picb_combined_clusters` (combined PICB clusters/strain) | `Fig_picb_combined_clusters.py` | source_data/SourceData_PICB_cluster_counts.csv |
| | `Fig_picb_rep_vs_combined` (single-rep vs combined, count) | `Fig_picb_rep_vs_combined.py` | source_data/SourceData_Fig_picb_rep_vs_combined.csv |
| | `Fig_picb_rep_combined_overlap` (overlap: *which* clusters — combined = reproducible rep clusters, <1% unique) | `Fig_picb_rep_combined_overlap.py` | source_data/SourceData_PICB_rep_combined_overlap.csv |
| **03_picb_vs_zamore_SV** ⭐REFRAMED (pangenome) | `Fig_zamore_expression_pangenome`, `Fig_zamore_retention_pangenome`, `Fig_zamore_expr_divergence_tests`, `Fig_zamore_stage_time_heatmap`, `Fig_zamore_retention_heatmap` | `Fig_zamore_*.py`, `zamore_expression_tests.py`, `pangenome_fraction_lifted.py`, `annotate_zamore_loci.py`, `rebuild_{zamore,picb}_COMBINED.py` | verified combined-PICB run + cactus halLiftover (pangenome) + GRCm39 VCF/RM/gff3 — see `03_picb_vs_zamore_SV/PIPELINE.md`. Old UCSC-liftOver figures in `superseded_UCSC_liftover/` |
| **04_picb_method_qc** | `PICB_chr_vs_wholebam_comparison`, `PICB_chrbychr_vs_wholebam`, `mapping_complexity/01-12 + report` | `picb_chr_vs_whole_bam_comparison.py`, `picb_comparison_figure.py` | results/picb_result + picb_result_combined |
| **05_C57BL_6NJ_pangenome** | `Fig1-7` (developmental classification, FPM, SV content, TE, Zamore liftover, SV expression) | `C57BL_6NJ_pangenome_figures.py` etc. | C57BL_6NJ pangenome analysis |
| **06_zamore_coverage** | `Fig1-3` (PICB cluster architecture, Zamore gene coverage, coverage detail) | `claude_biomni_*` scripts | |

## Notes
- **NEW this session:** `Fig_picb_combined_clusters`, `Fig_picb_rep_vs_combined` (+ the pending overlap figure).
- **03_picb_vs_zamore_SV** was **REFRAMED 2026-06-17**: rebuilt on the verified combined-PICB run with **pangenome (cactus halLiftover)** cross-strain projection (UCSC liftOver used only for mm10→mm39 of the Zamore loci). This overturned the old "20% SV-driven disruption" claim (a UCSC chain-sensitivity artifact → true not-projected 0.6%; loci are 98.8% retained). Old UCSC figures kept in `superseded_UCSC_liftover/`.
- **Excluded** (exploratory / superseded): `phasing_test/`, `phasing_C57BL_6NJ/` and `phasing_C57BL_6NJ_exact/` (non-1-random earlier attempts).
- Full method + verification details: `analysis/claude_biomni_analysis/source_data/README_source_data.txt`.

## 07–10 — UNIQUE piRNA & cluster-pangenome analysis (2026-06-10)
- **07_unique_piRNA_identification** — length/timepoint, strain-specific DA, genuinely-unique split (STAR), PCA.
- **08_unique_piRNA_TE_origin_and_strand** — TE families of strain-private piRNAs; sense/antisense (silencing).
- **09_TE_driven_evolution** — pangenome TE-driven test + two real locus schematics (L1, IAP).
- **10_cluster_pangenome_PAV** — open piRNA-cluster pangenome (core/accessory/private) + Zamore stage conservation.
Source scripts + full reviewer notes (Steps 0–11): analysis/claude_biomni_analysis/unique_pirna/METHODS_review_notes.md

## 11–14 — per-locus pangenome figures + circos (2026-06-17)
New themes, each self-documenting in the format: **`PIPELINE.md`** (full raw-FASTQ → figure pipeline, every
tool + param, biology grounded, novel claims tagged per `VERIFICATION_QUEUE.md`) + **`code/`** (generating
scripts) + **`data/`** (data tables; the 212 MB `picb_pangenome_clusters.tsv` lives once in `_shared_data/`
and is symlinked) + **`figures/`** (each figure as PDF+SVG+PNG with a per-figure `*.note.md`).
- **11_locus_catalogue** — 60 figs (20 loci × main/multi/single), `make_pav_locus{,_multi,_single}.py`.
- **12_creation_source_loci** — 27 figs, `make_source_pav_multi.py` (incl. snap-to-cluster window fix).
- **13_divergence_loci** — 12 figs, `make_pav_locus_multi.py` (conserved locus, strain-restricted expression).
- **14_circos_pangenome_TE** — 10 figs (TE expression sense vs piRNA silencing antisense, genome-wide).
- **09_TE_driven_evolution** — extended with the CAST chr5 ERVK flagship + `PIPELINE.md` + `code/`.

## 15 — piRNA biology & genomic context (2026-06-17)
- **15_pirna_biology_misc** — biogenesis (1U/ping-pong), developmental program, genic-feature context,
  pachytene-cluster architecture, lncRNA-derived example loci (incl. Gm10505). Characterisation/QC that
  supports themes 07–14.

---

## ✅ Organization complete — all 15 themes self-documenting (2026-06-17)

Every theme **01–15** now carries the full format: **`PIPELINE.md`** (raw FASTQ → figure, every tool +
param, biology grounded, novel claims tagged per `VERIFICATION_QUEUE.md`) + **`code/`** (generating scripts)
+ **`data/`** (data tables; the 212 MB `picb_pangenome_clusters.tsv` lives once in `_shared_data/` and is
symlinked) + per-figure **`*.note.md`**. **225 unique figures** organized, none left unplaced.

| # | Theme | Figs | Re-derived |
|---|---|---|---|
| 01 | phasing | 5 | ✓ from `Fig_phasing_*.py` |
| 02 | picb_clusters | 3 | ✓ from `Fig_picb_*.py` |
| 03 | picb_vs_zamore_SV | 5 | ✓ **REFRAMED 2026-06-17** (pangenome halLiftover; UCSC artifact corrected; +tools table; old figs in superseded/) |
| 04 | picb_method_qc | 3 | ✓ (library-size QC) |
| 05 | C57BL_6NJ_pangenome | 7 | ✓ pilot |
| 06 | zamore_coverage | 5 | ✓ (+black6 QC; external data) |
| 07 | unique_piRNA_identification | 18 | ✓ BioMNI-verified routes |
| 08 | unique_piRNA_TE_origin_and_strand | 8 | ✓ |
| 09 | TE_driven_evolution | 24 (+17 locus gallery) | ✓ finding |
| 10 | cluster_pangenome_PAV | 6 | ✓ two-mechanism (BioMNI) |
| 11 | locus_catalogue | 60 | ✓ most detailed |
| 12 | creation_source_loci | 29 | ✓ snap-fix |
| 13 | divergence_loci | 12 | ✓ |
| 14 | circos_pangenome_TE | 10 | ✓ |
| 15 | pirna_biology_misc | 12 | ✓ characterisation |

Originals stay under `analysis/claude_biomni_analysis/` (and `results/`); these are curated, documented,
double-verified copies. Biological claims are tagged *established* vs **[finding]**; pending triple-BioMNI
items are tracked in `analysis/claude_biomni_analysis/VERIFICATION_QUEUE.md`.
