# Fig_srccreation_LP_J_LTR-ERVL_chr5_28215000 — LP_J TE-creation source locus

**Carrier:** LP_J  ·  own coords chr5:28215000–28220000  ·  TE family **LTR/ERVL**  ·  master FPM 163.3  ·  class=creation
**Shows:** a strain-PRIVATE TE insertion that CREATED a piRNA cluster — expression across the 16-strain pangenome
 × 3 timepoints. The carrier expresses the cluster; the other strains lack the orthologous locus (○). Panel B uses
 the antisense-to-TE silencing second bar (length = %on-TE); Panel C = base resolution (1U, antisense=silencing).

**How generated** — full pipeline in [`../../11_locus_catalogue/PIPELINE.md`](../../11_locus_catalogue/PIPELINE.md)
(shared FASTQ→cutadapt→STAR→PICB→pangenome→RepeatMasker→fetch_primary), with one creation-specific step:
the carrier window is **snapped to its actual PICB cluster** near the TE breakpoint (gap-merged, paralog-capped) —
the master/projection coords are the one-directional TE-breakpoint window and frequently miss the adjacent cluster
(fixed 2026-06-17; see `code/make_source_pav_multi.py` "(1b) SNAP"). See this theme's `PIPELINE.md`.

- **Code:** `code/make_source_pav_multi.py` + `code/pav_clusters.py`  (driver `code/render_srccreation.sh`)
- **Command:** `make_source_pav_multi.py LP_J "LP_J#1#chr5" 28215000 28220000 "LTR/ERVL" "+" Fig_srccreation_LP_J_LTR-ERVL_chr5_28215000 data/source_projection_creation.tsv`
- **Data:** `data/source_loci_master_creation_private.tsv` (carrier loci), `data/source_projection_creation.tsv`
  (pangenome projection), `data/picb_pangenome_clusters.tsv` (cluster table, for the window snap)
- **Raw input:** `results/STAR_srna_strain_wise/.../Aligned.sortedByCoord.out.bam`  ·  **Formats:** PDF+SVG+PNG.
