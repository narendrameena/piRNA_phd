# Fig_srccreation_NOD_ShiLtJ_LTR-ERVK_chr10_127226150 — NOD_ShiLtJ TE-creation source locus

**Carrier:** NOD_ShiLtJ  ·  own coords chr10:127226150–127231150  ·  TE family **LTR/ERVK**  ·  master FPM 1762.1  ·  class=creation
**Shows:** a strain-PRIVATE TE insertion that CREATED a piRNA cluster — expression across the 16-strain pangenome
 × 3 timepoints. The carrier expresses the cluster; the other strains lack the orthologous locus (○). Panel B uses
 the antisense-to-TE silencing second bar (length = %on-TE); Panel C = base resolution (1U, antisense=silencing).

**How generated** — full pipeline in [`../../11_locus_catalogue/PIPELINE.md`](../../11_locus_catalogue/PIPELINE.md)
(shared FASTQ→cutadapt→STAR→PICB→pangenome→RepeatMasker→fetch_primary), with one creation-specific step:
the carrier window is **snapped to its actual PICB cluster** near the TE breakpoint (gap-merged, paralog-capped) —
the master/projection coords are the one-directional TE-breakpoint window and frequently miss the adjacent cluster
(fixed 2026-06-17; see `code/make_source_pav_multi.py` "(1b) SNAP"). See this theme's `PIPELINE.md`.

- **Code:** `code/make_source_pav_multi.py` + `code/pav_clusters.py`  (driver `code/render_srccreation.sh`)
- **Command:** `make_source_pav_multi.py NOD_ShiLtJ "NOD_ShiLtJ#1#chr10" 127226150 127231150 "LTR/ERVK" "-" Fig_srccreation_NOD_ShiLtJ_LTR-ERVK_chr10_127226150 data/source_projection_creation.tsv`
- **Data:** `data/source_loci_master_creation_private.tsv` (carrier loci), `data/source_projection_creation.tsv`
  (pangenome projection), `data/picb_pangenome_clusters.tsv` (cluster table, for the window snap)
- **Raw input:** `results/STAR_srna_strain_wise/.../Aligned.sortedByCoord.out.bam`  ·  **Formats:** PDF+SVG+PNG.
