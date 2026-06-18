# Fig_srccreation_C3H_HeJ_LTR-ERVK_chr8_40687590 — C3H_HeJ TE-creation source locus

**Carrier:** C3H_HeJ  ·  own coords chr8:40687590–40692590  ·  TE family **LTR/ERVK**  ·  master FPM 775.1  ·  class=creation
**Shows:** a strain-PRIVATE TE insertion that CREATED a piRNA cluster — expression across the 16-strain pangenome
 × 3 timepoints. The carrier expresses the cluster; the other strains lack the orthologous locus (○). Panel B uses
 the antisense-to-TE silencing second bar (length = %on-TE); Panel C = base resolution (1U, antisense=silencing).

**How generated** — full pipeline in [`../../11_locus_catalogue/PIPELINE.md`](../../11_locus_catalogue/PIPELINE.md)
(shared FASTQ→cutadapt→STAR→PICB→pangenome→RepeatMasker→fetch_primary), with one creation-specific step:
the carrier window is **snapped to its actual PICB cluster** near the TE breakpoint (gap-merged, paralog-capped) —
the master/projection coords are the one-directional TE-breakpoint window and frequently miss the adjacent cluster
(fixed 2026-06-17; see `code/make_source_pav_multi.py` "(1b) SNAP"). See this theme's `PIPELINE.md`.

- **Code:** `code/make_source_pav_multi.py` + `code/pav_clusters.py`  (driver `code/render_srccreation.sh`)
- **Command:** `make_source_pav_multi.py C3H_HeJ "C3H_HeJ#1#chr8" 40687590 40692590 "LTR/ERVK" "-" Fig_srccreation_C3H_HeJ_LTR-ERVK_chr8_40687590 data/source_projection_creation.tsv`
- **Data:** `data/source_loci_master_creation_private.tsv` (carrier loci), `data/source_projection_creation.tsv`
  (pangenome projection), `data/picb_pangenome_clusters.tsv` (cluster table, for the window snap)
- **Raw input:** `results/STAR_srna_strain_wise/.../Aligned.sortedByCoord.out.bam`  ·  **Formats:** PDF+SVG+PNG.
