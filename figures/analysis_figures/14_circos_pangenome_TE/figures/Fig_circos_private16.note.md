# Fig_circos_private16 — circos pangenome view

**Shows:** strain-private piRNA loci on the pangenome, 16 strains.  Rings = the 16 strains (canonical thesis order); axis = GRCm39 genome in 2 Mb bins;
 sense/antisense are relative to the TE strand (NOT genomic +/−).

**How generated** — established piRNA/TE biology (sense-to-TE = TE transcript/expression; antisense-to-TE =
 piRNA silencing). Per-strain primary sRNA (cutadapt→STAR, see `../11_locus_catalogue/PIPELINE.md` §1–4) is
 summed over TEs (RepeatMasker) and split SENSE vs ANTISENSE (`antisense=(TE_strand=='-')!=read.is_reverse`)
 per family per 2 Mb bin (`build_te_sense_antisense.py`), aggregated to the `data/active_*` tables; PICB
 clusters from `data/picb_pangenome_clusters.tsv`. This script draws the rings.

- **Code:** `code/Fig_circos_private16.py`
- **Data:** `data/active_te_expression_sRNA.tsv` (+_tp) [SENSE=TE expr], `data/active_pirna_on_te_sRNA_tp.tsv`
  [ANTISENSE=piRNA], `data/active_1u_bias_tp.tsv`, `data/te_sense_antisense/{strain}.tsv` (per-strain raw),
  `data/picb_pangenome_clusters.tsv`
- **Raw input:** strain sRNA BAMs + RepeatMasker TE BEDs  ·  **Formats:** PDF+SVG+PNG.
