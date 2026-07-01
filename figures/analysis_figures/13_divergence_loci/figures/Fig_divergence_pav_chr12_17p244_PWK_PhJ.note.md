# Fig_divergence_pav_chr12_17p244_PWK_PhJ — divergence locus chr12:17244521

**Locus:** GRCm39 chr12:17244521–17247794  ·  carrier PWK_PhJ  ·  maxFPM 4175.6
**Genome-conserved 16/16** (orthologous locus present in 16 strains) but **expressed in only 1/16** (PWK_PhJ).
**Shows (DIVERGENCE-driven, not insertion-driven):** a piRNA cluster whose LOCUS is present across the pangenome
 yet is silenced in most strains and active in few — i.e. cluster gain/loss by **regulatory/sequence divergence
 at a conserved locus**, contrasted with the TE-INSERTION (creation) mechanism (theme 12). Panels A–C as in the
 catalogue (per-strain × timepoint FPM + PAV; coverage + antisense-silencing second bar; base resolution).

**How generated** — identical pipeline to the catalogue, full detail in
[`../11_locus_catalogue/PIPELINE.md`](../11_locus_catalogue/PIPELINE.md) and this theme's `PIPELINE.md`.
- **Code:** `code/make_pav_locus_multi.py` + `code/pav_clusters.py`  (driver `code/render_divergence.sh`)
- **Command:** `make_pav_locus_multi.py 12 17244521 17247794 "Divergence locus chr12:17244521 — genome-conserved (16/16), strain-restricted expression (1/16)" _ Fig_divergence_pav_chr12_17p244_PWK_PhJ`
- **Data:** `data/divergence_loci_final.tsv` (the strain-restricted-but-conserved loci), `data/picb_pangenome_clusters.tsv`, `data/locus_genome_pav.tsv`
- **Raw input:** strain sRNA BAMs `results/STAR_srna_strain_wise/…`  ·  **Formats:** PDF+SVG+PNG.
