# Fig_divergence_pav_chr2_83p82_C3H_HeJ — divergence locus chr2:83082781

**Locus:** GRCm39 chr2:83082781–83087250  ·  carrier C3H_HeJ  ·  maxFPM 2059.9
**Genome-conserved 16/16** (orthologous locus present in 16 strains) but **expressed in only 2/16** (C3H_HeJ,DBA_2J).
**Shows (DIVERGENCE-driven, not insertion-driven):** a piRNA cluster whose LOCUS is present across the pangenome
 yet is silenced in most strains and active in few — i.e. cluster gain/loss by **regulatory/sequence divergence
 at a conserved locus**, contrasted with the TE-INSERTION (creation) mechanism (theme 12). Panels A–C as in the
 catalogue (per-strain × timepoint FPM + PAV; coverage + antisense-silencing second bar; base resolution).

**How generated** — identical pipeline to the catalogue, full detail in
[`../11_locus_catalogue/PIPELINE.md`](../11_locus_catalogue/PIPELINE.md) and this theme's `PIPELINE.md`.
- **Code:** `code/make_pav_locus_multi.py` + `code/pav_clusters.py`  (driver `code/render_divergence.sh`)
- **Command:** `make_pav_locus_multi.py 2 83082781 83087250 "Divergence locus chr2:83082781 — genome-conserved (16/16), strain-restricted expression (2/16)" _ Fig_divergence_pav_chr2_83p82_C3H_HeJ`
- **Data:** `data/divergence_loci_final.tsv` (the strain-restricted-but-conserved loci), `data/picb_pangenome_clusters.tsv`, `data/locus_genome_pav.tsv`
- **Raw input:** strain sRNA BAMs `results/STAR_srna_strain_wise/…`  ·  **Formats:** PDF+SVG+PNG.
