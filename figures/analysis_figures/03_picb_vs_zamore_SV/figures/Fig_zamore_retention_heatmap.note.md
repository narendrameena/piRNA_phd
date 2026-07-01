# Fig_zamore_retention_heatmap
- **Shows:** Per-locus sequence divergence across 16 strains (rows = strains, columns = 214 loci by stage). Single-hue pale (retained) → solid (divergent); divergent loci concentrate in the wild-derived strains (SPRET highest, 5.7% mean divergence).
- **Method:** divergence = 100 − % span projecting via cactus halLiftover (PowerNorm spreads small values).
- **Code:** `code/Fig_zamore_retention_heatmap.py` · **Data:** `data/zamore_fraction_lifted.csv`, `data/zamore_locus_annotation.csv`
- **Pipeline:** [`PIPELINE.md`](../PIPELINE.md)
