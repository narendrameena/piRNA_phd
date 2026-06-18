# Fig_zamore_retention_pangenome
- **Shows:** Structural RETENTION of loci across strains. A = mean span retained per strain (+ loci-count atop); B = distribution (broken axis, stacked wild/classical) — 98.8% retained; C = per-strain QC, pangenome 99.4% vs UCSC 79.8% loci lifted (UCSC under-calls divergent strains; SPRET 52 vs 210).
- **Method:** fraction_lifted = % of a locus's 100-bp windows projecting via cactus halLiftover; UCSC liftOver shown only as the (artifactual) comparison.
- **Code:** `code/Fig_zamore_retention_pangenome.py` · **Data:** `data/zamore_fraction_lifted.csv`, `data/SourceData_ucsc_vs_pangenome_lifted.csv`
- **Pipeline:** [`PIPELINE.md`](../PIPELINE.md)
