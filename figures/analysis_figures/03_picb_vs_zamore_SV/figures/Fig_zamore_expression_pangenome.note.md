# Fig_zamore_expression_pangenome
- **Shows:** Expression of the 214 Zamore pachytene piRNA loci across 16 strains × 3 timepoints; rises into pachytene (E16.5 78% → P12.5 86% → P20.5 98%). Number atop each bar = loci expressed (of 214); Panel B = developmental program.
- **Method:** locus "expressed" if its pangenome projection overlaps the strain's combined-run PICB clusters (cactus halLiftover; Zamore loci mm10→mm39 UCSC liftOver only). Developmental "maturation heat" timepoint ramp (gold→pumpkin→brick, colourblind-safe).
- **Error bars (Panel A):** ±1 SD across the 3 biological replicates — % expressed recomputed per replicate from each replicate's PICB clusters (`code/build_replicate_pct_expressed.py` → `replicate_pct_expressed.csv`); bar height = combined run (combined ≈ replicate mean; largest spread at P12.5 where one low-depth rep undercounts; P20.5 SD≈0). Caveat: NZO_HlLtJ-12.5dpp rep3 clusters span 19/20 major chromosomes.
- **Code:** `code/Fig_zamore_expression_pangenome.py` (+ `code/build_replicate_pct_expressed.py`) · **Data:** `data/SourceData_Fig_zamore_expression_pangenome.csv`
- **Pipeline:** [`PIPELINE.md`](../PIPELINE.md)
