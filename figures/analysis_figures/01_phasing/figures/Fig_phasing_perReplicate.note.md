# Fig_phasing_perReplicate

**piRNA phasing, every replicate**

- **Shows:** Every replicate (143 samples) as a horizontal point, by strain×timepoint.
- **How:** Fig_phasing_perReplicate.py. STAR 1-random-coord/read (--outSAMmultNmax 1 --outMultimapperOrder Random) on unmasked strain genome → 24–32 nt reads → GenomicRanges::follow 3′→5′ adjacency → +1-nt fraction = phasing (Almeida Genome Biol 2025, PMID 39844208).
- **Data:** SourceData_Fig_phasing_perReplicate.csv
- **Provenance:** per-replicate.

Full raw→figure pipeline: [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.
