# Fig_phasing_allstrains

**piRNA phasing across spermatogenesis, 16 strains (small-multiples)**

- **Shows:** Per-strain +1-nt phasing % at E16.5→P12.5→P20.5; bars=mean, dots=reps; wild strains red.
- **How:** Fig_phasing_allstrains.py. STAR 1-random-coord/read (--outSAMmultNmax 1 --outMultimapperOrder Random) on unmasked strain genome → 24–32 nt reads → GenomicRanges::follow 3′→5′ adjacency → +1-nt fraction = phasing (Almeida Genome Biol 2025, PMID 39844208).
- **Data:** phasing_allstrains_ALL_summary.csv; SourceData_Fig_phasing_allstrains_*.csv
- **Provenance:** 1-random run phasing_allstrains_1random/; established phasing biology.

Full raw→figure pipeline: [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.

- **Error bars (replicates):** ±SD across the 3 sequencing replicates per strain × timepoint; the combined/mean bar value ≈ the replicate mean. See [PIPELINE.md](../PIPELINE.md).
