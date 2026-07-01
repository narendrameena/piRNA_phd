# Fig_phasing_C57BL_6NJ_timepoints

**C57BL/6NJ phasing + z-score**

- **Shows:** Panel A = +1 phasing %; Panel B = phasing z-score(+1), z=2 guide; by timepoint.
- **How:** Fig_phasing_timepoints.py. STAR 1-random-coord/read (--outSAMmultNmax 1 --outMultimapperOrder Random) on unmasked strain genome → 24–32 nt reads → GenomicRanges::follow 3′→5′ adjacency → +1-nt fraction = phasing (Almeida Genome Biol 2025, PMID 39844208).
- **Data:** phasing_C57BL_6NJ_ALL_summary.csv
- **Provenance:** C57BL/6NJ 1-random run.

Full raw→figure pipeline: [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.

- **Error bars (replicates):** ±SD across the 3 sequencing replicates per strain × timepoint; the combined/mean bar value ≈ the replicate mean. See [PIPELINE.md](../PIPELINE.md).
