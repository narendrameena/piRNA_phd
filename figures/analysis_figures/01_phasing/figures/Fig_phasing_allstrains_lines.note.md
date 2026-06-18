# Fig_phasing_allstrains_lines

**piRNA phasing, all strains overlaid**

- **Shows:** Same data as small-multiples, all 16 strains overlaid as lines across timepoints.
- **How:** Fig_phasing_allstrains_lines.py. STAR 1-random-coord/read (--outSAMmultNmax 1 --outMultimapperOrder Random) on unmasked strain genome → 24–32 nt reads → GenomicRanges::follow 3′→5′ adjacency → +1-nt fraction = phasing (Almeida Genome Biol 2025, PMID 39844208).
- **Data:** phasing_allstrains_ALL_summary.csv
- **Provenance:** same run.

Full raw→figure pipeline: [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.
