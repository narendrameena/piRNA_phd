# 09 — TE-driven piRNA evolution (pangenome)
Decisive test: do strain-private piRNAs arise from strain-PRIVATE TE INSERTIONS? (PLAN §4.2.1)
- **Fig_te_driven_coord** — coordinate-based test (minimap2 of strain-private insertions → own genome;
  intersect piRNA loci). Fold-enrichment over random-locus null: COMMON piRNAs 5–82× (conserved TE piRNAs
  map to all copies incl. new private insertions) ≫ strain-private 0.8–2.5×. CONCLUSION: new TE insertions
  mostly PROPAGATE conserved family piRNAs; they do NOT preferentially create novel strain-private sequences
  (those are divergence-driven). Minority TE-insertion-derived: CAST 873 / SPRET 801 unambiguous (uniquely-
  mapping). BioMNI genomics+general CONFIRM (HIGH conf); citations UNVERIFIED.
  Source: TE_driven_coord_foldenrichment.csv, TE_driven_uniquelymapping.csv.
- **Fig_locus_example** — real example schematic: SPRET-private LINE/L1 insertion (chr2) → 368 strain-private
  piRNAs, antisense-dominant + ping-pong. (example_locus.json)
- **Fig_locus_example_IAP** — real example: SPRET-private IAP/LTR-ERVK insertion (chr7) → 80 piRNAs, 100%
  antisense (clean silencing of the iconic active mouse ERV). TE-vs-piRNA arms race, strain-specific.

## Comprehensive + nucleotide-resolution (2026-06-10)
- **Fig_locus_full_L1 / Fig_locus_full_IAP** — comprehensive per-locus: cross-strain strip (strain-private) +
  TE content + strand-resolved coverage (antisense=silencing) + 1U + a REAL ping-pong pair at base
  resolution. NOTE: these supersede the earlier Fig_locus_example* — the L1 (chr2) is on the − strand, so it
  is SENSE-dominant (512k vs 311k), not antisense-dominant (earlier label corrected here). The IAP (chr7,
  + strand) is genuinely antisense-dominant (20.9k vs 2.5k = clean silencing).
- **Fig_pingpong_nucleotide** — real sense+antisense piRNA pair, 10-nt complementary 5′ overlap, 1U + 10A.
- **Fig_biogenesis** — ping-pong (10-nt) + phasing (~27-nt periodicity) measured from full read depth; 1U 84.5%.
