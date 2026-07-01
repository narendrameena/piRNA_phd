# Fig_concept_four_routes — conceptual schematic: the 4 routes to a strain-private piRNA

**Shows:** conceptual schematic: the 4 routes to a strain-private piRNA (expressed-elsewhere / SNP-variant / conserved-but-silent / strain-private locus).  (unique-piRNA pipeline concept)
**How generated** — full unique-piRNA pipeline in [`PIPELINE.md`](../PIPELINE.md) (sRNA reads → strain-specific
 presence/absence + edgeR DA → STAR genome-anchored uniqueness test → classification → this figure).
- **Code:** `code/Fig_concept_four_routes.py`
- **Data table (plotted points):** (see script inputs / PIPELINE.md)
- **Formats:** PDF + SVG + PNG.

- **Two-level uniqueness (STAR vs pangenome):** STAR genome-anchored (mismatch ≤3) tests whether the piRNA *sequence* is expressed in other strains; **cactus halLiftover (pangenome)** tests *locus* orthology → splits genuinely-unique into **conserved-but-silent** (locus present elsewhere = expression change) vs **strain-private locus** (locus gain). The pangenome resolves the divergence-vs-absence confound (same lesson as theme 03). See PIPELINE.md S4/S6.
