# Fig_step4_classification — DA candidates split by expression in other strains (STAR gen

**Shows:** DA candidates split by expression in other strains (STAR genome-anchored ≤3mm): not-unique (exact/SNP) vs genuinely-unique (conserved-silent / strain-private locus).  (unique-piRNA pipeline Step 4)
**How generated** — full unique-piRNA pipeline in [`PIPELINE.md`](../PIPELINE.md) (sRNA reads → strain-specific
 presence/absence + edgeR DA → STAR genome-anchored uniqueness test → classification → this figure).
- **Code:** `code/Fig_step4_classification.py`
- **Data table (plotted points):** (see script inputs / PIPELINE.md)
- **Formats:** PDF + SVG + PNG.

- **Two-level uniqueness (STAR vs pangenome):** STAR genome-anchored (mismatch ≤3) tests whether the piRNA *sequence* is expressed in other strains; **cactus halLiftover (pangenome)** tests *locus* orthology → splits genuinely-unique into **conserved-but-silent** (locus present elsewhere = expression change) vs **strain-private locus** (locus gain). The pangenome resolves the divergence-vs-absence confound (same lesson as theme 03). See PIPELINE.md S4/S6.
