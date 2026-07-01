# Fig_unique16_class_breakdown

- **Shows:** class breakdown of unique piRNAs (route categories), 16 strains
- **Code:** `code/Fig_unique16_class_breakdown.py`
- **Pipeline:** see [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.

- **Two-level uniqueness (STAR vs pangenome):** STAR genome-anchored (mismatch ≤3) tests whether the piRNA *sequence* is expressed in other strains; **cactus halLiftover (pangenome)** tests *locus* orthology → splits genuinely-unique into **conserved-but-silent** (locus present elsewhere = expression change) vs **strain-private locus** (locus gain). The pangenome resolves the divergence-vs-absence confound (same lesson as theme 03). See PIPELINE.md S4/S6.
