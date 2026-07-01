# Fig_nonreference_clusters — non-reference piRNA clusters (confounder-checked)

**What:** the **1,393 piRNA clusters (0.4%)** whose sequence is ABSENT from GRCm39 — caught by re-liftover
(`06`; a cluster that does not lift = absent). These are exactly the clusters a reference-anchored PAV misses.

**Panels**
- **A** per-strain counts (42–132; wild-derived toward the high end); **93% overlap a TE**.
- **B** *multimapping-corrected* expression. Non-ref clusters are **25% multimapping** (vs 0.2% reference — TE-driven);
  on all-primary FPM they look higher (12.4 vs 6.9) but on **UNIQUE reads non-ref ≈ reference (6.6 vs 6.4)** =
  genuinely expressed, NOT inflated "higher". (Fixes an earlier all-primary overclaim; cols are allFPM/uniqFPM.)
- **C** evolution: non-ref cluster count tracks the **DIRECT** genome-wide non-ref TE-insertion burden (VCF) —
  **Spearman ρ=0.61, p=0.012**; robust to multimapping (unique-read share ρ=0.51, p=0.044) and total output
  (partial r=0.50). SPRET 76k insertions ≫ C57BL_6NJ 2.3k (clean positive control).
- **D** young **LTR/ERVK + LINE-1**; BioMNI-triple-verified piRNA–TE arms race.

**Confounders checked:** multimapping (corrected ✓), total output (partial-controlled ✓), lift-artifact
(`13` flank-lift: clean-insertion fraction does not fall with divergence ✓).
**CAVEAT:** the effect is wild-vs-classical driven (within-classical n.s.; the 4 wild strains are phylogenetically
non-independent) — a group contrast, not a continuous gradient.

**Biology:** a young active TE inserts after strains diverge → host piRNA response → a new, strain-specific,
well-expressed piRNA cluster absent from the reference = the young leading edge of piRNA-cluster evolution
(Aravin 2007; Girard 2010; Gainetdinov 2015; Frazer 2011).
