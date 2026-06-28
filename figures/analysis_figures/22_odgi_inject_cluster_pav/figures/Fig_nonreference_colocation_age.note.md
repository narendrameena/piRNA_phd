# Fig_nonreference_colocation_age — non-reference clusters are YOUNG, strain-private TE insertions

Two reference-free tests that the non-reference clusters are the *young leading edge*:

**A — Co-location** (all-vs-all strain→strain halLiftover, `14`/`15`). Each non-reference cluster lifted into every
other strain: lifts = shared insertion, no lift = strain-private.
- **916 (66%) strain-private**; **477 (34%) shared**; a small tail (~117) present in 14–15 strains but absent from
  GRCm39 (likely C57BL/6J-lineage deletions / reference gaps).
- **Sharing follows phylogeny — 62% of shared-insertion edges are within the same wild/classical group** (inherited
  by descent). Wild strains slightly more private (69% vs 64%).

**B/C — TE age** (`16`; RepeatMasker % divergence from consensus, lower = younger).
- **Non-reference cluster TEs are younger than reference/conserved: median 14.3% vs 16.4%, p = 7×10⁻⁶⁸.**
- Holds **within** the driver families — LINE/L1 (12.4 vs 13.9), LTR/ERVL-MaLR (20.0 vs 20.7), LTR/ERV1 (11.2 vs 12.8),
  all p<10⁻¹¹ — so it is genuine *insertion age*, not family composition. LTR/ERVK already youngest in both (n.s.).

**Synthesis:** conserved clusters = OLD TEs, shared across strains (ancient pachytene factories); non-reference
clusters = YOUNG TEs, strain-private, recent post-divergence LTR/ERVK + LINE-1 insertions caught by the host piRNA
response = the leading edge of piRNA-cluster evolution. Complements `Fig_nonreference_clusters` (novelty tracks the
genome-wide TE-insertion burden, ρ=0.61).
