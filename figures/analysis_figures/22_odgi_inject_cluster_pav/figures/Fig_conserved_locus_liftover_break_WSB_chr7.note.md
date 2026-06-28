# Fig_conserved_locus_liftover_break_WSB_chr7 — why coordinate liftover fails (the pangenome case)

The "82-kb / 31k-FPM WSB chr7 private cluster" turned out to be a major CONSERVED cluster mis-called by halLiftover.
This figure shows why — and why a sequence-based pangenome (not coordinate liftover) is needed.

- **(A) Expression** — the cluster is expressed in **ALL 16 strains** at the homologous locus, **20,000–79,000 FPM**
  (a major pachytene-scale cluster, NOT WSB-specific). Below: minimap2 **presence (qcov ≈ 1.0) in all 16 strains +
  GRCm39 = CONSERVED**, but each genome's **aligned coordinate** spans **chr7:0.95 Mb (AKR) → 3.22 (WSB) → 6.58 (GRCm39)
  → 7.0 (PWK)** — a ~6 Mb range. halLiftover (local-synteny coordinate projection) cannot bridge that, so it failed to
  lift the WSB copy to GRCm39 and flagged it "non-reference / strain-private."
- **(B/C)** the WSB cluster detail — sRNA coverage + LTR/ERVK + gene tracks + base resolution.

**The point:** minimap2 (sequence homology, coordinate-independent) is the arbiter; a pangenome graph would normally
resolve this by sequence — except this region is a graph gap, so even the graph misses it. A cautionary, illustrative
case that coordinate-anchored methods (liftover) can mis-call a conserved, highly-expressed cluster as strain-private —
exactly the failure mode the whole pangenome/graph workflow (themes 21–24) is built to avoid.

Built with `make_conserved_locus.py` (theme 22); cross-strain alignments in theme-24 `data/wsb_conserved.tsv`.
