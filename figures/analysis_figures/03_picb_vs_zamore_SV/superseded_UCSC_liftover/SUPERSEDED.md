# SUPERSEDED — UCSC-liftOver theme 03 figures (2026-05-21)

These 5 figures (`Fig_PICB_vs_Zamore`, `Fig_SV_TE`, `Fig_SV_mechanism`,
`Fig_allstrains_all_timepoints`, `Fig_allstrains_SV_expression`) used **UCSC chain
liftOver** for cross-strain projection and concluded "~20% of conserved piRNA loci
are structurally disrupted, SV-driven."

That ~20% was an **artifact of UCSC liftOver's ≥95%-identity requirement**, which
scores sequence-*diverged* loci as *absent*. The pangenome (cactus halLiftover)
analysis (parent folder) shows only **0.6%** are truly not-projected and **98.8%**
are retained — consistent with pachytene loci being conserved in position but
divergent in sequence (Yu 2021 PMID 33397987).

**Kept for provenance, not deleted.** Canonical theme-03 figures are the
`Fig_zamore_*` pangenome figures in the parent folder; see `../PIPELINE.md`.
