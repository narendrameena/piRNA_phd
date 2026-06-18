# Fig_strain_specific_DA16_decomposition

- **Shows:** the two filters behind the strain-specific call, **separately** (companion to [[Fig_strain_specific_DA16]]). 3 stacked panels, 16 strains × 3 timepoints, log-y:
  - **A — edgeR DA only** (FDR<0.05 & logFC>0; strain X vs mean of the other 15): ~**100k–685k** per strain — most piRNAs are merely strain-**enriched**, NOT specific.
  - **B — presence/absence only** (present ≥1 read in ≥2/3 reps of X, absent in all 15 others): the small, specific set.
  - **C — strain-specific = A ∩ B** (the final set = `Fig_strain_specific_DA16`).
- **Code:** `code/Fig_strain_specific_DA16_decomposition.py` (+ DA-only parsed from `edger16/*.log` → `da_only_counts.csv`; presence-only from `code/presence_only_counts.R` → `<tp>.presence_only_counts.csv`). **Data:** `data/Fig_strain_specific_DA16_decomposition.csv`.
- **Pipeline:** see [`PIPELINE.md`](../PIPELINE.md).

- **[finding]** **Presence/absence alone ≈ the intersection** (B ≈ C; presence-only ≥ strain-specific for every strain×tp, **max diff = 123**). So **presence/absence is the specificity-defining filter** and the **edgeR DA step is nearly redundant given it** — a piRNA present in one strain and absent in the other 15 is almost always also DA-significant. DA alone (Panel A) is ~2–3 orders of magnitude larger and not specific.
- Wild-derived strains (red, shaded block) dominate the specific set (B/C); same caveats as [[Fig_strain_specific_DA16]].
