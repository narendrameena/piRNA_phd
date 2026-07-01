# Fig_depth_confound_check

- **Shows:** sequencing-depth confound check for the strain-private (novel) piRNA finding. 3 panels:
  - **A** — pooled per-strain depth (Σ STAR input reads over all libraries) vs novel piRNAs (Σ over 3 timepoints), log-y, wild vs classical (one point per strain).
  - **B** — depth-normalised novel/Mread per strain (wild median 77 vs classical 2.8 = 27×, MWU p=5e-4).
  - **C** — **per-library sequencing depth** (each dot = one library = strain×timepoint×rep; — = strain mean). Wild & classical libraries span overlapping depths, so the wild excess is not a per-library depth effect. Source: `data/Fig_depth_confound_check_perlibrary.csv`.
- **A/B are per-strain (pooled across replicates+timepoints), NOT replicate-level**; C shows the per-library granularity. Novel/strain-private is inherently a pooled presence/absence property, so per-library *novel* counts are not meaningful — only depth is shown per library.
- **Verdict:** depth range (2.5×) cannot generate the novel range (290×); the gap survives depth normalisation → tracks phylogenetic divergence, not depth.
- **Code:** `code/depth_confound_check.py` · **Data:** `data/Fig_depth_confound_check.csv` + `_perlibrary.csv`
- **Pipeline:** see [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.
