# Fig_SV_mechanism

**SVs predict cluster disruption**

- **Shows:** Stacked expressed/not_expressed/not_lifted SV vs no-SV by stage; disruption scatter; per-strain SV vs disruption (P20.5).
- **How:** Fig_SV_mechanism.py. halLiftover of each locus across the 16-strain pangenome (Cactus); not_lifted = structurally disrupted; SVs from pangenome VCF; RepeatMasker TE class of each SV.
- **Data:** all_strains_*_matrix.csv
- **Provenance:** 2026-05-21 matrices.

Full raw→figure pipeline: [`PIPELINE.md`](PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.
