# Fig_PICB_vs_Zamore

**PICB clusters vs Zamore conserved loci (SV)**

- **Shows:** 28k PICB clusters vs 214 Zamore loci: expression, disruption with/without SV, liftover failure, SV-vs-disruption (PICB r=0.33, Zamore r=0.80).
- **How:** Fig_PICB_vs_Zamore.py. halLiftover of each locus across the 16-strain pangenome (Cactus); not_lifted = structurally disrupted; SVs from pangenome VCF; RepeatMasker TE class of each SV.
- **Data:** all_strains_expression_matrix.csv, all_strains_SV_matrix.csv
- **Provenance:** 2026-05-21 PICB matrices.

Full raw→figure pipeline: [`PIPELINE.md`](PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.
