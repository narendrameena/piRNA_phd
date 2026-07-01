# Fig_heterochronic_vs_restricted_exact

**Heterochronic vs stage-restricted unique piRNAs — EXACT-sequence definition (SNP-variants kept as unique) — SAME conclusion as the SNP-aware Theme-18 version**

- **Shows:** the EXACT-sequence (`klass_exact`) genuinely-unique piRNAs split into **HETERO** = `stage-shifted (heterochronic)` (n=2,588 with NH) vs **RESTRICTED** = `conserved-but-silent` + `strain-private` (n=12,393 — INCLUDING the 4,394 SNP-alleles the exact definition keeps as unique, vs 7,999 in the SNP-aware version). Same 4 panels as the Theme-18 companion: A TE-derived & multimapping %; B TE-fraction across NH; C single- vs multi-locus split; D stage-peak log2FC.
- **Result:** **identical conclusion** to the SNP-aware version — heterochronic are **more TE-derived (45 % vs 32 %, Fisher p=4e-34)**, **more multimapping (NH>1: 21 % vs 7 %, p=1e-84)**, and **more strongly stage-peaked (log2FC 9.4 vs 9.0, MWU p=1e-19)**. TE and multimapping inseparable (TE-fraction 30 %→100 % across NH). **79 % of heterochronic map to a single locus** and stay TE-enriched (34 % vs 29 %, p=2e-6) → genuine cross-stage re-expression, not a mapping artifact. The finding is **robust to including the SNP-alleles** in RESTRICTED.
- **Why trustworthy:** same method as Theme-18 `Fig_heterochronic_vs_restricted` but on `klass_exact`; each claim independently TESTED (Fisher / Mann–Whitney); multimapping = `NH` tag from the `cand_self16` self-genome alignment; TE↔multimapping confound disentangled by stratifying on NH.
- **How:** `code/Fig_heterochronic_vs_restricted_exact.py` on `data/exact_stagepeak_classified.csv.gz` + `cand_self16/*.bam` (NH) + `sense_antisense/SourceData_sense_antisense16_percand.csv.gz`.
- **Data:** `data/SourceData_Fig_heterochronic_vs_restricted_exact.csv.gz`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md). Theme-18 companion (SNP-aware): `Fig_heterochronic_vs_restricted`.
