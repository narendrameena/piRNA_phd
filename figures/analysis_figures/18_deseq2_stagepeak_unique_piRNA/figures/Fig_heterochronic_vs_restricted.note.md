# Fig_heterochronic_vs_restricted

**Heterochronic (stage-redeployed) vs stage-restricted within-tp unique piRNAs — TE-engaging multicopy piRNAs re-used across development, but 79 % are genuine single-locus re-expression, NOT a multimapping artifact**

- **Shows:** the within-tp genuinely-unique piRNAs split into **HETERO** = `stage-shifted (heterochronic)` (n=2,725 — unique in its home stage but the SAME sequence is expressed in another strain at a DIFFERENT stage) vs **RESTRICTED** = `conserved-but-silent` + `strain-private` (n=7,999 — not expressed at any other stage). A both TE-derived % and multimapping (NH>1) %; B TE-fraction across NH bins; C single- vs multi-locus split of heterochronic; D stage-peak strength (DESeq2 log2FC).
- **Result:** heterochronic piRNAs are **more TE-derived (45 % vs 34 %, Fisher p=4e-21)**, **more multimapping (NH>1: 21 % vs 8 %, p=2e-59)**, and **more strongly stage-peaked (log2FC median 9.4 vs 9.0, MWU p=1e-17)**. *Are they TE or multimapping?* The two are **inseparable** — TEs are multicopy, so TE-fraction rises monotonically with NH (31 %→59 %→85 %→98 %→100 %). BUT it is **NOT a pure multimapping artifact: 79 % of heterochronic piRNAs map to a SINGLE locus (NH=1)** and stay TE-enriched there (34 % vs 30 %, p=2e-3) → genuine cross-stage re-expression, with only a 21 % TE-multicopy minority. **NOT different:** length (27/30 nt), 1U, ping-pong 10A, antisense-to-TE fraction.
- **Why trustworthy:** each claim independently TESTED (Fisher exact / Mann–Whitney); multimapping = `NH` tag from the `cand_self16` self-genome alignment (primary record); classification = the adopted Theme-18 within-tp scheme (DESeq2 stage-peak, well-calibrated per `Fig_benchmark_da`). The TE↔multimapping confound is explicitly disentangled by stratifying on NH (panels B/C).
- **How:** `code/Fig_heterochronic_vs_restricted.py` on `deseq16_lenfilt/deseq_stagepeak_classified.csv.gz` + `cand_self16/*.bam` (NH) + `sense_antisense/SourceData_sense_antisense16_percand.csv.gz`.
- **Data:** `data/SourceData_Fig_heterochronic_vs_restricted.csv.gz` (per-candidate grp, NH, TE, classical, log2FC).

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
