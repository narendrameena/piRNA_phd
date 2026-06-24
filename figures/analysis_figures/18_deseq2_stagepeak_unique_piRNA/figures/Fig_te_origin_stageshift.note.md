# Fig_te_origin_stageshift

**TE family of origin of the 3 within-tp unique mechanisms — the new stage-shifted class is TE-driven like insertion-driven**

- **Shows:** A TE-derived fraction per mechanism (per-strain dots); B TE-family spectrum (% of each mechanism's TE-derived loci); C stage-shifted TE-origin, classical vs wild.
- **Result:** TE-derived (lower bound): **strain-private 50 %** (721/1,439) · **stage-shifted 45 %** (1,153/2,588) · **conserved-but-silent 31 %** (2,014/6,560). Family spectrum: strain-private and **stage-shifted are both dominated by young, active TEs — LTR/ERVK (45 % / 39 %) + LINE/L1 (17 % / 18 %)** — whereas conserved-but-silent is more diverse/older (ERVK 23 %, SINE/Alu 15 %, LTR/ERVL-MaLR 14 %, L1 15 %). So **heterochronic (stage-shifted) uniqueness is TE-driven by the same young ERVK/IAP + L1 families as insertion-driven uniqueness**, not by the older/diverse families of the regulatory (conserved-but-silent) class. Stage-shifted TE-origin is comparable in classical and wild strains (slightly higher in wild).
- **Why trustworthy:** method matches theme-08 `Fig_TE_private_families16` (cand_self16 mm0 own-genome loci, PanSN `#`-prefix stripped → chrN, ∩ per-strain RepeatMasker `{X}_repeatmasker.bed`, largest-overlap family). **CAVEAT:** cand_self16 index = main chr + MT → fractions are a **lower bound**. Interpretation **BioMNI 3/3-verified** (young ERVK/L1 = main source of recent strain-variable piRNAs; heterochronic timing shift plausibly via LTR-promoter / DNA-methylation / A-MYB).
- **How:** `code/Fig_te_origin_stageshift.py` on `deseq16_lenfilt/deseq_stagepeak_classified.csv.gz` + `cand_self16/*.bam` + `resources/repeatMasker/`.
- **Data:** `data/SourceData_Fig_te_origin_stageshift.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
