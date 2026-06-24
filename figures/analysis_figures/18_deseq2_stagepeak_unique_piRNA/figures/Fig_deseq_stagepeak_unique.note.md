# Fig_deseq_stagepeak_unique

**DESeq2 stage-peak (27/30 nt) WITHIN-TIMEPOINT unique piRNAs — 10,724 of 15,967 (67 %); 3 mechanisms**

- **Shows:** A per-strain within-tp genuinely-unique by stage (log); B within-tp class composition; C the three unique mechanisms per stage.
- **Result:** **genuinely-unique within-tp = 10,724 (67 %)** = conserved-but-silent 6,560 (regulatory divergence) + **stage-shifted/heterochronic 2,725** (expressed elsewhere only at a different stage) + strain-private 1,439 (insertion). NOT unique: SNP-variant 4,394 (same-stage allele) + low-quality 849. EE-same-stage = 0. Strict-sequence subset (excl. heterochronic) = 7,999. Wild-derived strains dominate (≥10–100× more).
- **Why trustworthy:** uniqueness judged **strictly within developmental stage** (a piRNA is non-unique only if another strain makes it AT THE SAME stage) — **BioMNI 3/3** (stage-specific MILI/MIWI2 vs MIWI; heterochronic divergence is real + precedented). DESeq2 candidates are a 100 % subset of the edgeR set at strain|tp|seq → determinants reused, no recompute.
- **How:** `analysis/.../classify_deseq_stagepeak.py` (reuses `ee_withintp_diag`, `present_in_{Y}.bed`, `snp_variant_refinement_withintp`, `cand_self16`) → `deseq16_lenfilt/deseq_stagepeak_classified.csv.gz`; `code/Fig_deseq_stagepeak_unique.py`.
- **Data:** `data/SourceData_Fig_deseq_stagepeak_unique_{perstrain,byclass}.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
