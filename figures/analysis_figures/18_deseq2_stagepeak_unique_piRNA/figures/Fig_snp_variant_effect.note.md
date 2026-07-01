# Fig_snp_variant_effect

**Effect of SNP-variants on the within-tp unique set — strain alleles removed by SNP-refinement**

- **Shows:** A genuinely-unique WITHOUT vs WITH SNP-refinement per stage; B mismatch-count distribution; C SNP position along the piRNA; D transition/transversion.
- **Result:** SNP-refinement removes **4,394 strain alleles = 29 % of the naive unique set** (−24/−30/−35 % at E16.5/P12.5/P20.5). These align 1–3 mm to another strain's genome where that allele is expressed at the SAME stage — i.e. strain alleles of shared piRNAs, not novel unique. **88 % differ by a single SNP** (10 % 2 mm, 2 % 3 mm); SNP positions spread across the piRNA; **transition/transversion ≈ 1.37**.
- **Why trustworthy:** SNP set = `snp_variant_refinement_withintp.csv` (per-tp: 1–3 mm allele expressed at the candidate's own stage elsewhere), restricted to the DESeq2 stage-peak SNP-variant candidates.
- **How:** `code/Fig_snp_variant_effect.py` (closest allele per candidate; home_seq vs Y_allele mismatch positions + ti/tv).
- **Data:** `data/SourceData_Fig_snp_variant_effect.csv`, `data/SourceData_Fig_snp_variant_effect_perSNP.csv.gz`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
