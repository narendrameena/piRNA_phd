# Fig_exact_vs_snp_uniqueness

**EXACT-sequence vs SNP-aware strain-unique piRNAs — exact counts 41% more (the +4,394 SNP-alleles)**

- **Shows:** A genuinely-unique counts overall + per stage, EXACT-sequence (1–3 SNP variants kept) vs SNP-aware (Theme 18, SNP-variants excluded); B per-strain (canonical, classical/wild, log); C exact-unique composition per stage (clean vs SNP-allele).
- **Result:** **EXACT = 15,118 vs SNP-aware = 10,724 (+41 %)**; the difference = **4,394 SNP-alleles** (per stage +47 % P12.5, +24 % E16.5, +35 % P20.5). The exact definition treats a 1–3 SNP variant of a shared piRNA as a separate unique piRNA; the SNP-aware definition rejects it as a strain allele. Per strain, the SNP-alleles add unique calls across both classical and wild strains.
- **Why trustworthy:** exact = Theme-18 within-tp scheme minus the D4 SNP step (SNP-variant → conserved-but-silent); same candidate set; deterministic relabel.
- **How:** `analysis/.../make_exact_unique.py` → `data/exact_stagepeak_classified.csv.gz`; `code/Fig_exact_vs_snp_uniqueness.py`. Companion test: `Fig_snp_allele_test`.
- **Data:** `data/SourceData_Fig_exact_vs_snp_uniqueness.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
