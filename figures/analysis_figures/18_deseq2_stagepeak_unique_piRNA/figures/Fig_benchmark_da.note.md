# Fig_benchmark_da

**Data-based benchmark: edgeR vs DESeq2 for strain-specific piRNA DA → DESeq2 adopted**

- **Shows:** A real-data calls per stage; B permutation-null false-positives (log); C null p<0.05 vs target 0.05; D null p-value distribution (E16.5).
- **Result:** edgeR calls ~20–25 % more, but under the **permutation null** (strain labels shuffled → all calls false) edgeR gives **~10–21× more false-positives** and is **~2× anti-conservative** (null p<0.05 ≈ 0.09–0.10 vs 0.05; KS-to-uniform 0.14–0.16). **DESeq2 is well-calibrated** (≈0.05; KS ≈0.05–0.08). So edgeR's extra calls are mostly false positives; its null p-values are U-shaped (panel D), DESeq2's are flat.
- **Why trustworthy:** one-vs-rest per strain on a 300k-feature representative sample, 5 label permutations × 3 tp; both methods their own standard normalization; benchmark design **BioMNI 3/3-verified**.
- **How:** `analysis/.../benchmark_da_methods.R` (SLURM `run_bench_da.sh`) → `bench_da/{tp}.{concordance,permnull,nullpv_hist}.csv`; `code/Fig_benchmark_da.py`.
- **Data:** `data/SourceData_Fig_benchmark_da.csv`.
- **Refs (verified):** Soneson&Delorenzi 2013 (PMID 23497356); Schurch 2016 (PMID 27022035); Love 2014 DESeq2 (PMID 25516281); Robinson 2010 edgeR (PMID 19910308).

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
