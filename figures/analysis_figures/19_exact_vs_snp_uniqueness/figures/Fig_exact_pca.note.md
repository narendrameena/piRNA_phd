# Fig_exact_pca

**PCA of EXACT-sequence unique piRNA expression — wild strains separate (SPRET dominant), like Theme 18**

- **Shows:** per stage (E16.5/P12.5/P20.5), PCA of the exact-sequence genuinely-unique piRNA expression across the 48 libraries (log2 CPM by libsize_window); points coloured classical vs wild-derived; wild outliers labelled, 12 classical strains annotated collectively.
- **Result:** the same structure as the Theme-18 PCA — the exact-unique repertoire separates the wild-derived strains, with **SPRET/EiJ dominating PC1** (63 % at P12.5, 82 % at P20.5); classical strains cluster tightly. Including the 4,394 SNP-alleles does NOT change the qualitative structure (wild divergence dominates).
- **Why trustworthy:** reads the cached per-rep CPM (`extract_exact_expression.py`); per-tp SVD on the feature-centred sample matrix; identical method to Theme-18 `Fig_pca_stagepeak`.
- **How:** `code/Fig_exact_pca.py` on `data/exact_cpm_perrep.csv.gz`.
- **Data:** `data/SourceData_Fig_exact_pca.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md). Theme-18 companion: `Fig_pca_stagepeak` (SNP-aware set).
