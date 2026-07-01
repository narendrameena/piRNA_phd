# Fig_pca_stagepeak

**PCA of within-tp unique-piRNA expression — wild-derived strains separate by their divergent repertoire**

- **Shows:** per stage (E16.5/P12.5/P20.5), PCA of the within-tp genuinely-unique piRNA expression across the 48 libraries, **plus a 4th panel pooling ALL timepoints** (union of unique sequences × 144 libraries, markers by timepoint); points coloured classical vs wild-derived; wild outliers labelled (arrow-leadered to clear space), the 12 classical strains annotated collectively.
- **Result:** the unique-piRNA repertoire **separates the wild-derived strains** (CAST, SPRET, PWK as outliers on their own axes) while the 12 classical strains cluster tightly at the origin. **SPRET (M. spretus, most divergent) dominates PC1** (62 % at P12.5, 83 % at P20.5). Reflects the divergence-driven unique repertoire — wild strains carry the bulk of unique piRNAs.
- **Why trustworthy:** expression streamed from the native `edger16` count matrices (only unique-sequence rows), CPM by `libsize_window`, log2; PCA = SVD on the feature-centred sample matrix (analogous to theme-07 `combine_pca16`).
- **How:** `code/Fig_pca_stagepeak.py` (coords cached in the source-data CSV).
- **Data:** `data/SourceData_Fig_pca_stagepeak.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
