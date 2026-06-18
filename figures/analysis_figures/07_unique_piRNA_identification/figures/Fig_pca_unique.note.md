# Fig_pca_unique — per-timepoint PCA (DESeq2 size-factor norm

**Shows:** per-timepoint PCA (DESeq2 size-factor norm; thesis Fig 5.21 method) of all-expressed vs genuinely-unique piRNAs, 3 pilot strains.  (unique-piRNA pipeline Step 7)
**How generated** — full unique-piRNA pipeline in [`PIPELINE.md`](../PIPELINE.md) (sRNA reads → strain-specific
 presence/absence + edgeR DA → STAR genome-anchored uniqueness test → classification → this figure).
- **Code:** `code/Fig_pca_unique.py`
- **Data table (plotted points):** `data/SourceData_pca_unique.csv`
- **Formats:** PDF + SVG + PNG.
