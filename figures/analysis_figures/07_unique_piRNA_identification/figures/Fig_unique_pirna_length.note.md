# Fig_unique_pirna_length — length distribution of strain-specific UNIQUE piRNAs (distin

**Shows:** length distribution of strain-specific UNIQUE piRNAs (distinct seqs/length) + per-timepoint, with the data-driven FWHM window.  (unique-piRNA pipeline Step 2 (QC))
**How generated** — full unique-piRNA pipeline in [`PIPELINE.md`](../PIPELINE.md) (sRNA reads → strain-specific
 presence/absence + edgeR DA → STAR genome-anchored uniqueness test → classification → this figure).
- **Code:** `code/Fig_unique_pirna_length.py`
- **Data table (plotted points):** `data/SourceData_unique_pirna_length.csv`
- **Formats:** PDF + SVG + PNG.
