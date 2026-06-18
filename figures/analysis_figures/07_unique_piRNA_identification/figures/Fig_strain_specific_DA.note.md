# Fig_strain_specific_DA — strain-specific piRNA counts = edgeR QL-DA (FDR<0.05, logFC>

**Shows:** strain-specific piRNA counts = edgeR QL-DA (FDR<0.05, logFC>0) ∩ presence/absence (≥2/3 reps in strain, absent elsewhere).  (unique-piRNA pipeline Step 3)
**How generated** — full unique-piRNA pipeline in [`PIPELINE.md`](../PIPELINE.md) (sRNA reads → strain-specific
 presence/absence + edgeR DA → STAR genome-anchored uniqueness test → classification → this figure).
- **Code:** `code/Fig_strain_specific_DA.py`
- **Data table (plotted points):** (see script inputs / PIPELINE.md)
- **Formats:** PDF + SVG + PNG.
