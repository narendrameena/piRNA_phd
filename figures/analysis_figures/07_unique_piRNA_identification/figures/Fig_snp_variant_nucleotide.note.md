# Fig_snp_variant_nucleotide — nucleotide-level SNP-variant view

**Shows:** nucleotide-level SNP-variant view.  (unique-piRNA pipeline Step 4/5)
**How generated** — full unique-piRNA pipeline in [`PIPELINE.md`](../PIPELINE.md) (sRNA reads → strain-specific
 presence/absence + edgeR DA → STAR genome-anchored uniqueness test → classification → this figure).
- **Code:** `code/Fig_snp_variant_nucleotide.py`
- **Data table (plotted points):** (see script inputs / PIPELINE.md)
- **Formats:** PDF + SVG + PNG.

- **Sequence-level method:** STAR genome-anchored with mismatch relaxed 0→3 to detect SNP-variants at the orthologous locus (the *locus*-level test uses cactus halLiftover — see PIPELINE.md S4/S6).
