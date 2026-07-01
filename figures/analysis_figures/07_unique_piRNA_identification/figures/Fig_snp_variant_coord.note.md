# Fig_snp_variant_coord — coordinate-anchored SNP-variant: the same piRNA at its homol

**Shows:** coordinate-anchored SNP-variant: the same piRNA at its homologous locus in SPRET vs CAST (own coords), SNPs marked — a conserved piRNA that looks strain-specific.  (unique-piRNA pipeline Step 4/5)
**How generated** — full unique-piRNA pipeline in [`PIPELINE.md`](../PIPELINE.md) (sRNA reads → strain-specific
 presence/absence + edgeR DA → STAR genome-anchored uniqueness test → classification → this figure).
- **Code:** `code/Fig_snp_variant_coord.py`
- **Data table (plotted points):** (see script inputs / PIPELINE.md)
- **Formats:** PDF + SVG + PNG.

- **Sequence-level method:** STAR genome-anchored with mismatch relaxed 0→3 to detect SNP-variants at the orthologous locus (the *locus*-level test uses cactus halLiftover — see PIPELINE.md S4/S6).
