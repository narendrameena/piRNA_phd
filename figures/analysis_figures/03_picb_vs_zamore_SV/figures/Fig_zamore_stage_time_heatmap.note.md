# Fig_zamore_stage_time_heatmap
- **Shows:** Developmental expression WAVE. Rows = 3 timepoint blocks × 16 strains; columns = 214 loci grouped by Zamore stage (top bar). Prepachytene loci fire at E16.5, the 99 pachytene loci fire at P20.5 — consistently across strains. Single-hue (pale=low → solid=high).
- **Method:** log10(FPM+1) of overlapping combined-run PICB clusters (pangenome projection). Stage = Zamore annotation.
- **Code:** `code/Fig_zamore_stage_time_heatmap.py` · **Data:** `data/zamore_expression_matrix_PANGENOME.csv`
- **Pipeline:** [`PIPELINE.md`](../PIPELINE.md)
