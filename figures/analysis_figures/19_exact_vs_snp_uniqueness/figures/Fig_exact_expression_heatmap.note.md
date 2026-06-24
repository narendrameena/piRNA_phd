# Fig_exact_expression_heatmap

**EXACT-sequence genuinely-unique piRNA expression heatmap (15,118, incl. 4,394 SNP-alleles) — Theme-18 style**

- **Shows:** rows = all exact-sequence genuinely-unique piRNAs (no labels, ordered home stage → strain → category); columns = 16 strains (canonical) × E16.5/P12.5/P20.5 (tp-major); Nature-Genetics single-hue RED (white = null → dark red = high). Left strip = 4 row categories — strain-private, conserved-but-silent (clean), **SNP-allele (standing variation, magenta)**, stage-shifted; top strip = stage.
- **Result:** same structure as the Theme-18 heatmap (within-tp strain block-diagonal; wild-derived dominate, SPRET darkest), but on the larger EXACT set. The **magenta SNP-allele rows** (the +4,394 the exact definition adds vs SNP-aware) are interspersed and expressed in their home strain/stage — visually marking the standing-variation component included by exact-sequence uniqueness.
- **Why trustworthy:** expression streamed from the native edger16 matrices for the 14,514 distinct exact-unique sequences (`extract_exact_expression.py`), CPM by libsize_window, mean of 3 reps; same method/scale convention as `Fig_unique_expression_heatmap` (Theme 18).
- **How:** `code/Fig_exact_expression_heatmap.py` on `data/exact_cpm_perrep.csv.gz` + `data/exact_stagepeak_classified.csv.gz`.
- **Data:** `data/exact_cpm_perrep.csv.gz` (seq × 144 per-rep CPM).

Full pipeline: [`PIPELINE.md`](../PIPELINE.md). Theme-18 companion: `Fig_unique_expression_heatmap` (SNP-aware set).
