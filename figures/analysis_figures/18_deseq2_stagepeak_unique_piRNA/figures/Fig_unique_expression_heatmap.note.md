# Fig_unique_expression_heatmap

**Expression heatmap of ALL within-tp genuinely-unique stage-peak (27/30 nt) piRNAs across strain × timepoint**

- **Shows:** rows = all 10,724 within-tp genuinely-unique piRNA sequences (NO labels), ordered home stage → home strain (canonical) → mechanism; columns = 48 strain × timepoint (16 strains × E16.5/P12.5/P20.5, tp-major). Colour = expression log2(CPM+1) (CPM by libsize_window, mean of 3 reps), Nature-Genetics-style single-hue RED ramp (white = null → dark red = high). Left strip = mechanism (strain-private / conserved-but-silent / stage-shifted); top strip = stage; strain x-labels coloured classical/wild.
- **Result:** **within-tp uniqueness appears as a clean tp block-diagonal** — each stage's unique piRNAs are expressed in that stage's column block, and within each block the signal is a strain-diagonal. Brightness is **concentrated in the wild-derived strains** (right of each block = CAST/PWK/SPRET), confirming their dominant unique repertoire. The **stage-shifted (green strip) rows additionally light up in another strain at a DIFFERENT stage (off-diagonal)** — the heterochronic signature made visible.
- **Why trustworthy:** expression streamed directly from the native edger16 count matrices (only the unique-sequence rows), CPM by `libsize_window`, mean over the 3 replicates; cached matrix for re-rendering. 10,263 distinct sequences across the 10,724 entries.
- **How:** `code/Fig_unique_expression_heatmap.py` on `deseq16_lenfilt/deseq_stagepeak_classified.csv.gz` + `edger16/{tp}.{counts,seqs,samples}`.
- **Data:** `data/SourceData_Fig_unique_expression_heatmap.csv.gz` (seq × 48 log2-CPM matrix).

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
