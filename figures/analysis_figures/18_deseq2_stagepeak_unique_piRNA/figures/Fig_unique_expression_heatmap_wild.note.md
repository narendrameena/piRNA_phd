# Fig_unique_expression_heatmap_wild

**Within-tp genuinely-unique stage-peak piRNA expression — WILD-derived strains only (SPRET dominates)**

- **Shows:** rows = the **10,170** piRNAs unique to a wild-derived strain (no labels, ordered home stage → strain → mechanism); columns = the 4 wild strains (WSB/CAST/PWK/SPRET, canonical order) × E16.5/P12.5/P20.5 (tp-major). Colour = log2(CPM+1), single-hue red (white = null → dark red = high), **shared scale with the classical heatmap**. Left strip = mechanism; top strip = stage.
- **Result:** wild-derived strains carry **10,170 of 10,724 (95 %)** of all genuinely-unique piRNAs. Clear within-tp **strain block-diagonal**, with **SPRET/EiJ (M. spretus, most divergent) dominating** — the large dense red blocks in the SPRET columns (esp. P12.5 and P20.5). The stage-shifted (green-strip) rows show off-diagonal expression (heterochronic — same piRNA at a different stage in another wild strain).
- **Why trustworthy:** reuses the cached seq×48 expression matrix (CPM by libsize_window, mean of 3 reps); rows subset to wild-strain home; same method as `Fig_unique_expression_heatmap`.
- **How:** `code/Fig_unique_expression_heatmap_split.py` (clade = wild).
- **Data:** `data/SourceData_Fig_unique_expression_heatmap.csv.gz` (shared seq × 48 matrix).

Full pipeline: [`PIPELINE.md`](../PIPELINE.md). Companion: `Fig_unique_expression_heatmap_classical`.
