# Fig_unique_expression_heatmap_classical

**Within-tp genuinely-unique stage-peak piRNA expression — CLASSICAL strains only (companion to the wild split)**

- **Shows:** rows = the **554** piRNAs unique to a classical strain (no labels, ordered home stage → strain → mechanism); columns = the 12 classical strains (canonical order) × E16.5/P12.5/P20.5 (tp-major). Colour = log2(CPM+1), single-hue red (white = null → dark red = high), **shared scale with the wild heatmap** for direct comparison. Left strip = mechanism; top strip = stage.
- **Result:** a clean within-tp **strain block-diagonal** — each classical strain's unique piRNAs are expressed in that strain's column at its home stage, sparse elsewhere. Classical strains contribute only **554 of 10,724 (5 %)** of all genuinely-unique piRNAs (vs 10,170 wild) — visually far sparser than the wild panel at the same colour scale.
- **Why trustworthy:** reuses the cached seq×48 expression matrix (CPM by libsize_window, mean of 3 reps); rows subset to classical-strain home; same method as `Fig_unique_expression_heatmap`.
- **How:** `code/Fig_unique_expression_heatmap_split.py` (clade = classical).
- **Data:** `data/SourceData_Fig_unique_expression_heatmap.csv.gz` (shared seq × 48 matrix).

Full pipeline: [`PIPELINE.md`](../PIPELINE.md). Companion: `Fig_unique_expression_heatmap_wild`.
