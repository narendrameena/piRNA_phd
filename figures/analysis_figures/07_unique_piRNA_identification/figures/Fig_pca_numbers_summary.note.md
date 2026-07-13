# Fig_pca_numbers_summary

- **Shows (added 2026-06-18):** a standalone infographic that explains, intuitively, how all the strain-specific piRNA candidate numbers fit together (≥2-read adopted). Built to answer the reader questions that came up about `Fig_pca_classes16` / `Fig_pca_unique16` (why bar totals ≠ per-panel n; what each class means; how "unique" is defined).
- **Panels:**
  - **A · donut** — class composition of the 404,769 pooled candidates.
  - **B · timepoint split** — each class total = the sum of its 3 per-timepoint PCA panels (E16.5+P12.5+P20.5 stacked).
  - **C · funnel** (`composition_cascade.draw_cascade`) — pooled candidates → distinct sequences (de-dup) → expressed at all 3 timepoints (= the pooled/'Combined' PCA set), **stacked by class**.
  - **D · slopegraph** — how much of each class survives into the 'Combined' PCA (pooled → all-3-tp, % retained, log y); genuinely-unique tie-in (CBS 7,343 + strain-private 2,091 = 9,434).
  - **Hierarchy tree + glossary** — the class/sub-class logic: **UNIQUE is defined by EXPRESSION** (expressed in 1 strain, absent in the other 15 — `classify_unique16.py`); the **locus** then splits it into **shared → conserved-but-silent** (unique by divergence) vs **new → strain-private** (unique by new locus + expression), via the genome-PAV locus lift (`classify_unique16_locus.py` + `present_in_{Y}.bed`). Plain-English card per class.
- **Key totals:** 404,769 pooled = NOT-unique 297,954 (74%) + genuinely-unique 106,815 (26%); all-3-tp intersection 55,442 (14%).  *(lift-anchored canonical klass5, 2026-07; was 297,808 / 106,961 under the pre-adoption delivered SNP set.)*
- **Code:** `code/Fig_pca_numbers_summary.py` (+ shared `analysis/claude_biomni_analysis/composition_cascade.py`).
- **Source data:** `data/source_data/SourceData_pca_numbers_summary.csv` (per class: pooled, E16.5/P12.5/P20.5, distinct_sequences, all3_combined, % of pooled, % survive to combined) — written live by the producer from canonical klass5; a pre-`source_data/`-convention duplicate in `data/` (stale, 217,559/86,115) was removed 2026-07.
- **Notes:** no unicode glyphs in text (Liberation Sans lacks ✓/★/subscripts) — uses "·" and words. Deterministic add_axes layout (no tight_layout).
