# Fig_pilot_unique_pirna_counts

- **Shows:** pilot strain-private piRNA counts (C57BL_6NJ/CAST/SPRET × 3 tp)
- **Code:** `code/Fig_pilot_unique_pirna_counts.py`
- **Pipeline:** see [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.

## Why the unique-piRNA analysis includes a PCA (rationale, with numbers — verified raw→figure 2026-06-21)

The PCA (`Fig_pca_unique16` / `Fig_pca_classes16` / `Fig_pca_numbers_summary`; data `pca16/`) is the
unsupervised counterpart to the per-sequence strain-specific caller — it does four things the caller cannot:

1. **It proves strain identity structures the whole piRNA repertoire** (the premise of the entire analysis).
   16-strain, top-500-variable features: the between-strain centroid distance exceeds the within-strain
   replicate scatter by **37.6× (E16.5) / 46.4× (P12.5) / 18.6× (P20.5)** (within-replicate PC1–PC2 distance
   1.3–3.9 vs between-strain 58.5–72.7). Strain dominates piRNA variation ⇒ there IS strain-specific piRNA
   biology to dissect — the unique-piRNA calling is not chasing technical noise.
2. **It recapitulates the known phylogeny ⇒ the signal is biological, not batch/technical.** PC1 separates the
   wild-derived strains (*M. m. musculus*/*spretus*: CAST/PWK/SPRET/WSB) from the classical (*M. m. domesticus*);
   in the 3-strain pilot SPRET (*M. spretus*, most divergent) is the single most distinct cluster.
3. **Strain-distinctiveness RISES into pachytene ⇒ justifies the per-timepoint design + supports the core
   pachytene-divergence finding.** On the genuinely-unique feature set PC1 variance climbs **37.9 → 61.6 →
   76.9 %** (E16.5→P12.5→P20.5, 16-strain; n_features 50,078 / 33,502 / 23,381); pilot all-expressed 63 → 75 →
   87 %. Repertoires become MOST strain-divergent at pachytene — exactly where strain-private piRNAs accumulate.
4. **Replicate-reproducibility QC** for the ≥2/3-replicate caller: tight replicate clustering (within-rep
   1.3–3.9 ≪ 58.5–72.7 between strains) confirms no outlier/batch sample.

Verified numbers (`pca16/*.pca.csv`; `pca16/classes_pca.csv`): the pooled "Combined" PCA set = **55,442**
candidate sequences expressed in **all 3 timepoint matrices** (per-class 13,121 / 30,953 / 1,934 / 7,343 /
2,091); pooled strain-specific candidates **404,769**, of which **genuinely-unique = 106,961 (26%)**. Full
rationale: `analysis/claude_biomni_analysis/unique_pirna/METHODS_review_notes.md` (Step 7).
