# Fig_strain_specific_coverage16

- **Shows:** how much of **total piRNA expression (read-mass)** is **strain-private**, per strain — i.e. coverage, not counts. Read-mass-weighted % of each strain's total piRNA reads (all sequences, all 3 reps) that fall in the strain-private set.
  - **A** — per-strain coverage under the absence-rule ladder (loose ≥2/3 reps / **≥2 reads = adopted** / strict 0 reads); log-y; wild = red + shaded.
  - **B** — by-group + overall summary (classical / wild / all), same ladder (+ DA-only once available).
  - **C** — DA-only + the **DA-intersection under the absence ladder**: DA∩loose (current pipeline) → DA∩≥2-read (ADOPTED) → DA∩strict (0 reads); per strain, replicate error bars + labels. **DA-only ≈ 49% of total piRNA read-mass** (classical 45%, wild 57%) — ~half of expression is strain-ENRICHED, but that is mostly relative-abundance differences, NOT specific. The **DA-intersection is a small slice and shrinks as absence tightens** (overall 1.36% loose → 1.14% ≥2-read → 0.72% strict), wild-enriched throughout. DA∩strict computed exactly by restricting the ≥2-read set to 0-reads-elsewhere (no extra edgeR run). DA-only from `edger16_coverage.R` (recomputed intersection == saved ✓); DA∩≥2-read = `strain_specific_DA_2read` (`edger16_2read.R`).
- **Code:** `code/Fig_strain_specific_coverage16.py` (+ `coverage_probe.R`, `edger16_coverage.R`) · **Data:** `data/Fig_strain_specific_coverage16.csv`; `edger16/*.coverage_probe.csv` / `*.coverage_full.csv`.
- **Pipeline:** see [`PIPELINE.md`](../PIPELINE.md).

- **Key numbers [finding]:** strain-private piRNAs are numerous but a **small slice of total piRNA read-mass** — overall **1.14%** (≥2-read; 1.36% loose, 0.72% strict). **Wild-enriched ~70×:** classical ≈0.05% vs wild ≈3.5% (≥2-read); **SPRET 9.4%** of its entire piRNA read-mass is strain-private (CAST 2.7%, PWK 1.9%, WSB 0.28%). All classical ≤0.12%.
- **≥1-rep vs ≥2-rep:** identical for coverage (filterByExpr forces 3 reps for private calls). Absence definition + BioMNI sign-off: see [[Fig_strain_specific_DA16]] note + [[project_unique_pirna_verification]].
