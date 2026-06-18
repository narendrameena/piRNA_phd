# 07 — Unique piRNA identification (sequence level)
Pilot: C57BL_6NJ, CAST_EiJ, SPRET_EiJ × {E16.5,P12.5,P20.5} × 3 reps. piRNA identity = the unique sequence.
Scripts + full reviewer notes: analysis/claude_biomni_analysis/unique_pirna/METHODS_review_notes.md.

- **Fig_unique_pirna_length** — length distribution OF the unique piRNAs + timepoint structure. Mode 27 nt;
  data-driven window **24–32 nt** (captures 96.7%); developmental shift E16.5/P12.5 mode 27–28 → P20.5 mode 30.
  Source: SourceData_unique_pirna_length.csv. (`Fig_unique_pirna_length.py`)
- **Fig_strain_specific_DA** — strain-specific piRNAs = edgeR filterByExpr + QL-DA (FDR<0.05) ∩ presence/absence.
  E16.5 305k / P12.5 256k / P20.5 274k (835k total). Source: SourceData_strain_specific_DA.csv.
- **Fig_step4_classification** — STAR genome-anchored expression test splits candidates: expressed-elsewhere
  (exact/SNP-variant = NOT unique) vs genuinely-unique (conserved-but-silent / strain-private locus).
  Genuinely unique: C57 110k, CAST 143k, SPRET 202k. Source: SourceData_step4_classification.csv.
- **Fig_pca_unique** — per-timepoint PCA (DESeq2 1.42.1, thesis Fig 5.21 method); strains separate, SPRET most
  distinct (PC1 var 63→75→87%). Source: SourceData_pca_unique.csv.
