# Fig_sense_antisense_stageshift

**Sense/antisense-to-TE of the 3 unique mechanisms — the stage-shifted class is antisense-biased (silencing-competent)**

- **Shows:** A % antisense-to-TE per within-tp unique mechanism (16 strains pooled, dots = per-strain; 50% = no strand bias); B antisense fraction by top TE family × mechanism. Orientation is relative to the TE strand: antisense = silencing-competent (base-pairs the TE transcript).
- **Result:** **stage-shifted 56% antisense** and **conserved-but-silent 57%** are antisense-biased (silencing-competent), whereas **strain-private (new insertions) 47%** is ~50/50 (dual-strand). So the heterochronic stage-shifted class is **silencing-competent like conserved-but-silent**, and distinct from the dual-strand de-novo signature of brand-new insertions. Combined with the TE-origin figure, the stage-shifted class is a **hybrid**: young-TE ORIGIN (ERVK/L1) like insertion-driven, but the antisense SILENCING signature of conserved loci — i.e. shared young-TE silencing loci expressed at different developmental stages in different strains.
- **Why trustworthy:** canonical theme-08 method — piRNA strand (cand_self16) vs STRANDED TE annotation (RM `.out`, `$9=="C"→"-"`), orientation per production locus; reused from the classification-independent cache `sense_antisense/SourceData_sense_antisense16_percand.csv.gz` (DESeq2 stage-peak candidates are a subset by cand_id). Lower bound (cand_self16 = main chr+MT). **Interpretation BioMNI 3/3-verified** (new insertions dual-strand vs conserved loci antisense-biased; stage-shifted antisense ⇒ silencing-competent).
- **How:** `code/Fig_sense_antisense_stageshift.py` (join DESeq2 within-tp classes to the orientation cache).
- **Data:** `data/SourceData_Fig_sense_antisense_stageshift.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
