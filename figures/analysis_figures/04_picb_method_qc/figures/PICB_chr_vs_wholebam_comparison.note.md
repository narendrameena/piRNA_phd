# PICB_chr_vs_wholebam_comparison

**PICB chr-by-chr vs whole-BAM (LIBRARY.SIZE)**

- **Shows:** Three PICB modes: correct chr-by-chr (genome-wide LIBRARY.SIZE), naive chr-by-chr (per-chr size, WRONG), whole-BAM; FPM must use genome-wide library size.
- **How:** picb_chr_vs_whole_bam_comparison.py. PICB run per-chr vs whole-BAM with explicit vs implicit LIBRARY.SIZE.
- **Data:** results/picb_result + picb_result_combined (xlsx)
- **Provenance:** method-QC; see memory project_picb_concurrency_bug.

Full raw→figure pipeline: [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.
