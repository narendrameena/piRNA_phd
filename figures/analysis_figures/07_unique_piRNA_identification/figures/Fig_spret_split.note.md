# Fig_spret_split

- **Shows:** SPRET SNP-variant vs distinct split (Step 5 mismatch cutoff)
- **Code:** `code/Fig_spret_split.py`
- **Pipeline:** see [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.

- **Sequence-level method:** STAR genome-anchored with mismatch relaxed 0→3 to detect SNP-variants at the orthologous locus (the *locus*-level test uses cactus halLiftover — see PIPELINE.md S4/S6).
