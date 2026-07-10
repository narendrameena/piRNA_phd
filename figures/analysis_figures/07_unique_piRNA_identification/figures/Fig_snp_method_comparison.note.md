# Fig_snp_method_comparison

**SNP-variant producer: the delivered numbers are correct, reproducible (100 %), and set by the lift-anchored method.**

- **Shows (2×2):** (A) reference determines correctness — 129S1→C57 recall (n=330): genomic-STAR proxy (`classify_step416`) **84.5 %** vs direct expressed-pool **100 %**; (B) why the genome proxy misses — 84 % of the 51 misses don't align to the genome (presence ≠ expression); (C) full-scale reproduction of the delivered SNP set (n=217,559) by each producer — genomic proxy 86 %, pure-pool bowtie 99.98 %, STAR co-location 85.4 %, **lift-anchored 100 % (0 missed)**; (D) resulting genuinely-unique — the method band collapses: pure-pool 79.5k (over-counts coincidental), STAR co-location 123.9k (under-counts, STAR blind-spot), **lift-anchored 106.8k == delivered 106,961**.
- **Conclusion:** the SNP-variant class = piRNAs 1–3 substitutions from a piRNA **expressed at the ORTHOLOGOUS locus** in another strain. `classify_step416`'s genome reference was wrong; pure-pool over-counts; STAR co-location under-counts. Anchoring the ortholog with the **cactus lift** (`build_snp_variant_lift.py`) reproduces the delivered set exactly (100 %, 0 missed) → **genuinely-unique = 106,961, definitively and reproducibly.**
- **How:** `code/make_snp_method_fig.py` (self-contained). Full method, numbers, diagnostics: `SNP_VARIANT_METHOD_TEST.md`.
- **Data:** A/B measured on 129S1→C57 (C57 pool 34.9 M seqs); C/D full-scale (all 16 strains); producers in `analysis/.../unique_pirna/build_snp_variant_{bowtie,colocation,lift}.py`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
