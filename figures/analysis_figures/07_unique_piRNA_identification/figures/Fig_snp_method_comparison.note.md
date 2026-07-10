# Fig_snp_method_comparison

**SNP-variant producer method test — the reference (expressed pool vs genome), not the tool, determines correctness.**

- **Shows:** (A) recall of the delivered SNP-variant set (129S1→C57, n=330) — genomic-STAR proxy (`classify_step416`) **84.5 %** vs direct expressed-pool search by **bowtie -v3 and STAR both 100 %**; (B) why the genome proxy misses (84 % of the 51 misses don't align to the genome — genomic presence ≠ expression); (C) bowtie exhaustive (635 found) + substitution-only vs STAR non-exhaustive (631, missed 4) + 31 % gapped / 22 % reverse records.
- **Conclusion:** the delivered SNP-variant numbers (217,559 = 54 % of klass5) are **correct** — produced by direct expressed-pool matching; `classify_step416.py`'s **genome reference** was the defect (only ~50–86 % reproducible). The producer should be a bowtie1 `-v3 --norc` direct-pool search.
- **How:** `code/make_snp_method_fig.py` (self-contained). Full method, numbers, and diagnostics: `SNP_VARIANT_METHOD_TEST.md`.
- **Data:** measured on the 129S1→C57 task; C57BL/6NJ expressed pool = 34.9 M sequences; bowtie1 -v3, STAR 2.7.10a.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
