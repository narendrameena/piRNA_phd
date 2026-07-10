# SNP-variant producer — method test (2026-07)

**Question.** The delivered SNP-variant class (`snp_variant_refinement.csv`, 217,559 rows = 54 % of klass5) has
no committed producer; the only committed determinant, `classify_step416.py`, reproduces it only ~50 %
(forward-only) → ~86 % (minus-strand-fixed). Is the delivered set reproducible, and by which method?

**Root cause of the gap (confirmed).** `classify_step416.py` asks the *wrong question*. It uses **STAR to align
each candidate to another strain's GENOME**, then checks whether the *genomic* sequence is expressed. The
correct criterion for "SNP-variant" is: **is the candidate within ≤3 substitutions of a piRNA EXPRESSED in
another strain** — a direct comparison to the **expressed pool**, not the genome. Diagnostics on the delivered
"missed" candidates: 100 % have their variant allele present in the current pool (not pool drift), 100 % are
valid same-length ≤3-substitution variants, but the genome proxy cannot reach them — they either don't align to
the genome, or align to a locus whose genomic sequence ≠ the expressed allele.

**Distance metric = Hamming (same-length substitutions).** SNP = substitution; verified 100 % of delivered rows
have `len(home_seq)==len(Y_allele)` and `hamming==mm`. Edit distance / indels would not reproduce the delivered
`mm`. (Length-isoform variation is a separate axis — theme 20 — not "SNP".)

## Head-to-head test

Task: which of 129S1_SvImJ's 1,161 candidates are ≤3-substitution variants of a C57BL/6NJ-expressed piRNA?
Ground truth = the 330 delivered 129S1→C57 SNP-variants (C57 holds 84 % of 129S1's delivered set). 34.9 M C57
pool sequences.

| method | reference | recall of delivered (n=330) | notes |
|---|---|---|---|
| genomic STAR (`classify_step416`) | **genome** | **279 / 330 = 84.5 %** | misses 51: **43 don't align to genome**, 8 genome≠expressed |
| **bowtie1 `-v3 --norc -a`** | expressed **pool** | **330 / 330 = 100 %** | 34.9 M reads in **3 min**; exhaustive (found 635 SNP-variants); output is directly the answer |
| STAR-to-pool | expressed **pool** | **330 / 330 = 100 %** | non-exhaustive (found 631, **missed 4** of bowtie's 635); adds **31 % gapped-indel + 22 % reverse** records that are wrong for SNP (bowtie emits neither); 0.15 % clean map rate; needs 1.4 GB read conversion + custom Hamming re-parse; no `NM` by default |

Record composition (same task): bowtie 51,836 records = 9,931 same-length (19 %) + 41,905 substring (81 %),
**0 gapped, 0 reverse**. STAR 108,234 records = 9,549 same-length (9 %) + 41,072 substring (38 %) +
**34,086 gapped (31 %) + 23,527 reverse (22 %)**. Both filter substring (length-isoforms); only STAR adds the
model-inappropriate gapped/reverse records.

Global genomic-proxy reproduction of the *full* delivered set: **50 % (forward-only) → 86 % (minus-strand-fixed)**.

## Conclusions

1. **The reference, not the tool, is what matters.** With the expressed pool as reference, BOTH bowtie and STAR
   recover the delivered set 100 % — so the delivered producer used direct pool matching, and the delivered
   numbers are the correct ones. The genome reference in `classify_step416.py` was the defect (84.5 %).
2. **Tool = bowtie1 `-v3 --norc -a`** (not STAR): exhaustive for ≤3 mismatches (STAR missed 4 by preferring
   gapped alignments); substitution-only and same-strand by construction (0 gapped, 0 reverse records — STAR
   adds 31 % gapped-indel + 22 % reverse that are wrong for SNP); right architecture for short end-to-end reads
   (STAR fights "read ≈ same-length reference", 0.15 % clean map rate); fast, smallRNA-standard, citable, emits
   mismatch positions directly. STAR is not wrong in principle (100 % recall with the pool reference) — it is the
   wrong instrument for exhaustive short k-mismatch set-matching.
3. **`classify_step416.py` should be retired to a documented genomic PROXY**; the committed producer should be a
   bowtie `-v3` direct-pool search.

## Recommended producer (next step)

One all-candidate bowtie `-v3` index (402,937 candidates), stream each of the 16 pools through it once
(`--norc -v3 -a`), per-candidate min-mm across other strains (0 = expressed-exact, 1–3 = SNP-variant), emit
`snp_variant_refinement.csv`. ~1–2 h. Then validate global reproduction of the delivered set (expected ~100 %).

Figure: `figures/snp_method_comparison.{pdf,png,svg}` (rendered by `code/make_snp_method_fig.py`, self-contained).
Producer scripts live in `analysis/claude_biomni_analysis/unique_pirna/`: `classify_step416.py` (genomic proxy,
to retire), `build_snp_variant_refinement.py` (genomic-proxy reconstruction ~50 %), `build_snp_variant_bowtie.py`
(the direct-pool bowtie producer), `make_klass5.py` (consumes `snp_variant_refinement.csv`).
