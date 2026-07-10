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

## Committed producer + full-scale result

`build_snp_variant_bowtie.py`: one all-candidate bowtie `-v3` index (402,937 candidates), each of the 16 pools
streamed through it once (`--norc -v3 -a`), per-candidate min-mm across OTHER strains (0 = expressed-exact
excluded, 1–3 = SNP-variant). Ran in ~15 min on 128 cores → `unique16/snp_variant_refinement.bowtie.csv`
(295,493 rows; does NOT overwrite the delivered file).

**RESULT — reproduces the delivered set at 99.98 %:** 217,526 / 217,559 recovered, only **33 missed**. This
CONFIRMS the delivered SNP-variant numbers are correct and now reproducible from committed code — resolving the
original ~50 % concern, which was purely the wrong method (genomic proxy), not a data problem. The delivered set
is a near-exact subset of the exhaustive bowtie output; `make_klass5` uses only base-CBS candidates and the
delivered set is 100 % within base-CBS.

**NUANCE — the 27,465 extras are a MIX (delivered numbers KEPT):** exhaustive pool-matching also flags 27,465
extra base-CBS candidates as 1–3mm variants (mm=2,3-enriched: 49/28/23 %). A GENOMIC CO-LOCATION test (align the
candidate to the variant strain's genome via the committed cand_to_Y BAMs; is the matched Y_allele the sequence
at the candidate's OWN orthologous locus?) splits them: **45.4 % orthologous** (real strain-SNPs at the
candidate's own locus — the delivered's NON-exhaustive search missed them; note bowtie also finds CLOSER matches
than delivered on the shared set, 86 % vs 76 % mm=1) vs **41.1 % aligns-but-different-locus** (coincidental — the
matched expressed sequence is not the ortholog) + 13.5 % no-align (ambiguous). For calibration the delivered's
OWN calls are 81.0 % orthologous by the same test — a stringent-but-imperfect approximation, not a clean anchor.

So neither set is exactly "correct": bowtie (pure pool-match) OVER-calls (~half the extras are coincidental);
delivered (217,559) is a reasonable stringent middle estimate; a STRICT genomic-locus anchor (expression AT the
orthologous locus) would give ~188 k (but under-counts by the same genome blind-spot that limited classify_step416
to 84.5 %). The true SNP/CBS boundary is thus uncertain within roughly [188 k strict-genomic, 217.6 k delivered,
245 k pure-pool] — i.e. genuinely-unique is ~106,961 with a ~±15 % method-dependent band. The delivered numbers
sit sensibly in the middle and are KEPT; the bowtie producer is the reproducibility check (99.98 %) not a
replacement. A definitive set would need a bowtie(expression) + genomic-co-location(locus) producer.

Figure: `figures/snp_method_comparison.{pdf,png,svg}` (rendered by `code/make_snp_method_fig.py`, self-contained).
Producer scripts live in `analysis/claude_biomni_analysis/unique_pirna/`: `classify_step416.py` (genomic proxy,
to retire), `build_snp_variant_refinement.py` (genomic-proxy reconstruction ~50 %), `build_snp_variant_bowtie.py`
(the direct-pool bowtie producer), `make_klass5.py` (consumes `snp_variant_refinement.csv`).
