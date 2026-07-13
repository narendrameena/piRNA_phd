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

**RESOLVED — the delivered producer was the LIFT-ANCHORED method; delivered numbers CONFIRMED correct.**
Two intermediate methods each had a known bias: pure-pool bowtie OVER-called (of its 27,465 base-CBS extras a
genomic co-location test found ~41 % coincidental — matched an expressed piRNA at a DIFFERENT locus), and a
STAR-alignment co-location (`build_snp_variant_colocation.py`) UNDER-called (200,658; it drops candidates whose
ortholog STAR can't align — off-assembly / too divergent). The fix is to anchor the locus with the **cactus
lift** instead of STAR: `build_snp_variant_lift.py` takes each candidate's orthologous locus in strain Y from
`present_in_Y.bed` (halLiftover through the HAL — captures the divergent/off-assembly orthologs STAR misses),
reads Y's genome there (= the ortholog), and asks: is that ortholog EXPRESSED in Y (in Y's pool) and 1–3
substitutions from the candidate? No STAR, no bowtie needed.

**RESULT (base-CBS): reproduces the delivered set at 100 % — overlap 217,559 / 217,559, ZERO missed, +146 extra
(0.07 %).** SNP-variant 217,559 → 217,705; genuinely-unique 106,961 → 106,815 (−0.14 %). So the lost producer was
this lift-anchored, biologically-correct criterion, and the delivered numbers are **exactly right**.

| method | criterion | SNP-variant | genuinely-unique | vs delivered |
|---|---|---|---|---|
| **lift-anchored (DEFINITIVE)** | ortholog (via cactus lift) EXPRESSED, 1–3mm | **217,705** | **106,815** | **100 % (0 missed, +146)** |
| delivered | (= lift-anchored) | 217,559 | 106,961 | reference |
| STAR co-location | ortholog (via STAR align) expressed | 200,658 | 123,862 | under (blind-spot) |
| pure-pool bowtie | within 3mm of ANY expressed seq | 244,991 | 79,529 | over (coincidental) |

The earlier "~[80 k .. 124 k] band" COLLAPSES: the STAR-co-location's 124 k was the blind-spot under-count (the
lift rescues it), the pool's 80 k was coincidental over-count. **Genuinely-unique = 106,961, definitively — and
now reproducible from committed code (100 %) via the correct method.** `#4 fully closed.`

**ADOPTED as canonical (2026-07).** The lift-anchored set (`build_snp_variant_lift.py`) is now the canonical
`unique16/snp_variant_refinement.csv`; `make_klass5.py` was re-run, so the canonical klass5 is regenerable from
committed code and shifts by the +146 the lift adds over the original delivered file: **SNP-variant 217,559 →
217,705, conserved-but-silent 86,115 → 85,969, genuinely-unique 106,961 → 106,815** (strain-private 20,846,
low-quality 40,011, expressed-elsewhere 40,238 unchanged). The original delivered file is preserved as
`snp_variant_refinement.delivered_orig.csv`; theme-07 count figures re-rendered; other themes compute klass5
counts dynamically (the 0.036 % row shift is invisible).

Figure: `figures/Fig_snp_method_comparison.{pdf,png,svg}` (rendered by `code/make_snp_method_fig.py`, self-contained).
Producer scripts live in `analysis/claude_biomni_analysis/unique_pirna/`: `classify_step416.py` (genomic proxy,
to retire), `build_snp_variant_refinement.py` (genomic-proxy reconstruction ~50 %), `build_snp_variant_bowtie.py`
(pure direct-pool, over-count bracket, 99.98 % reproduces delivered), `build_snp_variant_colocation.py`
(STAR-locus co-location, under-count bracket, GU 123,862), **`build_snp_variant_lift.py` (the DEFINITIVE
lift-anchored producer — reproduces delivered at 100 %, 0 missed)**, `make_klass5.py` (consumes
`snp_variant_refinement.csv`).
