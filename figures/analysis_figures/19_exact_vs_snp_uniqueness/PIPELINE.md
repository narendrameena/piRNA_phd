# 19 — EXACT-sequence vs SNP-aware definitions of strain-unique piRNAs

**What this theme is.** A head-to-head test of two ways to call a piRNA "strain-unique", on the Theme-18 DESeq2
stage-peak (27/30 nt) within-tp candidate set. Additive — Theme 18 and all others untouched. Built 2026-06-23
on user request ("unique based on exact sequence; 1–3 SNP variants kept as distinct unique piRNAs — then test
against the version where unique does not allow 1–3 SNP").

- **EXACT-sequence** (this theme): a piRNA is unique iff its EXACT sequence is not expressed in another strain at
  the same stage. 1–3 SNP variants are treated as DIFFERENT sequences → KEPT as unique (NO SNP-refinement).
- **SNP-aware** (Theme 18): a 1–3 mm variant of a same-stage-expressed allele elsewhere is a strain ALLELE of a
  shared piRNA → EXCLUDED ("SNP-variant" class).

---

## STEP-BY-STEP

**S1 · exact-only classification.** `make_exact_unique.py` = the Theme-18 within-tp scheme WITHOUT the D4
SNP-refinement step: candidates classed "SNP-variant" revert to conserved-but-silent (their exact sequence is
still strain-specific). Genuinely-unique (exact) = conserved-but-silent + strain-private + stage-shifted.
**Result: EXACT = 15,118 vs SNP-aware = 10,724 (+4,394 SNP-alleles, +41 %).** (`Fig_exact_vs_snp_uniqueness`)

**S2 · the test — are the +4,394 innovation or standing variation?** For each of the 4,394 piRNAs that EXACT
calls unique but SNP-aware rejects, count how many OTHER strains express a 1–3 mm allele (same stage) and the
mismatch distance (`Fig_snp_allele_test.py`). **Result: they are shared (±1–3 SNP) with a MEAN of 11.1 of the 15
other strains (the distribution peaks at all 15), and 88 % differ by a SINGLE SNP** → they are SNPs in
widely-shared / conserved piRNAs (STANDING GENETIC VARIATION), not strain-specific INNOVATION.
(`Fig_snp_allele_test`)

**S3 · verdict (BioMNI 3/3-verified).** A piRNA that is a 1–3 SNP variant of one expressed in other strains is
better described as a strain ALLELE of a shared piRNA, not a novel unique piRNA; the **SNP-aware definition is
biologically preferred** because it separates true piRNA INNOVATION (new loci/sequences) from STANDING GENETIC
VARIATION (SNPs in conserved piRNAs). The EXACT-sequence definition over-counts unique by 41 %, of which **29 %
is standing variation**. Use exact-sequence only if the question is "all distinct strain-specific sequences
incl. SNP alleles"; use SNP-aware (Theme 18) for piRNA innovation.

## TOOLS
| Tool | What/why |
|---|---|
| pandas | reclassify (drop D4) + join SNP-variant records | 
| BioMNI (3 agents) | interpretation (innovation vs standing variation) — 3/3 confirm SNP-aware preferred |
| matplotlib | figures (Liberation Sans, vector) |

## INPUTS
`unique_pirna/deseq16_lenfilt/deseq_stagepeak_classified.csv.gz` (Theme-18 within-tp 6-class) ·
`unique_pirna/unique16/snp_variant_refinement_withintp.csv` (per-tp 1–3 mm SNP-allele records).

## OUTPUTS (`figures/`, PDF+SVG+PNG + `.note.md`)
`Fig_exact_vs_snp_uniqueness` (counts: exact vs SNP-aware, per stage/strain/composition) ·
`Fig_snp_allele_test` (the test: SNP-alleles shared with ~11 strains, 88 % single-SNP → standing variation) ·
`Fig_exact_expression_heatmap` (Theme-18-style red heatmap of all 15,118 exact-unique, SNP-alleles flagged) ·
`Fig_exact_pca` (per-stage PCA, wild strains separate). Data: `data/exact_stagepeak_classified.csv.gz`,
`data/exact_cpm_perrep.csv.gz`, `data/SourceData_*.csv`.

## VERDICT
**SNP-aware (Theme 18) is the biologically preferred definition of strain-unique piRNAs** (BioMNI 3/3): it
captures piRNA innovation. The exact-sequence definition inflates "unique" by 41 % with strain SNP-alleles of
widely-shared piRNAs (standing variation), so it answers a different question (distinct strain-specific
sequences) and should not be reported as innovation.

## SCRIPTS
- `analysis/claude_biomni_analysis/unique_pirna/make_exact_unique.py` — exact-only classification (drop D4)
- `figures/analysis_figures/19_exact_vs_snp_uniqueness/code/Fig_exact_vs_snp_uniqueness.py`
- `figures/analysis_figures/19_exact_vs_snp_uniqueness/code/Fig_snp_allele_test.py`
- `figures/analysis_figures/19_exact_vs_snp_uniqueness/code/extract_exact_expression.py` — stream edger16 → per-rep CPM cache
- `figures/analysis_figures/19_exact_vs_snp_uniqueness/code/Fig_exact_expression_heatmap.py` — heatmap (Theme-18 style)
- `figures/analysis_figures/19_exact_vs_snp_uniqueness/code/Fig_exact_pca.py` — per-stage PCA
