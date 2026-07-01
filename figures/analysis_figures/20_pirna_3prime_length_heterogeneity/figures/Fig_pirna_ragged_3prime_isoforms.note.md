# Fig_pirna_ragged_3prime_isoforms

**piRNA 3′ ends are RAGGED, not cut at an exact length — imprecise PNLDC1 trimming around the per-stage window (27 nt pre-pachytene → 30 nt pachytene)**

- **Shows:** A per-tp piRNA length distribution (broad spread; the 27→30 nt developmental shift; per-tp window circled); B % of window piRNAs with a ragged-3′ isoform per tp×window; C the 3′-offset **magnitude** profile (±1/2/3 nt shorter/longer); D the developmental **MILI(27 nt)→MIWI(30 nt) trajectory** with P12.5 carrying both classes.
- **Result:** a shorter read is the **EXACT 5′ prefix** of a longer piRNA (defined 5′, ragged 3′). Raggedness **rises with development**: **39 %** (E16.5·27) → 44 % (P12.5·27) → 56 % (P12.5·30) → **68 %** (P20.5·30, pachytene/MIWI = 1.7× E16.5). Offsets are mostly **±1 nt**, decaying to ±3 nt. At P12.5 the 30 nt (MIWI) class is already raggeder than 27 nt (MILI), 44→56 %.
- **Why trustworthy:** **exact-prefix matching, zero mismatch, no SNP** (the data are an exact-sequence set with no SNP-variant class; class-exclusion identical). **Genomically validated** — 99–100 % of prefix-isoforms map to the **same 5′ genomic locus** (one molecule, ragged 3′ trim), seq-prefix raggedness ≈ 5′-anchored raggedness. Biology **BioMNI 3/3-verified** (PNLDC1 trims, HENMT1 caps; **Pnldc1-KO → longer untrimmed piRNAs**).
- **How:** `code/Fig_pirna_ragged_3prime_isoforms.py`; validation `code/verify_genomic_5prime.py`, `code/verify_nosnp_exact.py`.
- **Data:** `unique_pirna/unique16/final_classified_clean_2read.csv.gz`; window `make_stage_peak_unique.py` (E16.5→27; P12.5→27&30; P20.5→30; range 25–32 nt). SourceData `data/SourceData_Fig_pirna_ragged_3prime_isoforms.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
