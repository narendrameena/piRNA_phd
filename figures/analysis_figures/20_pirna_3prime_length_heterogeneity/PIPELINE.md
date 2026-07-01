# Theme 20 — piRNA 3′ length heterogeneity (ragged 3′ ends / imprecise trimming)

## Question / hypothesis
PIWI proteins do **not** cut piRNAs at a single exact length (27 or 30 nt). 3′-end trimming is imprecise, so a piRNA
has a **defined 5′ end** but a **variable (ragged) 3′ end** — a shorter read is the **EXACT 5′ prefix** of a longer
isoform of the same molecule. Test whether this holds per timepoint, for the stage-characteristic 27 nt and 30 nt
classes, and find the biology behind it.

## Biology (BioMNI triple-verified — 3/3 agents, genomics + literature + general; twice)
- Length is set by the **bound-PIWI footprint** + **3′→5′ trimming by PNLDC1** (Trimmer, PARN-like), terminated by
  **HENMT1** 2′-O-methylation — trimming is **not single-nt precise** → ragged 3′ ends.
- PIWI-length classes: **MILI/PIWIL2 ~26–27 nt** (pre-pachytene; fetal/perinatal, E16.5), **MIWI2/PIWIL4 ~28 nt**
  (fetal nuclear), **MIWI/PIWIL1 ~29–30 nt** (pachytene; postnatal, P20.5). Developmental **27→30 shift**.
- **30 nt class double-verified in mouse (3/3):** MIWI binds, PNLDC1 trims to the footprint, HENMT1 caps, **A-MYB/MYBL1**
  drives postnatal onset; **Pnldc1-knockout → longer, untrimmed ~31–35 nt piRNAs** — the experimental proof that
  trimming sets the length and is otherwise imprecise (= the ragged-3′ phenomenon).
- **E16.5:** MILI(PIWIL2) + MIWI2(PIWIL4); the **ping-pong** cycle (MILI–MILI + MILI→MIWI2) is active (10-nt 5′-overlap).

## Data
- `unique_pirna/unique16/final_classified_clean_2read.csv.gz` — ≥2-read, 25–32 nt, within-tp classified piRNAs
  (399,812 rows). **Exact-sequence set** (classes: conserved-but-silent / strain-private / expressed-elsewhere-exact
  — **no SNP-variant class**); SNP variants are already excluded by construction.
- Per-tp window (data-driven, `make_stage_peak_unique.py`): **E16.5 → 27 nt; P12.5 → 27 & 30 nt; P20.5 → 30 nt**.
  24 nt excluded by the 1U-peak test → piRNA range **25–32 nt**.
- `cand_self16/{strain}.cand_self16.bam` — per-strain self-genome coords (5′ position, NH); `cluster_pav/{strain}.clusters_fpm.bed`
  — pachytene clusters (FPM ± strand, timepoint); `sense_antisense/SourceData_sense_antisense16_percand.csv.gz` — TE family + orientation;
  `edger16/{tp}.{seqs,counts}` — 48-library read abundance.

## Method — ragged-3′ isoform detection (exact, no SNP)
A window-length-L piRNA has a **shorter isoform** if its `[:L−k]` prefix is in the (L−k) sequence set, and a **longer
isoform** if it is the `[:L]` prefix of an (L+k) sequence (k = 1,2,3). Matching is **exact set membership = zero
mismatch** (a SNP breaks the match). **Genomically validated** (`verify_genomic_5prime.py`): seq-prefix raggedness ≈
5′-genomic-anchored raggedness, and **99–100 %** of prefix-isoforms sit at the **same 5′ genomic locus** (one molecule,
ragged 3′ trim) — not coincidental cross-locus / SNP matches. Class-exclusion check (`verify_nosnp_exact.py`): identical.

## Findings
1. **Ragged, rising, peaking in pachytene** — 39 % (E16.5·27) → 44 % (P12.5·27) → 56 % (P12.5·30) → **68 %** (P20.5·30)
   of window piRNAs have a length-isoform; offsets mostly ±1 nt, decaying to ±3 nt. (Fig 1)
2. **Dominates the expressed pool** — read-weighted raggedness 59→86 %; ragged isoforms **more abundant** (median 189 vs
   108 reads at P20.5) — not low-count noise. (Fig 2A)
3. **Universal across 16 strains** — 34±6 % → 68±7 %; developmental rise in **8/8** strains; wild ≈ classical. (Fig 2B)
4. **Secondary (ping-pong) piRNAs are raggeder** — E16.5 10-nt overlap **z = 7.2**; ping-pong piRNAs **66 % vs 39 %**
   (Fisher p=3e-22). Slicer-set 5′ end, trimmed (ragged) 3′. (Fig 2C)
5. **TE↔genic origin flips** — TE-derived raggeder pre-pachytene (43 vs 37 %), genic/cluster raggeder in pachytene
   (70 vs 60 %) — across the MILI→MIWI transition. (Fig 2D)
6. **Pachytene-cluster-derived** — **90 %** of pachytene piRNAs are cluster-derived; cluster 67 % vs non-cluster 57 %
   (p=9e-14); the dominant **high-expression clusters are raggedest (73 %)**. (Fig 2E)

## Scripts (code/)
- `Fig_pirna_ragged_3prime_isoforms.py` → Fig 1 (core: length spread, ragged-any %, ±offset magnitude, MILI→MIWI trajectory).
- `Fig_pirna_ragged_3prime_drivers.py` → Fig 2 (drivers: expression, per-strain, ping-pong, TE/orientation, cluster, model).
- `compute_{expression_ragged,perstrain_ragged,te_orientation,pingpong_e16,cluster_origin}.py` → SourceData_*.csv.
- `verify_genomic_5prime.py` (5′-locus concordance), `verify_nosnp_exact.py` (class-exclusion robustness).
- `biomni_verify_length.py`, `biomni_verify_30nt.py` (BioMNI triple-verification of the biology).

## Key references
- Gainetdinov et al. 2018 *Mol Cell* — PIWI-footprint + PNLDC1 trimming + HENMT1 methylation (universal piRNA biogenesis).
- Ding et al. 2017 *Nat Commun*; Zhang et al. 2017 *PNAS*; Nishimura et al. 2018 — **Pnldc1-KO → longer untrimmed piRNAs**.
- Aravin et al. 2006/2007; Girard et al. 2006 — MILI/MIWI2/MIWI length classes; pre-pachytene→pachytene shift.
- Brennecke et al. 2007; Gunawardane et al. 2007 — **ping-pong** 10-nt 5′-overlap signature.
- Li et al. 2013 *Mol Cell* — **A-MYB** initiates the pachytene piRNA program.
- Lim et al. 2015; Horwich et al. 2007 — HENMT1/Hen1 2′-O-methylation of piRNA 3′ ends.
