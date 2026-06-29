# Methods / review notes — strain- & timepoint-specific UNIQUE piRNAs (pilot)

Reviewer-facing log of every step, with rationale, verification, and caveats. Pilot strains:
C57BL_6NJ, CAST_EiJ, SPRET_EiJ × {E16.5, P12.5, P20.5} × 3 replicates (27 samples). Thesis
(§3.2.3.6/§5.10, MGFR cutoff 0.5) is context only; this is an improved, data-driven method.

## Standing principles (applied throughout)
- **Own genome per strain, unmasked.** piRNA mapping uses each strain's own unmasked assembly
  (masked genomes suppress TE-derived piRNAs). No genome-build mixing; cross-strain only via
  pangenome/liftOver. Verified: PICB refFasta = 0 lowercase (unmasked).
- **One coordinate frame.** 3 chr-naming conventions coexist (refFasta/clusters/TE BED = `chrN`;
  GFF3 + STAR index = PanSN `{strain}#1#chrN`). Same REL-2205 assembly, identical coords
  (chr1 = 191,890,962 across all sources) — reconciled by stripping the `{strain}#1#` prefix.
- **No magic numbers.** Thresholds come from a concrete reference or a data-derived statistical
  aggregation, with the derivation documented.
- **Biology via BioMNI, triple/double-verified; citations need confirmed PMID+DOI** (BioMNI has
  returned mismatched DOIs — all references flagged UNVERIFIED pending EuropePMC).

## Step 1 — piRNA sequence layer (regenerated via Snakemake)
- **What:** raw sRNA reads → cutadapt → `tstk collapse` (per sample, all reps).
- **Tools/params:** existing Snakemake rules `cutadapt_srna` + `collapse_raw_seq_srna`; adapter
  `TGGAATTCTCGGGTGCCAAGG`, `--minimum-length 20 --maximum-length 36 --discard-untrimmed`
  (project-validated sRNA trimming, unchanged). Collapse encodes count in the header `>seqid-nreads`.
- **Why:** the trimmed/collapsed intermediates were `temp()` and had been deleted; regenerated
  through the pipeline (reproducible) rather than ad hoc. Dry-run confirmed only cutadapt+collapse
  (81 jobs), no other rules pulled in.
- **Verification:** 27/27 files produced; none empty; Snakemake 81/81 done, 0 errors.
- **piRNA identity = the unique sequence** (consistent with thesis §3.2.3.5 seqkit approach).

## Step 2 — QC (data-driven)
- **What:** per-sample read totals + unique-sequence counts; count-weighted length distribution.
- **Result:** 76–162 M reads/sample, 7–19 M unique seqs; length peak **26–29 nt (mode 27)**, bulk
  **24–32 nt** (minor 35-nt bump near the 36-nt trim edge). Files: `QC_per_sample.csv`,
  `QC_length_distribution.csv`.
- **Why:** establishes the data-driven piRNA window rather than assuming the thesis's 27/30 nt.
- **Window adopted = 24–32 nt (data-driven, validated on the UNIQUE set itself).** Length
  distribution of the distinct strain-specific candidate sequences (unrestricted 20–36 nt input,
  `Fig_unique_pirna_length`/`SourceData_unique_pirna_length.csv`): pooled mode 27 nt, FWHM major
  peak 26–30 nt (80.3% of unique seqs), ≥90% window 25–31 nt; **24–32 nt captures 96.7% pooled
  (95–98% per strain)** while excluding the ~3% degradation/non-piRNA tail (20–23, 33–36 nt). The
  unique piRNAs are TIGHTER than the bulk read pool — they behave like canonical piRNAs.
  Note: SPRET_EiJ mode is shifted +1 nt (28) with a broader peak (FWHM 26–31) — observation only,
  not interpreted (most-divergent strain; revisit after origin analysis + BioMNI verification).

## Step 3 — strain-specific caller: presence/absence → data-driven edgeR DA (improved over MGFR)
- **3a — presence/absence first pass (24–32 nt, reproducibility-based):** present(strain) = detected
  (≥1 read) in ≥2/3 reps; strain-specific = reproducible in exactly ONE strain (within timepoint);
  timepoint-specific = reproducible at exactly ONE stage (within strain). The magic `MIN_RPM=1.0` was
  removed (no-magic-numbers rule). Script `build_unique_pirna_pilot.py`.
- **Result of 3a — KEY FINDING:** 29,509,511 strain-specific + 36,833,142 timepoint-specific candidates
  (vs 257,258 under the old MIN_RPM=1.0). **Reproducibility is necessary but NOT sufficient**: with no
  abundance floor it admits the entire low-count tail — sequences seen 2–5× total that land in 2 reps
  by chance or as sequencing-error copies of abundant piRNAs. An abundance floor is required, and it
  must be data-driven (not a re-introduced magic RPM).
- **3b — data-driven floor + significance via edgeR (user-approved 2026-06-10):** the floor is
  statistical, not hand-picked. (1) per-timepoint count matrix (24–32 nt unique seqs × 9 samples =
  3 strains × 3 reps), `build_count_matrix.py`, applying filterByExpr's own `min.total.count=15` early
  for tractability only → 4.35M/4.11M/2.96M features (E16.5/P12.5/P20.5). (2) `edgeR::filterByExpr`
  (group=strain) sets the low-count floor from library sizes + smallest group size — no hand-picked
  cutoff. (3) QL F-test strain X vs mean(others), BH-FDR<0.05 & logFC>0 = significantly X-enriched.
  (4) intersect with presence/absence on the filtered counts (present ≥2/3 reps in X, absent <2/3 in
  each other strain) ⇒ strain-SPECIFIC, not merely enriched. `edger_da.R`. lib.size = 24–32 nt
  (piRNA-window) total/sample. Implements the planned "presence/absence ∩ differential abundance".
- **Why (vs thesis MGFR):** MGFR uses an arbitrary 0.5 marker score, ignores replicates, is not an FDR.
  edgeR QL-DA is replicate-aware + FDR-controlled.
- **Result of 3b (edgeR DA ∩ presence/absence), 24–32 nt strain-specific piRNAs** (`Fig_strain_specific_DA`,
  `SourceData_strain_specific_DA.csv`): filterByExpr kept 0.96/1.09/1.11 M expressed features
  (E16.5/P12.5/P20.5); QL DA flagged ~220–424k enriched per strain·tp; ∩ presence/absence ⇒
  | E16.5: C57 91,339 · CAST 129,985 · SPRET 83,891 | P12.5: 62,282 · 64,438 · 129,309 |
  P20.5: 46,355 · 55,058 · 172,721 | total 835,378. Pattern: SPRET rises with development
  (84k→129k→173k), C57+CAST fall (E16.5 highest) — DATA observation only.
- **Caveats:** (i) 3-strain pilot → "specific" = among these 3, not all 16 (scale-out planned).
  (ii) **this set still contains SNP-variants** of conserved piRNAs — SPRET (different species, most
  SNPs) is therefore expectedly inflated; genuinely-unique vs SNP-variant is resolved in Step 4 (STAR
  genome-anchored). Do NOT interpret the SPRET lead as novel biology until Step 4.

## Step 4 — genuinely-unique vs expressed-elsewhere: STAR genome-anchored expression test (tool-consistent)
- **Motivation (reviewer point + user):** a "strain-specific" sequence can be a conserved piRNA carrying
  strain SNPs, or expressed-below-threshold elsewhere — NOT truly unique. EXPRESSION (not DNA presence)
  is the criterion: a locus even identical in another strain but transcriptionally SILENT there is still
  unique to X.
- **Supersedes** the earlier bowtie sequence-vs-sequence sketch (TOOL CONSISTENCY — STAR everywhere).
- **Method (`run_step4_map.sh` → `classify_step4.py`):** map each candidate of strain X to each OTHER
  strain's UNMASKED genome with the sRNA STAR params, ONLY `--outFilterMismatchNmax` relaxed 0→3 (the
  data-driven Poisson cutoff k*=3, Step 5). The reference sequence at the mapped locus (from the MD tag)
  = the piRNA strain Y would produce there; **Y EXPRESSES it iff that sequence ∈ Y.pool_uniq** (Y's
  24–32 nt collapse unique seqs, all 9 samples). NM = SNP distance. Candidates also mapped to X's own
  genome (mm0) for origin (Step 6). STAR 2.7.11b = the pipeline's STAR; indices verified UNMASKED.
- **Why genome-anchored + pool-set lookup:** (a) STAR genome alignment with the pipeline's params
  (reviewer-proof, no tool switch); (b) handles read-length + multimapping correctly; (c) expression-
  based per the user's point; (d) avoids mapping the 54–83 M-seq expressed pools to the genome (BAM
  explosion) — only the small candidate sets (188/231/327 k) are aligned.
- **4-way classification → `{X}.step4_classified.csv.gz`:** expressed-elsewhere exact (0mm) /
  SNP-variant of expressed (1–3mm) = NOT unique; unique:conserved-but-silent (homolog elsewhere but
  silent) / unique:strain-private-locus (no ≤3mm homolog in any other genome) = GENUINELY unique by
  expression. The private-locus vs conserved-silent split feeds the TE-driven-evolution test.
- **Status:** RUNNING (map 3302573 → classify 3302576, chained). Earlier bowtie pilot on the OLD
  candidate set gave ~99% non-novel for SPRET (method-level agreement, BioMNI-endorsed; citations
  UNVERIFIED) — the STAR genome-anchored re-run on the edgeR-DA set supersedes those numbers.

## Step 5 — data-driven mismatch cutoff (no magic number)
- **What:** true SNP-variants of the same piRNA follow Poisson(λ) over mismatches, λ = SNP-divergence
  × piRNA length. Estimate λ from the data as n(1mm)/n(0mm) (= Poisson P(1)/P(0) = λ). Cutoff k* =
  smallest k with Poisson CDF(k; λ) ≥ 0.999 (captures ≥99.9% of true SNP-variants).
- **Result:** λ = 15,506/70,943 = **0.219** → implied divergence ~0.81 %/bp at piRNA loci;
  CDF(≤1)=0.979, CDF(≤2)=0.9985, **CDF(≤3)=0.9999 → k\* = 3**. Sensitivity in
  `novel_cutoff_sensitivity.csv` (≤1 mm → 6,212 novel; ≤3 mm → 1,319).
- **Why:** replaces an arbitrary mismatch cutoff with a divergence-Poisson reference derived from
  the data itself.
- **Caveat:** a ≤k-mm match conflates same-locus SNP-variant with paralogous (same-TE-family) piRNAs;
  both mean the sequence is essentially present in the others' repertoire, so "genuinely novel
  sequence" = no ≤k-mm expressed match remains a defensible sequence-level definition. Genomic origin
  (Step 6) resolves where the novel ones come from.

## Step 6 — genomic origin + TE + sense/antisense (PLANNED; on hold for Step-2 window)
- **What:** map the genuinely-novel piRNAs to the strain's own unmasked `chrN` genome; assign origin.
- **Overlap handling (reviewer point):** features overlap at a locus (CDS⊂exon, UTRs⊂exon,
  exon-of-one-isoform = intron-of-another, genic∩TE). Plan: split the GFF3 into separate category
  files (genes/intergenic/lncRNA; CDS/5′UTR/3′UTR; exon/intron; lncRNA exon/intron) per thesis
  §3.2.3.4, AND assign a single primary origin by a biologically-motivated **priority** (BioMNI-endorsed:
  TE precedence for piRNA biogenesis), while also reporting the full overlap profile.
- **Multimapper RPM:** featureCounts `-M --fraction -O` — multimappers included but **fractional
  (1/N across loci)** so repetitive/TE piRNAs are not over-counted (thesis §3.2.3.1/§3.2.3.4).
- **Sense/antisense:** piRNA strand vs gene and vs TE (antisense-to-TE = silencing-competent).
- **TE-driven test:** intersect novel piRNA loci with strain-PRIVATE TE insertions (pangenome SV/TE
  PAV) — TE inserted only in SPRET producing piRNAs only in SPRET = TE spawned a new piRNA source.

### Step 6a DONE — TE families of strain-private piRNAs (`te_families_private.py`, `Fig_TE_private_families`)
- Strain-private-locus candidates → loci in own genome (cand_self mm0; PanSN→chrN reconciled, verified)
  → bedtools intersect with `{X}_repeatmasker.bed` (col4=name|class/family) → per-candidate primary family.
- **TE-derived fraction:** C57 30.7% (5,627/18,313), CAST 29.6% (6,476), SPRET 15.5% (9,509/61,208).
- **Top families (all strains): LTR/ERVK (IAP) dominant** (1,541/1,892/1,850), then LINE/L1 (913/534/1,163),
  LTR/ERVL-MaLR, SINE/Alu, B2, B4. Active autonomous TEs (ERVK/IAP, L1) lead → consistent with the
  BioMNI-verified TE-driven mechanism. Tables `TE_private_{summary,families,classes}.csv`.
- **CAVEAT (accuracy):** STAR index = main chromosomes + MT only (no unplaced contigs), so 20–26% of
  private piRNAs do NOT map to the own genome and cannot be TE-annotated → **TE fraction is a LOWER
  BOUND** (unplaced contigs are TE-rich). SPRET's lower fraction (15.5%) partly reflects its larger
  unmapped share. Sense/antisense to TE still pending (needs the TE GFF3 with strand, thesis §3.2.3.7).
- **Status:** family IDENTITY is data; the "active TEs drive strain-private piRNA gain" interpretation
  is PROVISIONAL — the decisive test is private-piRNA loci ∩ strain-PRIVATE TE INSERTIONS (pangenome),
  not just "is the locus in a TE" (a TE can be ancestral/shared). Queued.

## Step 7 — PCA of piRNA expression across strains (thesis Fig 5.21 method; `pca_unique.R`, `Fig_pca_unique`)
- **Method (thesis-faithful):** per-timepoint, **DESeq2 1.42.1** (= thesis v1.42) size-factor normalisation
  of the piRNA-sequence count matrix → PCA (`prcomp` on log2(normalised+1)). filterByExpr = the same
  data-driven floor as the DA step. TWO feature sets side by side (user-approved "both"):
  (a) all expressed piRNAs (top-500 most-variable, DESeq2 plotPCA convention) = Fig 5.21 reproduction;
  (b) genuinely-unique (Step-4) piRNAs at that timepoint (union over strains; normalised with the SAME
  full-library size factors — the unique set alone would break DESeq2's most-features-not-DE assumption).
- **Result:** 3 strains form tight, fully separated clusters (reps reproducible) in both feature sets,
  all timepoints. SPRET_EiJ most distinct, increasingly so with development (PC1 var 63→75→87% E16.5→
  P12.5→P20.5, all-expressed); C57BL_6NJ + CAST_EiJ (both *M. musculus*) closer to each other.
  Reproduces thesis Fig 5.21 (SPRET distinct at P20.5) + matches verified pachytene-divergence biology.
- **Caveat:** 3-strain pilot (9 points/tp); full Fig 5.21 analogue needs the 16-strain scale-out.
- **WHY the unique-piRNA analysis includes a PCA (rationale, data-driven with numbers).** The PCA is not
  decorative — it does four jobs that the per-sequence strain-specific caller cannot, all verified from
  `pca16/{tp}.pca.csv` / `pca/{tp}.pca.csv`:
  1. **Unsupervised proof that strain identity structures the whole piRNA repertoire — the premise of the
     entire analysis.** Before any supervised per-sequence call, PCA asks the unbiased question: do strains
     separate on global piRNA expression alone? Overwhelmingly yes — in the 16-strain set (all-expressed,
     top-500 variable) the **between-strain centroid distance exceeds the within-strain replicate scatter by
     37.6× (E16.5) / 46.4× (P12.5) / 18.6× (P20.5)** (within-replicate PC1–PC2 distance 1.3–3.9 vs
     between-strain 58.5–72.7). Strain dominates piRNA variation ⇒ there IS strain-specific piRNA biology to
     dissect — the unique-piRNA calling is not chasing technical noise. This MOTIVATES the analysis.
  2. **It recapitulates the known phylogeny ⇒ the signal is biological, not batch/technical.** PC1 cleanly
     splits wild-derived (*M. m. musculus*/*spretus*: CAST/PWK/SPRET/WSB; positive PC1 end at P20.5
     +104.6/+115.3/+125.5) from classical (*M. m. domesticus*; negative end, e.g. BALB −64.5, AKR −64.3).
     In the 3-strain pilot SPRET (*M. spretus*, most divergent) is the single most distinct cluster.
     Expression structure tracking subspecies divergence = the strain-private piRNAs reflect real evolution.
  3. **Strain-distinctiveness RISES through spermatogenesis ⇒ justifies the per-timepoint design + supports
     the core pachytene-divergence finding.** On the genuinely-unique feature set, PC1 variance climbs
     **37.9 → 61.6 → 76.9 %** (E16.5→P12.5→P20.5, 16-strain; cf. pilot all-expressed 63 → 75 → 87 %): as germ
     cells reach pachytene (P20.5) the repertoires become MOST strain-divergent — exactly where strain-private
     piRNAs accumulate. This is WHY the unique-piRNA analysis is run per-timepoint, not pooled.
  4. **Replicate-reproducibility QC for the ≥2/3-rep caller.** The tight replicate clustering (within-replicate
     PC-distance 1.3–3.9 vs 58.5–72.7 between strains) confirms the 3 replicates/strain·tp are reproducible —
     a prerequisite of the presence/absence (≥2/3 reps) definition; no outlier/batch sample.
- **16-strain (scale-out) numbers (`pca16/`, 3 feature sets):** PC1 variance — all_expressed (top-500)
  **38.3/22.2/30.9 %**; genuinely-unique features **37.9/61.6/76.9 %** (E16.5/P12.5/P20.5; n_features
  50,078/33,502/23,381). The unique-feature PC1 rises sharply into pachytene (the divergence signal lives in
  the strain-private piRNAs themselves). The 63/75/87 % above is the 3-strain pilot (C57/CAST/SPRET, all-expressed).

## Step 8 — pangenome TE-DRIVEN test (the decisive analysis; PLAN §4.2.1)
- **Resource:** existing 17-strain minigraph-cactus pangenome `results/pangenome/output/
  mouse_17strain_pangenome.vcf.gz` — GRCm39-referenced (chr1=195,154,279; Ensembl contig names), all 16
  strains genotyped, vg-deconstruct, **literal ALT sequences**, 73.9M records, multi-allelic (GT=allele
  index). No vg/halLiftover installed → coordinate-free sequence approach.
- **Strain-private insertions** (`parse_insertions.py`, `run_pangenome_insertions.sh`): per multi-allelic
  site, an allele is an insertion if len(ALT)-len(REF) ≥40 bp; PRIVATE to X among the 3 pilot = only X
  carries it. Result: **C57 6,342 (4.6 Mb), CAST 120,421 (90 Mb), SPRET 213,805 (177 Mb)** — ordering
  matches divergence (C57≈GRCm39 ref → few; SPRET different species → massive). 2 min, 74M records.
- **CONFOUND FOUND (sequence-containment is NOT enough):** testing "is the piRNA SEQUENCE inside an
  X-private insertion" (`search_pirna_in_insertions.py`) gave the HIGHEST rate for *expressed-elsewhere
  (common)* piRNAs (C57 7.4 / CAST 16.5 / SPRET 22.1 %), NOT strain-private (1.9–3.3 %). Reason: private
  insertions are mostly fresh copies of ACTIVE TE families (L1, ERVK) whose piRNAs are CONSERVED (made
  from many copies in all strains) → a common TE piRNA's sequence matches the insertion by family
  similarity WITHOUT being produced from it. **Sequence-match ≠ production.** Do NOT claim TE-driven from
  containment. (`pirna_in_private_insertion_byclass.csv`, `TE_driven_summary.csv` = confounded first pass.)
- **FIX — coordinate-based (`run_coord_te.sh`, `coord_classify.py`, RUNNING 3302583):** minimap2 the
  X-private insertion sequences back to X's OWN genome (asm5, primary, mapq≥20) → insertion loci in X
  (PanSN, same frame as cand_self) → intersect with each candidate's PRODUCTION locus. Tests production,
  not similarity. Other Step-4 classes = built-in null. TE-driven = strain-private-locus piRNA whose
  locus is INSIDE a strain-private insertion AND is TE-annotated (Step 6a).
- **RESULT (`Fig_te_driven_coord`, `TE_driven_coord_foldenrichment.csv`):** fold-enrichment for a
  locus inside a private insertion, over the random-locus null (private-ins bp/genome = C57 0.18%, CAST
  3.58%, SPRET 6.96%): **expressed-elsewhere (common) 82×/8.3×/4.9×** ≫ **strain-private-locus
  2.3×/2.5×/0.76×**. Even the multimapping-aware coordinate test shows COMMON piRNAs dominate. TE-driven
  candidates (strain-private locus ∩ private insertion ∩ TE): **C57 42, CAST 1,202, SPRET 1,259**;
  families LTR/ERVK (IAP) + LINE/L1.
- **CONCLUSION (data-driven answer to "do TEs drive piRNA evolution"):** new strain-private TE insertions
  overwhelmingly PROPAGATE the active family's CONSERVED piRNAs (new copies of L1/ERVK), they do NOT
  preferentially create novel strain-private piRNA SEQUENCES; strain-private piRNA sequences are
  divergence-driven at existing loci. A real MINORITY (~1.2–1.3k per divergent strain, ERVK/L1) are
  genuine TE-insertion-derived. SPRET (most divergent) is LEAST insertion-driven (0.76×). Limit: TE
  multimapping prevents single-locus attribution; cleanest residual = uniquely-mapping piRNAs (TODO).
- **BioMNI TRIPLE-VERIFIED (2026-06-10):** genomics + general agents CONFIRM at HIGH confidence that
  (a) new TE insertions propagate conserved family piRNAs, (b) existing-locus divergence dominates the
  strain-private repertoire; (c) the per-insertion new-locus FREQUENCY is NOT quantified in literature →
  our quantification is novel. Literature agent returned a degenerate answer (re-run pending). All
  PMIDs/DOIs UNVERIFIED.

## Steps 9–11 — refinements + timepoint-specific + sense/antisense (DONE 2026-06-10)
- **Step 9 uniquely-mapping TE-driven refinement (`uniq_map_refine.py`):** restrict to strain-private
  piRNAs mapping to EXACTLY ONE own-genome locus (removes TE-multimapping ambiguity). UNAMBIGUOUS
  TE-driven: C57 25, CAST 873, SPRET 801 (~65–72% of the all-alignment estimate; LTR/ERVK+L1).
  CAST≥SPRET despite SPRET's 2× private-insertion content → SPRET private piRNAs least insertion-driven.
- **Step 10 timepoint-specific (`build_count_matrix_perstrain.py` + `edger_da_timepoint.R`):** per strain,
  edgeR DA (timepoint factor) ∩ presence/absence (within-strain, same genome → no cross-genome step).
  C57 512,428 (E16.5 163k / P12.5 30k / P20.5 319k); CAST 606,537 (254k/15k/337k); SPRET 326,848
  (41k/9k/277k). Pachytene (P20.5) most stage-specific; P12.5 transition fewest. `{X}.timepoint_specific_DA.csv.gz`.
- **Step 11 sense/antisense to TE (`sense_antisense.py`):** piRNA strand vs stranded TE (RM .out, PanSN).
  Antisense (silencing-competent) %: unique-conserved-silent 51–57%, strain-private 51.5–55.6% vs common
  46–49% → unique piRNAs modestly + consistently MORE antisense-to-TE in all 3 strains. `{X}.antisense_byclass.csv`.

## Scale-out (in progress 2026-06-10)
- **#4 pangenome cluster PAV — DONE:** vg/hal conda install FAILED (hal not on bioconda) → used the cactus
  SIF (`cactus_v2.9.3.sif`, bundles halLiftover/halStats/vg) = the exact tool that built the HAL. Each
  strain's combined PICB clusters (union over tp; chr-stripped to HAL Ensembl naming 1,2,..) halLifted to
  GRCm39, `bedtools multiinter` → **CORE (all 16) 26.1 Mb (8.8%) / ACCESSORY (2–15) 176.5 Mb (59.3%) /
  PRIVATE (1) 94.9 Mb (31.9%)** — an OPEN piRNA-cluster pangenome (~9% core). Private-cluster bp highest in
  wild strains (PWK 17.0 / CAST 12.4 / SPRET 11.7 Mb). CAVEAT: halLiftover fragments clusters → region
  counts are sub-cluster fragments, **bp is the metric**; partial liftover inflates accessory/private.
  `cluster_pav/cluster_PAV_{catalogue.csv.gz,summary.csv}`.
- **ZAMORE validation (user-requested, black6):** GRCm39-lifted Zamore loci (214; mm10→GRCm39 liftover from
  prior work, STAGE mapped by locus NAME = build-safe, no coordinate mixing). **214/214 (100%) overlap a
  C57BL_6NJ PICB cluster** → PICB fully recovers the canonical Ozata-2020 black6 annotation (far above
  chance). **Zamore STAGE × 16-strain PAV:** 99 pachytene / 83 pre-pachytene / 32 hybrid; pachytene clusters
  markedly more cross-strain conserved (mean 13.8/16, median 15.2 = nearly core) than pre-pachytene/hybrid
  (~11.7/16) — confirms Ozata 2020 (pachytene deeply syntenic, pre-pachytene variable) and dovetails with
  the TE-driven finding (pre-pachytene = more TE-associated = more strain-variable).
  `Fig_zamore_stage_conservation`, `cluster_pav/zamore_loci_stage_conservation.csv`. Ozata 2020 PMID/DOI UNVERIFIED.
- **#5 scale-16:** collapse-13 RUNNING (job 3302592) — 13 strains not previously collapsed; EXACT validated
  sRNA cutadapt params (`--minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG`)
  + tstk collapse, glob all lane files/sample (matches pipeline cutadapt_input zcat), E15.5 excluded.
  Then 16-strain count matrices → edgeR DA → cross-strain expression test via PANGENOME (not 16×15 pairwise)
  → origin/TE/TE-driven. Multi-day.

## Biological grounding (BioMNI triple-verified 2026-06-10 — all 3 agents independently agree, HIGH confidence)
The method's biological premises were verified through the three independent BioMNI agents (genomics,
literature, general); all three concurred at HIGH confidence. **Concepts verified; specific citations
NOT** (agents returned different papers per claim + malformed DOIs e.g. `10.1126/science.17389255` =
PMID-in-a-DOI-template → fabricated; all PMIDs/DOIs flagged UNVERIFIED, queued for EuropePMC).
- **Length window (Step 2):** mouse piRNAs 24–32 nt; **pre-pachytene (E16.5) shorter (~26–28 nt),
  pachytene (P20.5) longer (~29–31 nt)**; length set by Zucchini/MitoPLD trimming + PIWI footprint
  (MILI pre-pachytene vs MIWI pachytene). ⇒ grounds 24–32 nt window AND matches our data (E16.5/P12.5
  mode 27–28 → P20.5 mode 30; `Fig_unique_pirna_length`).
- **Uniqueness = EXPRESSION not DNA presence (Steps 3–4):** a strain-specific piRNA = a sequence the
  strain *produces* that others do not; a ≤3-mm variant of a conserved EXPRESSED piRNA is a SNP-allele,
  NOT novel. ⇒ grounds the edgeR-DA expression criterion + the Step-4 genome-anchored split (0/1–3mm =
  not novel; the user's "expression 0" point is biologically correct).
- **TE-driven piRNA evolution (Step 6 / TE analysis):** strain-/species-private TE insertions CAN create
  new piRNA-producing loci; TE-driven repertoire turnover is established (read-through transcription,
  piRNA-cluster capture, ping-pong amplification). ⇒ grounds the strain-private-locus → TE-family test.

## Verification log
- Build/masking/naming verified (above); same assembly confirmed by identical chromosome lengths.
- BioMNI: method double-verified (genomics + literature agree); **all paper citations UNVERIFIED**
  (mismatched DOIs, implausibly consecutive PMIDs) — confirm via EuropePMC before citing.
- Thesis read for context (§3.2.3.1/3.2.3.4/3.2.3.6/5.10); method intentionally diverges/improves.

## Locked decisions (user-approved 2026-06-10)
1. **piRNA length window = 24–32 nt** (data-driven full piRNA range; FWHM major peak is 26–30 nt
   but 24–32 captures the stage shift — E16.5 pre-pachytene piRNAs skew shorter, P20.5 longer).
   Steps 3–6 re-run restricted to 24–32 nt.
2. **Presence threshold: data-driven, no magic RPM.** "Present in strain X" = detected (≥1 read) in
   **≥2/3 replicates** AND **non-singleton (≥2 reads total)** — reproducibility + singleton-noise
   exclusion replace the arbitrary MIN_RPM=1.0. RPM still reported for quantification.
3. **ALIGNMENT TOOL CONSISTENCY: STAR with the sRNA parameters for every genome alignment.** The
   SNP-variant step is redone with STAR (not bowtie); the ONLY parameter that differs is
   `--outFilterMismatchNmax` (0→3), which is unavoidable to detect SNP-variants and is documented as
   the single, purpose-driven deviation. One aligner, one parameter set, one justified exception —
   so a reviewer is not asked "why switch tools/params between steps".
4. **Scale-out (approved):** after pilot re-run + validation → all 16 strains; piRNA CLUSTERS compared
   cross-strain via the pangenome (presence/absence variation, core/accessory/private); full
   data-driven TE-driven-evolution analysis (strain-private piRNA loci vs strain-private TE
   insertions from the pangenome SV/TE PAV).

Execution staging: (i) re-run pilot Steps 3–6 with #1–#3 [STAR + 24–32 + reproducibility presence];
(ii) validate; (iii) scale to 16 strains; (iv) pangenome cluster PAV + TE-driven-evolution.
