# Methods — strain- and timepoint-specific UNIQUE piRNAs across 16 inbred mouse strains

Reproducible protocol (every tool, version, parameter, threshold and command). Citations carry PMID + DOI;
each is **[✓ EuropePMC]** (citation existence confirmed via the EuropePMC REST API) and each biological claim
carries its BioMNI verification status (genomics G / literature L / general Gn). Tool versions are pinned in
`workflow/envs/*.yaml`; alignment/trimming parameters are in `config/config_para.json`. This file is the
clean protocol; the running decision log with rationale is `unique_pirna/METHODS_review_notes.md`.

> Verification convention: **[✓ VERIFIED]** = independently confirmed by ≥2 BioMNI agents (the literature
> agent has been degenerate; double G+Gn confirmation is the project standard where it fails). Citations are
> database-confirmed (EuropePMC). Anything not yet confirmed is tagged **[⏳ PROVISIONAL]**.

---

## 0. Samples, genomes, exclusions

- **16 inbred strains**, sRNA-seq (single-end) at three spermatogenic timepoints: **E16.5 (16.5dpc)**,
  **P12.5 (12.5dpp)**, **P20.5 (20.5dpp)**, up to 3 replicates each (144 sample-replicates total).
- Pilot = C57BL_6NJ, CAST_EiJ, SPRET_EiJ; full run = all 16.
- **E15.5 mPGC samples EXCLUDED** (Black6-hybrid, not pure inbred).
- **Genomes:** each strain's own **unmasked** assembly (REL-2205); masked genomes suppress TE-derived
  piRNAs, so `dna` (not `dna_sm`/`dna_rm`) is used. Three chromosome-naming conventions coexist and are
  reconciled explicitly: PICB refFasta / RepeatMasker BED / cluster BEDs = `chrN`; GFF3 + STAR index =
  PanSN `{strain}#1#chrN`; pangenome HAL = Ensembl `1..19,X,Y,MT`. Identical assembly confirmed by
  identical chromosome lengths (chr1 = 191,890,962 across sources).
- Cross-strain comparison is done only through the **pangenome / liftOver**, never by mixing coordinates.

## 1. piRNA sequence layer — trimming + collapse

- **Tool/params (`config/config_para.json → params.srna`):** cutadapt
  `--minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG`
  (project-validated 3′ adapter), then `tstk collapse` to unique sequences with read counts encoded in the
  header `>seqid-nreads`.
- **piRNA identity = the unique sequence.** Outputs `results/collapse/{strain}-{tp}.{rep}.raw.fasta.gz`.
- Versions: cutadapt and tstk pinned in `workflow/envs/cutadapt.yaml`, `workflow/envs/tstk.yaml`.

## 2. QC and the piRNA length window (data-driven)

- Per-sample read totals + unique-sequence counts; count-weighted length distribution. 76–162 M reads/sample,
  7–19 M unique sequences; pooled length mode 27 nt.
- **Window = 24–32 nt**, validated on the strain-specific candidate set itself (captures 96.7% pooled;
  excludes the ~3% degradation/non-piRNA tail). The unique piRNAs are tighter than the bulk pool.
- **Biology [✓ VERIFIED, G+Gn, 2026-06-11]:** mouse piRNAs are ~24–32 nt with a PIWI-partner length split —
  **MILI-bound ~24–28 nt, MIWI-bound (pachytene) ~29–31 nt** — set by Zucchini/MitoPLD trimming and the PIWI
  footprint; this grounds both the window and the observed per-timepoint shift (E16.5/P12.5 mode 27–28 →
  P20.5 mode 30). Refs: Aravin et al. 2006 *Nature* 442:203 **[✓ EuropePMC]** PMID **16751777** DOI
  **10.1038/nature04916**; Gainetdinov et al. 2018 *Mol Cell* 71:775 **[✓ EuropePMC]** PMID **30193099** DOI
  **10.1016/j.molcel.2018.08.007**.

## 3. Strain-specific caller — presence/absence ∩ edgeR differential abundance

Improves on the thesis MGFR (arbitrary 0.5 marker score, replicate-blind, not FDR-controlled).

- **3a presence/absence:** present(X) = detected (≥1 read) in **≥2/3 replicates** AND non-singleton
  (≥2 reads total). Reproducibility alone is necessary but not sufficient (it admits a 29–37 M low-count
  noise tail), so it is combined with a data-driven abundance floor + significance.
- **3b data-driven floor + DA (`build_count_matrix*.py` → `edger_da.R`/`edger16.R`; edgeR):**
  per-timepoint count matrix (24–32 nt unique sequences × samples), `filterByExpr(group=strain)` sets the
  low-count floor from library sizes + smallest group size (no magic RPM); `glmQLFit` + `glmQLFTest`
  contrast **strain X vs mean(others)**, BH-FDR < 0.05 & logFC > 0 = significantly X-enriched; intersected
  with presence/absence (present ≥2/3 reps in X, absent <2/3 in every other strain) ⇒ strain-**specific**.
  `lib.size` = 24–32 nt window total per sample.
- **Why edgeR:** replicate-aware, FDR-controlled, library-size-normalised — none of which MGFR provides.

## 4. Genuinely-unique vs expressed-elsewhere — genome-anchored expression test (STAR, tool-consistent)

- **Principle [✓ VERIFIED, G+Gn]:** uniqueness is defined by **EXPRESSION, not DNA presence** — a sequence a
  strain *produces* that others do not. A ≤3-mm variant of a conserved expressed piRNA is a SNP allele, not a
  novel piRNA; a locus identical in another strain but transcriptionally silent there is still unique to X.
- **Method (`run_step4_map.sh` → `classify_step4.py`):** map each candidate of strain X to every other
  strain's unmasked genome with the **sRNA STAR parameters**, the only change being `--outFilterMismatchNmax`
  0→3 (the data-driven Poisson cutoff, §5). The reference sequence at the mapped locus (MD tag) is the piRNA
  strain Y would produce there; **Y expresses it iff that sequence ∈ Y's expressed pool** (Y's 24–32 nt
  collapse unique sequences). NM = SNP distance. Tool consistency: STAR 2.7.11b everywhere (supersedes an
  earlier bowtie sketch); indices verified unmasked.
- **4-way classification (`{X}.step4_classified.csv.gz`):** expressed-elsewhere exact (0 mm) / SNP-variant
  (1–3 mm) = NOT unique; unique:conserved-but-silent / unique:strain-private-locus = genuinely unique.
- **sRNA STAR parameters (`config_para.json → params.srna.STAR`):**
  `--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10
  --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000`.
- **Mapping-parameter rationale [✓ VERIFIED w/ CORRECTION, G+Gn, 2026-06-11]:** 0 mismatches + permissive
  multimapping (≤800) is correct for piRNAs, but the rationale is **(i) exact genomic origin** — piRNAs are
  genome-templated, so a mismatch is a SNP/sequencing error — **and (ii) piRNAs derive from repetitive
  TE/cluster sequence** so they legitimately map to many loci. It is *not* "exact target complementarity"
  (that concerns silencing, not mapping).

## 5. Data-driven mismatch cutoff (Poisson; no magic number)

- True SNP-variants of one piRNA follow Poisson(λ) over mismatches, λ = SNP-divergence × piRNA length;
  estimate λ = n(1mm)/n(0mm). Cutoff k* = smallest k with Poisson CDF(k;λ) ≥ 0.999.
- **Result:** λ = 0.219 (≈0.81%/bp at piRNA loci); CDF(≤3) = 0.9999 ⇒ **k\* = 3**.

## 6. Genomic origin, TE families, sense/antisense

- Map genuinely-novel piRNAs to the strain's own unmasked `chrN` genome; assign origin by category-split GFF3
  (genes/intergenic/lncRNA; CDS/5′UTR/3′UTR; exon/intron) with TE-precedence priority, reporting the full
  overlap profile. Multimapper RPM via featureCounts `-M --fraction -O` (fractional 1/N).
- **TE families of strain-private loci:** bedtools-intersect own-genome loci with per-strain RepeatMasker
  BED (`col4 = name|class/family`). LTR/ERVK(IAP) dominant, then LINE/L1, ERVL-MaLR, SINE. TE-derived
  fraction is a **lower bound** (STAR index = main chromosomes + MT only; ~20–26% of private piRNAs on
  unplaced contigs are unmappable/untyped).
- **Sense/antisense to TE** (stranded RepeatMasker .out, PanSN): antisense-to-TE = silencing-competent.
  Unique piRNAs are modestly but consistently more antisense-to-TE (51–57%) than common (46–49%).

## 7. PCA of piRNA expression (thesis Fig 5.21 method)

- Per timepoint, **DESeq2 1.42.1** size-factor normalisation → `prcomp` on log2(normalised+1); `filterByExpr`
  floor as in §3. Two feature sets: all-expressed (top-500 variable) and genuinely-unique. Strains form
  tight separated clusters; SPRET most distinct, increasingly with development.

## 8. Pangenome TE-driven-evolution test (decisive analysis; coordinate-based)

- Resource: existing 17-strain **minigraph-cactus** pangenome (`results/pangenome/output/`,
  GRCm39-referenced VCF + `.full.hal` + `.gbz`).
- Strain-private insertions: per multi-allelic site, ALT with len(ALT)−len(REF) ≥ 40 bp, private to X.
- **Confound found + fixed:** sequence-containment of a piRNA in a private insertion is confounded (active-TE
  family piRNAs are conserved across strains → match by family similarity without being produced there).
  **Fix (coordinate-based):** minimap2 the private insertion sequences back to X's own genome (asm5,
  mapq≥20) → insertion loci in X; intersect with each candidate's **production** locus. Tests production,
  not similarity.
- **Result/conclusion [✓ VERIFIED, G+Gn]:** new strain-private TE insertions overwhelmingly **propagate the
  active family's CONSERVED piRNAs** (new L1/ERVK copies); they do **not** preferentially create novel
  strain-private piRNA sequences. A real minority (uniquely-mapping: CAST 873, SPRET 801; ERVK/L1) are
  genuine TE-insertion-derived. Strain-private piRNA sequences are predominantly **divergence-driven at
  existing loci**. The per-insertion new-locus frequency is not quantified in prior literature → this
  quantification is novel.
- **TE-driven biology grounding [✓ VERIFIED, G+Gn]:** strain/species-private TE insertions *can* create new
  piRNA-producing loci (read-through transcription, cluster capture, ping-pong amplification). Refs:
  Aravin et al. 2008 *Mol Cell* 31:785 **[✓ EuropePMC]** PMID **18922463** DOI **10.1016/j.molcel.2008.09.003**;
  Brennecke et al. 2007 *Cell* 128:1089 **[✓ EuropePMC]** PMID **17346786** DOI **10.1016/j.cell.2007.01.043**.

## 9–11. Refinements, timepoint-specific, sense/antisense

- **9 uniquely-mapping refinement:** restrict TE-driven calls to strain-private piRNAs mapping to exactly one
  own-genome locus (removes TE-multimapping ambiguity).
- **10 timepoint-specific:** per strain, edgeR DA across timepoints ∩ within-strain presence/absence (same
  genome, no cross-genome step). Pachytene (P20.5) most stage-specific.
- **11 sense/antisense:** as §6.

## Phasing (developmental dynamics across spermatogenesis)

- **Method [✓ EuropePMC; method source]:** Almeida et al. 2025 *Genome Biol* 26:14 PMID **39844208** DOI
  **10.1186/s13059-025-03475-z** — load 24–32 nt alignments into R as GenomicRanges; for same-strand piRNAs
  use `GenomicRanges::follow()` to find the next 3′→5′ adjacent piRNA; phasing = enrichment of directly
  adjacent (+1 nt) pairs (`workflow/scripts/R/phasing_analysis.R`, per-chromosome streaming).
- **Alignment for phasing/peterplot/PCA/1U:** STAR 2.7.11b with **1 random coordinate per read**
  (`--outSAMmultNmax 1 --outMultimapperOrder Random --runRNGseed 777`); cutadapt 5.0. (The PICB BAM, in
  contrast, is used for PICB only.)
- **Mechanism [✓ VERIFIED, G]:** Zucchini (PLD6) cleaves the precursor; the loaded piRNA's 3′ end defines
  the next cleavage ~26–28 nt downstream → a phased same-strand series; distinct from ping-pong (same-strand
  ~26–28 nt 5′–5′ spacing vs cross-strand 10-nt 5′ overlap). Refs: Han et al. 2015 *Science* 348:817
  **[✓ EuropePMC]** PMID **25977554** DOI **10.1126/science.aaa1264**; Mohn et al. 2015 *Science* 348:812
  **[✓ EuropePMC]** PMID **25977553** DOI **10.1126/science.aaa1039**.
- **Developmental timing [✓ VERIFIED, G+Gn, 2026-06-11]:** phasing strongest at pachytene (~P20.5,
  MIWI-dominant), intermediate P12.5, weakest E16.5 (ping-pong-dominant). Our measured +1 phasing fraction:
  **E16.5 18.5% → P12.5 30.5% → P20.5 58%** (monotonic, peaks at pachytene). PIWI staging
  **[✓ VERIFIED, G+Gn]:** MIWI2/PIWIL4 nuclear ~E12.5–E16.5 (de novo TE methylation via DNMT3L/3A/3C);
  MILI/PIWIL2 perinuclear E12.5→adult; MIWI/PIWIL1 from pachytene (~P14). Refs: Carmell et al. 2007
  *Dev Cell* 12:503 **[✓ EuropePMC]** PMID **17395546** DOI **10.1016/j.devcel.2007.03.001**; Watanabe et al.
  2011 *Science* 332:848 **[✓ EuropePMC]** PMID **21566194** DOI **10.1126/science.1203919**.

## PICB clusters; proTRAC

- **PICB** (piRNA Cluster Builder): Ahrend, Konstantinidou, …, Haase 2025 *STAR Protoc* 6:103759
  **[✓ EuropePMC]** PMID **40220304** DOI **10.1016/j.xpro.2025.103759**. Run per strain×timepoint
  (replicates combined); completeness verified by clusters spanning chr1–19 + X.
  *Concurrency note:* the chunked driver must isolate its per-task temp directory (a shared temp +
  EXIT-trap caused silent xlsx corruption — fixed by per-task isolation).
- **proTRAC:** Rosenkranz & Zischler 2012 *BMC Bioinformatics* 13:5 **[✓ EuropePMC]** PMID **22233380** DOI
  **10.1186/1471-2105-13-5**.

## Scale-out to 16 strains (`unique16/`, `edger16/`)

1. **Collapse** all 16 strains (same cutadapt + tstk collapse as §1) → `results/collapse/` (144 files).
2. **Per-timepoint count matrix** (`build_count_matrix16.py`): 24–32 nt unique sequences × up to 48 samples;
   `min.total.count ≥ 15` early filter for tractability → `edger16/{tp}.counts.tsv.gz/.seqs/.samples`.
3. **edgeR DA** (`edger16.R`): the §3 logic, strain-specific among all 16 (present ≥2/3 reps in X, absent in
   all 15 others) → `edger16/{tp}.strain_specific_DA.csv.gz`.

## Pangenome cross-strain expression test = the 16-strain Step 4 (`PANGENOME_PHASE_README.md`)

Avoids 16×15 pairwise mapping by separating the two questions.

1. **Expressed pools** (`build_pools16.py`): per strain, 24–32 nt sequences with ≥2 reads (non-singleton),
   pooled over its 9 collapse files → `unique16/pools/{strain}.pool.txt.gz`.
2. **Exact cross-strain expression** (`classify_unique16.py`): load the 16 pools as a per-sequence 16-bit
   mask; a candidate is `expressed elsewhere (exact)` iff its sequence is in another strain's pool (0 mm,
   build-independent). Remainder = novel-sequence candidates emitted per strain/timepoint as FASTA.
3. **Production loci** (`cand_loci16.sh`): STAR-map novel candidates to the own genome with the §4 sRNA
   parameters (0 mm, 800 multimappers); BED of loci, chroms stripped PanSN→Ensembl. *Validated by smoke
   test (50 SPRET seqs → 96 loci; candidate ids + Ensembl chroms intact).*
4. **Cross-strain locus homology** (`lift_cand16.sh`, `lift_presence16.sh`): halLiftover each strain's loci
   X→GRCm39 through the cactus pangenome HAL, then lift the union GRCm39→each strain Y. A candidate's
   orthologous locus is present in Y iff it lifts to Y.
5. **Final classification** (`classify_unique16_locus.py`): novel-sequence + locus present in no other strain
   → `unique: strain-private locus`; + locus present elsewhere (not exact-expressed) → `unique:
   conserved-but-silent`. *(SNP-variant ≤3-mm refinement via per-strain allele fetch is documented as an
   add-on.)* → `unique16/final_classified.csv.gz`.

## Figures

Publication style (Liberation Sans; vector PDF/SVG + 300-dpi PNG; colourblind-safe; per-figure source data
in `SourceData_*.csv`; canonical thesis-Fig-4.4 strain order from `strain_order.py`). Narrative/mechanistic
figures (`unique_pirna/pangenome_te/`) are data-rich at genomic-coordinate + nucleotide resolution:
`Fig_biogenesis` (real ping-pong pair + phased run at chr2 coords), `Fig_locus_full_{IAP,L1}` (coverage +
TE + 1U + nucleotide pair with coordinate ruler + zoom callout), `Fig_locus_example_{IAP,L1}` (real coverage
across the 3 timepoints + verified cross-strain absence + antisense-piRNA inset), `Fig_concept_four_routes`
(real Step-4 sequence per route), `Fig_snp_variant_{nucleotide,coord}`, `Fig_pca_unique`, etc.

## Tool / environment versions
Pinned in `workflow/envs/*.yaml`. Key: STAR 2.7.11b; cutadapt 5.0; DESeq2 1.42.1; edgeR (envs/…);
samtools/bedtools (envs/ccTE, envs/samtools); minimap2 (envs/minimap); Trinity (envs/trinity);
Progressive Cactus / halLiftover / vg via `cactus_v2.9.3.sif`; PICB (R); tstk (collapse, peterplot).

## Verification status (this document)
All 11 cited references are **EuropePMC-confirmed** (PMID + DOI above). All biological claims used here are
**BioMNI-verified** (genomics + general; the literature agent has been degenerate, so double G+Gn confirmation
is used per the project standard) — see `analysis/claude_biomni_analysis/VERIFICATION_QUEUE.md` for the
per-claim log (M1–M7 VERIFIED; M6 with the mapping-rationale correction recorded above). The PICB-vs-SV
cluster-disruption conclusions (D1–D7) are outside this unique-piRNA methodology and remain partially
verified in the queue.
