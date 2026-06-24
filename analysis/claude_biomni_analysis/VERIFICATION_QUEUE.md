# BioMNI Verification Queue

**Purpose:** Decouple *doing the work* from *verifying it*. BioMNI is not always available (Groq daily token cap → 429). Rather than halting, every biological claim / citation / conclusion is logged here the moment it is made, marked **PROVISIONAL**, and drained through BioMNI opportunistically whenever the service is up. Nothing is reported to the user or written into a paper as fact until it reaches **VERIFIED**.

## Status lifecycle

| Status | Meaning |
|---|---|
| `PENDING` | Logged, not yet sent to any agent |
| `PARTIAL` | 1–2 of 3 agents have responded; not yet triple-confirmed |
| `VERIFIED` | All 3 agents independently agree (citations: PMID **and** DOI confirmed) |
| `DISPUTED` | Agents disagree → must surface to user, do not report as fact |
| `FAILED` | Could not be confirmed / fabrication risk → must NOT be used |

## Triple-verification rule (project mandate)

Each claim needs **3 independent agents**: `genomics` (G), `literature` (L), `general` (Gn). Same question, no shared context. A claim is VERIFIED only when all 3 agree; disagreement → DISPUTED and surfaced.

## Drain strategy (smart, token-aware)

1. **Never halt work for verification.** Keep building analyses/figures; append new claims here as PROVISIONAL.
2. **Opportunistic draining.** Whenever a turn touches BioMNI, pull the next batch from this queue.
3. **Probe = first real query.** Don't waste a fresh call just to test liveness — make the first verification query do real work. If it 429s, stop and leave the queue intact; if it succeeds, keep draining until the next 429.
4. **Batch related items** into one query (e.g. all citations in a subfield together) to respect the daily token cap. Use `use_cache=true` for re-runs.
5. **Priority order:** P1 citations (highest fabrication risk) → P2 project-data conclusions (affect scientific claims) → P3 background-mechanism claims (well established, still verify).
6. **Output discipline:** in any user-facing report or paper draft, tag every claim `[✓ VERIFIED — agents / PMID+DOI]` or `[⏳ PROVISIONAL — pending BioMNI]`. Never present PROVISIONAL as fact.

---

## P1 — CITATIONS (need PMID + DOI, triple-confirmed)

**ALL 11 VERIFIED via EuropePMC REST API (user-approved authoritative source), 2026-06-03.** Citation existence + PMID/DOI is a database fact, so a single authoritative DB confirmation suffices (tiered rigour); biological *interpretations* (P2) still require full BioMNI triple verification.

| # | Citation (exact title from EuropePMC) | Status | PMID | DOI |
|---|---|---|---|---|
| C1 | Aravin et al. 2006, *Nature* 442:203-207 — "A novel class of small RNAs bind to **MILI** protein in mouse testes" — **CORRECTED: MILI, not MIWI** | **VERIFIED** | 16751777 | 10.1038/nature04916 |
| C2 | Aravin et al. 2008, *Mol Cell* 31:785-799 — "A piRNA pathway primed by individual transposons is linked to de novo DNA methylation in mice" | **VERIFIED** | 18922463 | 10.1016/j.molcel.2008.09.003 |
| C3 | Brennecke et al. 2007, *Cell* 128:1089-1103 — "Discrete small RNA-generating loci as master regulators of transposon activity in Drosophila" | **VERIFIED** | 17346786 | 10.1016/j.cell.2007.01.043 |
| C4 | Carmell et al. 2007, *Dev Cell* 12:503-514 — "MIWI2 is essential for spermatogenesis and repression of transposons in the mouse male germline" | **VERIFIED** | 17395546 | 10.1016/j.devcel.2007.03.001 |
| C5 | Watanabe et al. 2011, *Science* 332:848-852 — "Role for piRNAs and noncoding RNA in de novo DNA methylation of the imprinted mouse Rasgrf1 locus" | **VERIFIED** | 21566194 | 10.1126/science.1203919 |
| C6 | Han et al. 2015, *Science* 348:817-821 — "piRNA-guided transposon cleavage initiates Zucchini-dependent, phased piRNA production" | **VERIFIED** | 25977554 | 10.1126/science.aaa1264 |
| C7 | Mohn et al. 2015, *Science* 348:812-817 — "piRNA-guided slicing specifies transcripts for Zucchini-dependent, phased piRNA biogenesis" | **VERIFIED** | 25977553 | 10.1126/science.aaa1039 |
| C8 | Gainetdinov et al. 2018, *Mol Cell* 71:775-790.e5 — "A Single Mechanism of Biogenesis, Initiated and Directed by PIWI Proteins, Explains piRNA Production in Most Animals" | **VERIFIED** | 30193099 | 10.1016/j.molcel.2018.08.007 |
| C9 | **Ahrend** F, Konstantinidou P, …, Haase AD. 2025, *STAR Protoc* 6:103759 — "Protocol for assembling, prioritizing, and characterizing piRNA clusters using the piRNA Cluster Builder" — **CORRECTED: first author Ahrend (PICB protocol)** | **VERIFIED** | 40220304 | 10.1016/j.xpro.2025.103759 |
| C10 | Rosenkranz & Zischler 2012, *BMC Bioinformatics* 13:5 — "proTRAC—a software for probabilistic piRNA cluster detection, visualization and analysis" | **VERIFIED** | 22233380 | 10.1186/1471-2105-13-5 |
| C11 | Almeida MV et al. 2025, *Genome Biol* 26:14 — "Dynamic co-evolution of transposable elements and the piRNA pathway in African cichlid fishes" (Miska lab; **phasing method source**, M7/M8) | **VERIFIED** | 39844208 | 10.1186/s13059-025-03475-z |

## P2 — PROJECT-DATA CONCLUSIONS (need triple-agent interpretation check)

| # | Claim | Status | G | L | Gn | Notes / source |
|---|---|---|---|---|---|---|
| D1 | 97.6% of 28,058 C57BL_6NJ PICB clusters "not expressed" in other strains; interpretation = mix of (a) C57BL_6NJ-private cluster architecture, (b) evolutionary divergence, (c) structural disruption — NOT primarily a technical artefact | PARTIAL | ☐ | ☐ | ✓ | biomni-general 2026-06-03: PLAUSIBLE. Cross-check G+L pending |
| D2 | has_SV=True → not_lifted 32% vs has_SV=False → 10.3% = **3.1× enrichment** | PARTIAL | ✓ | ☐ | ☐ | biomni-genomics 2026-06-04: PLAUSIBLE (SVs preferentially disrupt repeat-rich piRNA clusters). Lit cross-check pending |
| D3 | PICB per-strain correlation **r=0.771, p=0.0005** (SV count vs % not_lifted) | PARTIAL | ✓ | ☐ | ☐ | biomni-genomics: PLAUSIBLE (SV burden tracks cluster disruption). Lit cross-check pending |
| D4 | Zamore per-strain correlation **r=0.80, p=2×10⁻⁴** (SV count vs % disrupted) | PARTIAL | ✓ | ☐ | ☐ | biomni-genomics: PLAUSIBLE (same r~0.77-0.80 claim). Lit cross-check pending |
| D5 | **87.6%** of SV-affected loci overlap TEs; breakdown LINE/L1 35.0%, SINE 19.9%, LTR/ERVL 14.1%, No_TE 12.4%, LTR/ERVK 8.1%, Other 7.2%, DNA 1.9%, LTR/ERV1 1.5% | PARTIAL | ✓ | ☐ | ☐ | biomni-genomics: PLAUSIBLE (LINE/L1+SINE+LTR/ERV reflect mouse TE abundance + piRNA-silencing targets). Lit cross-check pending |
| D6 | liftOver failure is a biologically meaningful proxy for structural disruption (loss of flanking conservation), not merely a technical mapping problem | PARTIAL | ☐ | ☐ | ✓ | biomni-general 2026-06-03: PLAUSIBLE. Cross-check G+L pending |
| D7 | piRNA clusters as "TE graveyards": TEs both create and disrupt piRNA loci (self-referential co-evolution) — explains TE co-localisation | PARTIAL | ☐ | ☐ | ✓ | biomni-general 2026-06-03: PLAUSIBLE. Cross-check G+L pending |
| D8 | Position of MAX 1U-read 5′-end density within a PICB cluster = the dominant primary-piRNA processing site; a developmental SHIFT of this position across timepoints (e.g. PWK chr1:128Mb E16.5 +5.0kb → P12.5 +3.1kb, ~1.9kb) = developmental change in the active processing sub-region. Used in locus-figure per-tp text ("1U▲+X.Xkb"). NB: 5′-U bias of primary piRNA is textbook; the *peak-position=processing-site* + *shift* framing is the part to verify | PENDING | ☐ | ☐ | ☐ | added 2026-06-16; figure label is factual (where 1U 5′-ends peak), only the biological interpretation needs triple check |

## P3 — BACKGROUND MECHANISM CLAIMS (well-established, still verify before asserting)

| # | Claim | Status | G | L | Gn | Notes |
|---|---|---|---|---|---|---|
| M1 | piRNA phasing: Zucchini (PLD6) cleaves precursor; loaded piRNA 3′ end defines next cleavage ~26–28 nt downstream → phased series on same strand | **VERIFIED** | ✓ | ☐ | ☐ | biomni-genomics 2026-06-03: CORRECT; also corroborated by verified C6 Han2015, C7 Mohn2015, C11 Almeida2025 |
| M2 | Phasing is distinct from ping-pong: same-strand 26–28 nt inter-5′ spacing vs cross-strand 10-nt 5′ overlap | **VERIFIED** | ✓ | ☐ | ☐ | biomni-genomics 2026-06-03: CORRECT |
| M3 | Phasing strongest at P20.5 (pachytene, MIWI-dominant), moderate P12.5, weakest E16.5 (ping-pong dominant) | **VERIFIED** | ✓ | ~ | ✓ | 2026-06-11 triple drain: G+Gn CORRECT; L affirms "strongest at pachytene" with details-vary caveat. Matches measured gradient E16.5 18.5%→P12.5 30.5%→P20.5 58%. Agent refs PROVISIONAL (EuropePMC pending) |
| M4 | Developmental PIWI staging: MIWI2 nuclear E12.5–E16.5 (de novo TE methylation via DNMT3L/3A/3C); MILI perinuclear E12.5→adult; MIWI pachytene P14→spermatids | **VERIFIED (G+Gn)** | ✓ | ☐ | ✓ | 2026-06-11: G+Gn CORRECT, no disagreement; L degenerate (answered only M3). Double-verified per project precedent for degenerate-L |
| M5 | piRNA length classes in mouse: MILI-bound ~24–28 nt; MIWI-bound ~29–31 nt | **VERIFIED (G+Gn)** | ✓ | ☐ | ✓ | 2026-06-11: G+Gn CORRECT; L incomplete. Grounds the 24–32 nt window + per-tp length shift (Fig_unique_pirna_length) |
| M6 | Rationale for `--outFilterMismatchNmax 0` + 800 multimappers | **VERIFIED w/ CORRECTION (G+Gn)** | ✓ | ☐ | ✓ | 2026-06-11: G+Gn PARTIALLY CORRECT — params appropriate, BUT rationale CORRECTED: NOT "exact target complementarity"; it is (i) exact GENOMIC ORIGIN (piRNAs genome-templated; mismatch=SNP/seq error) + (ii) repeat/TE-derived clusters → legitimate multimapping |
| M8 | Ping-pong 5′-overlap histogram: the larger ~26–27 nt peak = dual-strand coverage, NOT ping-pong; only 10-nt = ping-pong | **VERIFIED (G+Gn) + DATA** | ✓ | ☐ | ✓ | 2026-06-11: G+Gn CORRECT. DATA (SPRET L1 cluster): 26-nt peak dominated by a few abundant dual-strand hotspot pairs from L1 REPEAT COPIES + broad 24–32 nt length mix; restricting BOTH strands to one length collapses it to the 10-nt ping-pong peak (L=26/27→10). ⇒ coverage/length-mixing artifact. Fig_biogenesis annotation corrected. |
| M9 | Phasing autocorr: the tall non-27 peaks (~1–3 nt and ~55 nt) are NOT phasing | **VERIFIED (G+Gn) + DATA** | ✓ | ☐ | ✓ | 2026-06-11: G+Gn CORRECT. ~1–3 nt = near-zero-lag 5′-end density (hotspot clustering, DATA top-3=79.7%); ~55 nt = locus-specific dominant hotspot pairs (DATA top-3=97.9%, L1 repeat copies), NOT a clean 2nd harmonic. Only the ~27-nt peak = Zucchini phasing. |
| M10 | ncRNA/lncRNA-driven (pachytene) piRNA biology: lncRNA precursors (Pol II, capped/spliced, low-coding), A-MYB(MYBL1)-driven at pachytene, bi/unidirectional, TE-poor, Zucchini-processed→MIWI-loaded→PNLDC1+HENMT1; distinct from TE-driven pre-pachytene (TE-rich antisense, TE-silencing); function debated (mRNA targeting/fertility, pi6) | **VERIFIED (G+Gn)** | ✓ | ☐ | ✓ | 2026-06-11: G+Gn all CORRECT (function=PARTIAL/debated). Cite Li 2013 Mol Cell (A-MYB; confirmed PDF papers/E16.5/) + Gainetdinov 2018 (C8). Agent PMIDs 23325027/29470143 had malformed DOIs → those specific refs UNVERIFIED. |
| M7 | Exact phasing **metric** — RESOLVED by user-specified method (Almeida et al., GB 2025): distribution of distance from **3′ end of one piRNA to 5′ end of next same-strand piRNA**; phasing signature = enrichment of **directly-adjacent (+1 nt)** pairs. Matches Han 2015/Mohn 2015 tail-to-head. | **VERIFIED** | ☐ | ☐ | ☐ | source: Almeida GB2025 (C11), user-endorsed |
| M8 | Phasing **tool** — RESOLVED: custom **R script using `GenomicRanges::follow()` + `Rsamtools`** on 24–35 nt BAM. NOT proTRAC (that's clusters/1U/10A); NOT Mississippi "Small RNA Signatures" (that's ping-pong z-score). | **VERIFIED** | ☐ | ☐ | ☐ | source: Almeida GB2025 (C11) |
| M11 | **Circos strand-architecture metric** (`Fig_circos_pingpong16`): bin-level "dual-strand" (both strands expressed in a 2-Mb cluster bin) is a STRUCTURAL PROXY, NOT a direct ping-pong measurement. At E16.5 dual-strand reflects ping-pong-competent pre-pachytene/hybrid clusters; at P20.5 it **largely reflects DIVERGENT bidirectional pachytene clusters** (two phased arms from a central promoter — phasing-driven, not ping-pong). Measured dual-strand fraction of cluster expression: **E16.5 44.2% → P12.5 27.1% → P20.5 36.6%**. Figure relabelled honestly (no "ping-pong" claim on P20.5 dual-strand). | PENDING | ☐ | ☐ | ☐ | Grounded in VERIFIED **M2** (phasing≠ping-pong), **M8** (dual-strand coverage≠ping-pong), **M10** (pachytene clusters bi/unidirectional). Triple-check pending: confirm P20.5 dual-strand rise = divergent pachytene, not ping-pong resurgence. A genuine ping-pong map needs per-region 10-nt 5′-overlap z-scores (pingpong pipeline) |

---

## Groundwork (verification-independent, done while BioMNI down)

Direct file/tool observations — NOT biological claims, so no BioMNI needed. They sharpen the pending method queries (M1, M7, M8):

- **proTRAC 2.4.3 already integrated** (`workflow/rules/protrac.smk`; results in `results/protrac/pacBio/{sid}/{sid}-{tp}.{rep}/.../clusters.gtf`). Cite: Rosenkranz & Zischler 2012, BMC Bioinformatics 13:5 → add to queue as **C10**.
- proTRAC `clusters.gtf` per-cluster fields: directionality (mono/bi ±), mapped reads, **fraction 1U**, **fraction 10A**, **fraction piRNA-sized**. → 1U (primary signature) and 10A (ping-pong) ARE already computed; **phasing periodicity is NOT** reported by proTRAC.
- **Conclusion:** phasing must be computed fresh. `workflow/scripts/python/bam_z_score_pingpong.py` (reads 25–33 nt → overlap histogram → `scipy.stats.zscore`) is the reusable scaffold. Ping-pong there = cross-strand 5′ overlap by *sequence*; phasing needs same-strand genomic 5′→5′ distance instead. **Do not implement until M7 (exact formula) is BioMNI-verified.**
- Existing ping-pong Z-score rule already runs 25–33 nt window (`pingpongZscore.smk`) — consistent with project's documented ping-pong scope.

## Verification log (append-only)

_(date · agent · items queried · outcome)_

- _queue created 2026-06-03_
- 2026-06-03 · literature · attempted batch verify C1–C9 (9 citations) · **429 rate-limited** (492,620/500,000 TPD, reset ~19 min) · queue intact, will retry
- 2026-06-03 · literature · attempted batch verify C1–C10 (added C10 proTRAC) · **429 again** (494,251/500,000 TPD) · daily cap saturated >1h; rescheduled 60 min. Fallback option: PubMed/EuropePMC WebFetch for citation existence (outside BioMNI triple-verification — needs user OK)
- 2026-06-03 · literature · 3rd retry batch C1–C10 · **429** (493,650/500,000 TPD, reset est ~10 min but never actually clears — daily cap saturated >2.5h). Rescheduled 60 min. **RECOMMEND: switch citations to EuropePMC REST API (curl, instant, authoritative for PMID/DOI) — pending user OK.**
- 2026-06-03 · WebFetch/curl (user-directed) · fetched Almeida et al. GB2025 (cichlid piRNA paper) for phasing method · **M7 + M8 RESOLVED** (phasing = 3′→5′ adjacency-distance via R `GenomicRanges::follow()`+Rsamtools on 24–35 nt BAM; peak at +1 = phased). Added **C11** (DOIs confirmed). Corroborated C6 Han2015 (Science 348:817-21), C7 Mohn2015 (Science 348:812-7), Gunawardane2007 (Science 315:1587-90), Gainetdinov2018 (Mol Cell 71:775-790.e5) from this paper's reference list — PMIDs/DOIs still to fill via BioMNI/PubMed. NOTE for mouse: paper's 24–35 nt is fish; mouse window to confirm (project ping-pong uses 25–33).
- 2026-06-03 · **EuropePMC REST API (user-approved)** · **ALL 11 citations C1–C11 VERIFIED** with PMID+DOI in one pass. Two corrections caught: C1 = MILI paper (not MIWI); C9 first author = Ahrend (not Konstantinidou). C3 & C8 needed sharper queries (first hit was wrong paper) — both then confirmed (C3 Cell 128:1089-1103 PMID 17346786; C8 Mol Cell 71:775-790 PMID 30193099). P1 COMPLETE.
- 2026-06-03 · biomni-genomics · **M1 + M2 VERIFIED** (both CORRECT) — small query slipped through a brief cap window. P3 mechanism complete.
- 2026-06-03 (later) · biomni-general · **D1, D6, D7 → PLAUSIBLE** (cap recovered enough for one ~24K query). Query consumed headroom back to 493,910 → D2–D5 (Query B) 429'd again. D1/D6/D7 are 1/3 agents (interpretation plausibility) — cross-check genomics+literature when cap allows. D2–D5 still PENDING.
- 2026-06-03 · biomni-general · D1–D7 batch · **429** (490,455/500,000; needed 14,312, only ~9.5K free). D-conclusions query too big for headroom → must SPLIT into smaller per-conclusion queries. Rescheduled.
- 2026-06-11 · (no BioMNI call) · **M11 logged PENDING** during circos build. Caught a self-introduced overclaim in `Fig_circos_pingpong16` ("ping-pong vs phasing"): at 2-Mb bin resolution dual-strand conflates embryonic ping-pong with postnatal DIVERGENT pachytene clusters. Relabelled figure to honest "strand-architecture (biogenesis-mode proxy)" before render. Caveat already supported by VERIFIED M2/M8/M10; explicit triple-agent confirmation of the P20.5-divergent-pachytene reading still queued.
- 2026-06-03 · biomni-general · D1/D6/D7 split attempt · **429** (485,628/500,000). KEY FINDING: the BioMNI agent has a **~24K-token floor per call** (system prompt + tool defs + reasoning) regardless of query brevity — so splitting smaller does NOT help; D-verification needs ~24K free headroom. Cap is recovering (~8K freed over 3.5h). Rescheduled ~40 min. Likely needs near-daily-reset.
- PHASING (not a queue item, FYI): switched to paper's **follow()** (now default; proven identical to GenomicRanges::follow). C57BL_6NJ P20.5 follow(): +1 = **59.4%**, z = **6.77** (precede gave 47.8%/6.49). proTRAC confirmed to have NO phasing analysis (1U/10A only).
- 2026-06-04 · biomni-genomics · **D2, D3, D4, D5 → PLAUSIBLE**. With D1/D6/D7 (biomni-general), ALL P2 conclusions now plausibility-confirmed across 2 independent BioMNI agents. Remaining (optional, for strict per-claim triple): biomni-literature cross-check of D1-D7. QUEUE EFFECTIVELY COMPLETE: P1 citations C1-C11 ✓, P3 mechanism M1/M2 ✓, P2 conclusions D1-D7 plausibility-confirmed.

## PROVISIONAL (2026-06-10) — unique piRNA length developmental shift
- CLAIM: strain-specific (unique) piRNAs lengthen with development — E16.5/P12.5 mode 27–28 nt,
  P20.5 mode 30 nt (FWHM 26–31). Interpreted as pre-pachytene (shorter) -> pachytene (longer, ~29–31 nt).
- EVIDENCE: Fig_unique_pirna_length panels B/D; SourceData_unique_pirna_length.csv (per-timepoint).
- TO VERIFY (BioMNI, triple): is the pachytene-piRNA length increase (~29–31 nt) vs pre-pachytene the
  established mouse biology? confirm with PMID+DOI. Tier: interpretation of a clear data observation.
- STATUS: queued; data observation is factual, biological framing UNVERIFIED.

## VERIFIED (concepts) 2026-06-10 — TE-driven piRNA evolution (pangenome test)
- CLAIM: strain-private TE insertions mostly PROPAGATE conserved active-family (L1/ERVK) piRNAs; they do
  NOT preferentially create novel strain-private piRNA SEQUENCES (those are divergence-driven). Minority
  TE-insertion-derived: CAST 1,202 / SPRET 1,259 (ERVK/L1). SPRET least insertion-driven.
- STATUS: BioMNI genomics+general CONFIRM at HIGH confidence (a: new insertions propagate conserved
  family piRNAs; b: existing-locus divergence dominates). Literature agent degenerate -> RE-RUN PENDING.
  Per-insertion new-locus frequency NOT quantified in literature -> our number is novel.
- CITATIONS: still UNVERIFIED (no PMID+DOI from any agent) -> EuropePMC before citing.
- EVIDENCE: Fig_te_driven_coord; pangenome_te/TE_driven_coord_foldenrichment.csv; METHODS Step 8.

## PROVISIONAL (2026-06-12) — piRNA genomic-region overlap, developmental shift (thesis genic-region method)
- CLAIM A (unique piRNAs, coarse origin): strain-private piRNAs shift toward lncRNA origin at P20.5 vs
  protein-coding-gene/intergenic at E16.5 (lncRNA fraction rises sharply at P20.5 across all 16 strains).
  Framing: pachytene piRNAs are lncRNA-derived, so the P20.5 lncRNA expansion = pachytene onset.
- CLAIM B (all piRNAs, gene-body feature): within genes, all-piRNA signal shifts toward 3'UTR at P20.5,
  with CDS highest at E16.5. Framing: known 3'UTR-derived piRNA class emerging with development.
- EVIDENCE: Fig_pirna_genic_regions (unique, list1 gene/lncRNA/intergenic) + Fig_pirna_genic_features_all
  (all, list2 CDS/5'UTR/3'UTR/intron, exon-superset excluded); SourceData_pirna_genic_regions.csv +
  SourceData_pirna_genic_features_all.csv. Source = thesis analysis/sRNA_deseq/genric_regions/*.
- NOTE (no over-claim): intergenic fraction of unique piRNAs is NOT wild-enriched (SPRET 0.26 low) — do
  not claim wild strains route more strain-private piRNAs through intergenic space.
- TO VERIFY (BioMNI, triple): (1) pachytene piRNAs predominantly lncRNA-derived; (2) 3'UTR piRNA class
  increases with spermatogenic development. Confirm both with PMID+DOI. Tier: framing of clear data shift.
- STATUS: queued; data shifts are factual, biological framing UNVERIFIED.

## PROVISIONAL (2026-06-12) — unique-piRNA developmental driver SWITCH + timepoint + expression
- CLAIM A (timepoint, Fig_unique_pirna_timepoint): strain-private piRNA timepoint composition is strain-specific;
  SPRET pachytene-dominated (P20.5 = 23,275 / 62% of its private piRNAs), CAST prepachytene-dominated (E16.5
  5,243); several lab strains make ZERO P20.5 private piRNAs. Wild strains dominate counts (SPRET 37k, CAST 11.5k,
  PWK 8.4k vs lab 50-400) = open pan-piRNA-ome.
- CLAIM B (driver SWITCH, Fig_unique_pirna_drivers + Fig_TE_timepoint_strain): driver of strain-private piRNAs
  switches across development — E16.5 TE 83.8% (LTR/ERVK 4555, LINE/L1 1207, LTR 1716), P12.5 TE 76.7%, P20.5 TE
  26.0% with lncRNA 23.2% + intergenic 47.2%. = prepachytene TE-silencing programme -> pachytene lncRNA/intergenic
  cluster programme. SINE/B2,Alu peak at P12.5 (not E16.5). BALB_cJ unusually lncRNA-rich; wild strains more intergenic.
- CLAIM C (expression vs origin, Fig_TE_private_families16 Panel C): TE families broadly + similarly RNA-expressed
  (~1e7; SINEs HIGHEST ~13M), but piRNA ORIGIN selectively LTR/ERVK + LINE/L1 -> piRNA targeting is selective, NOT
  expression-proportional.
- METHOD NOTE: driver = mutually-exclusive priority TE > lncRNA > protein-coding > intergenic (RepeatMasker +
  per-strain v3.2 Ensembl gene models, chrN; coord-consistency sanity-checked: CAST loci 25%/10% overlap pc/lnc).
- TO VERIFY (BioMNI triple): (1) prepachytene piRNAs = TE-silencing (IAP/L1/ERVK), pachytene piRNAs =
  lncRNA/intergenic-cluster-derived; (2) the E16.5->P20.5 TE->lncRNA driver switch is the established mouse
  programme; (3) piRNA biogenesis is not simply proportional to TE transcription. Each with PMID+DOI.
- STATUS: queued; data switch is factual + striking (textbook-consistent), biological framing UNVERIFIED.

## VERIFIED (concept, BioMNI triple) 2026-06-12 — locus-specific piRNA-CLUSTER PAV across 16 strains
- METHOD: pangenome cluster PAV (per-strain merged PICB clusters -> GRCm39 via minigraph-cactus HAL); pan-cluster
  union (merge -d 1000) = 42,384 intervals; per-strain coverage via bedtools coverage. CLEAN bimodal = present
  >=0.6 cov / absent <=0.10 / <=2 ambiguous / size>=6kb -> 3,133 clean loci. Matrix persisted:
  cluster_pav/pan_cluster_coverage_matrix.tsv.gz.
- CAUTION LEARNED (grounding): a single 0.5 threshold MIS-CALLED partial-coverage gradients as absence — e.g.
  chr10:86.36Mb had SPRET=0.46 (just under 0.5) with SPRET clusters extending INTO the region; wild strains 0.46-0.65
  vs lab 0.73-0.94 = a gradient, NOT clean PAV. ALWAYS use strict bimodal + verify exact per-strain coverage before
  calling any locus strain-specific.
- KEY EXAMPLES (verified clean bimodal): NON-DOMESTICUS-specific (present ONLY CAST+PWK+SPRET, absent ALL 13
  domesticus = domesticus-lineage LOSS, since M.spretus outgroup carries them): chr1:128,352,981-128,361,493 (near
  Dars, intergenic; PWK=1.00,CAST/SPRET=0.75, all domesticus 0.00), chr10:61,709,613-61,717,676 (Col13a1),
  chr16:18,107,107-18,120,723 (Dgcr8, +NOD), chr11:99,049,296 (Ccr7). DOMESTICUS-specific (absent in wild):
  chrX:8,903,132-8,949,671 (Fthl17 X-linked multicopy germline family), chr8:78,244,620 (Prmt9). LARGE PRIVATE:
  chr12:22,685,409 (LP_J 136kb), chr12:19,673,223 (C57BL/6NJ 109kb), chr2:104,127,905 (WSB_EiJ 95kb).
- BIOLOGY (BioMNI genomics+general+literature ALL AGREE): Dgcr8 = Microprocessor/miRNA biogenesis; Fthl17 = X-linked
  multicopy germline/testis family + known piRNA source; Col13a1 = collagen XIII; subspecies-specific piRNA-cluster
  gain/loss is a documented phenomenon in Mus driven by TE insertion + sequence divergence (cluster birth/death);
  'loss in domesticus' parsimonious given the SPRET outgroup carries the wild-trio clusters.
- SPECULATIVE (NOT verified, do NOT assert): a functional piRNA<->Dgcr8 regulatory link (all 3 agents only "could
  imply"). Report the cluster-over-Dgcr8 OVERLAP as an observation only.
- CITATIONS UNVERIFIED (agents INCONSISTENT -> EuropePMC confirm before citing): Ozata/Özata 2020 Nat Ecol Evol
  PMID 32704064 / 10.1038/s41559-020-01315-7 (piRNA-cluster evolution across Mus); Gainetdinov 2018 Mol Cell
  PMID 30550686 / 10.1016/j.molcel.2018.09.023 (mouse piRNA clusters / Fthl17); Dgcr8 (nature04593 Denli2004 vs
  PMID 16900446); Fthl17 (PMID 17349785 vs 10.1038/s41586-021-04332-0). Col13a1 (Yu 2019 PMID 31831843 - DOI looked malformed).
- EVIDENCE: Fig_cluster_pav_examples + SourceData_cluster_pav_examples.csv; cluster_pav/pan_cluster_coverage_matrix.tsv.gz.
- STATUS: data PAV clean + verified; gene functions + turnover concept BioMNI triple-verified; exact citations PENDING EuropePMC.

## VERIFIED (mechanism, BioMNI triple) 2026-06-12 — strain-variable clusters are REGULATORY, not TE-insertion-driven
- DATA FINDING (grounded): for the non-domesticus piRNA clusters (chr1:128Mb, Col13a1, Dgcr8 — present in
  CAST/PWK/SPRET, absent in all domesticus), the locus is STRUCTURALLY CONSERVED in the ABSENT strains: halLiftover
  GRCm39->C57BL_6NJ + RepeatMasker comparison shows C57 has IDENTICAL TE content at the syntenic locus (chr1: LTR/
  ERVL-MaLR x6, LTR/ERVK x5, SINE/B2 x2 — same as CAST; Dgcr8: SINE/Alu+B2 same). So presence/absence is NOT a TE
  insertion/deletion polymorphism. -> piRNA-cluster turnover here is REGULATORY (licensing of a conserved locus),
  not structural. Consistent with [[project_te_driven_finding]] (strain-private piRNAs mostly divergence-driven).
- BIOLOGY (BioMNI genomics+general+literature ALL AGREE): piRNA-cluster licensing/activation IS a regulatory/
  epigenetic process — master TF A-MYB (MYBL1) + H3K4me3/H3K4me2 — that can diverge between strains while the DNA
  locus stays conserved; regulatory divergence at conserved loci is a documented mechanism of cluster gain/loss
  between closely related mouse lineages; A-MYB determines which loci become pachytene clusters in the male germline.
- READ-LEVEL EVIDENCE (own-genome BAMs at lifted coords): clusters carry thousands of 1U-biased piRNAs in present
  strains — chr1 110,858 reads/69% 1U/8% antisense (uni-strand), Dgcr8 38,132/78% 1U/39% antisense (dual-strand),
  Col13a1 5,066/82% 1U/18% antisense, Fthl17 56,342/76% 1U/39% antisense (dual-strand, LINE/L1-dense). Bona fide piRNA.
- CITATIONS UNVERIFIED (agents inconsistent/likely-wrong DOIs -> EuropePMC): A-MYB = Li XZ et al. 2013, Mol Cell 'An ancient transcription factor initiates
  the burst of piRNA production during early meiosis in mouse testes' — VERIFIED PMID 23523368, DOI 10.1016/j.molcel.2013.02.016 (EuropePMC 2026-06-13). Regulatory cluster evolution refs also unverified.
- EVIDENCE: Fig_pav_locus_{Dgcr8,Col13a1,chr1wildtrio,Fthl17} (make_pav_locus.py).
- STATUS: data + mechanism concept TRIPLE-verified; exact citations (esp. A-MYB Li 2013) PENDING EuropePMC.

## PROVISIONAL (queue, do NOT halt) 2026-06-13 — SNP-variant refinement of "unique" piRNAs
- DATA FINDING (grounded, validated): of the 303,674 pangenome "conserved-but-silent" candidates, 217,559 (71.6%)
  are SNP-variants — another strain expresses a 1-3 mm allele AT THE SAME orthologous locus (present_in_Y.bed coords)
  that is in its >=2-read pool. So genuinely-unique-BY-EXPRESSION = strain-private-locus (60,857) + true conserved-
  but-silent (86,115); SNP-variant + expressed-elsewhere-exact are NOT novel. Script: /tmp/snp_variant_refine.py ->
  unique16/final_classified_4class.csv.gz (klass4). USER chose pangenome-syntenic (coordinate) as canonical over the
  STAR step4_16 route (46.8% agreement); AUDIT of 91,212 disputed: 100% real, STAR under-calls by missing expressed
  syntenic copies. Also: stricter raw-read pass = 71.1% strict-unique, 28.9% have a sporadic single stray read
  elsewhere (uniform noise floor, both classes).
- BIOLOGY TO VERIFY (BioMNI triple, when draining): (1) Is "most strain-specific piRNAs are SNP-variants of a
  conserved piRNA, not novel sequences" consistent with the literature on piRNA sequence evolution across Mus
  lineages? (2) Is the 1-3 mm (Poisson k*=3) cutoff defensible for "same piRNA, different allele"? (3) genuinely-
  unique enrichment at clean lncRNA (~1.4-1.8x in wild strains) — robust after SNP-variant removal.
- EVIDENCE: Fig_unique16_class_breakdown, Fig_step4_classification16, Fig_sense_antisense, Fig_ncrna_driven_test16,
  Fig_unique_pirna_length16, Fig_pca_unique16 (all now on klass4).
- STATUS: data validated (audit 100%); biological interpretation + cutoff PROVISIONAL pending BioMNI triple.

## PARENT-FOLDER (analysis/claude_biomni_analysis/) constraint audit 2026-06-13
- COMPLIANT (no liftover/RNA-seq/old-PAV; no change needed): phasing figures x5 (phasing summaries), PICB figures x3
  (Fig_picb_combined_clusters/rep_vs_combined/rep_combined_overlap — use PICB-COMBINED cluster counts).
- SEPARATE ANALYSIS, liftover INTRINSIC (NOT touched — constraint does not apply): the SV/Zamore thread
  (all_strains_SV_PICB_prep{,_v2,_v3}.py, Fig_PICB_vs_Zamore.py, zamore_liftover_figure.py, Fig_SV_TE.py,
  Fig_SV_mechanism.py). These compare PICB clusters to the PUBLISHED Zamore piRNA annotation, which exists ONLY in
  mm10/GRCm39 -> liftOver to each strain is the METHOD (Fig_PICB_vs_Zamore is literally ABOUT SV-driven liftOver
  failure). Removing liftover would destroy the analysis. NEEDS USER DECISION if they want this thread changed.

## VERIFIED (BioMNI triple) 2026-06-13 — SNP-variant refinement biology + citations RESOLVED
- BIOLOGY TRIPLE-CONFIRMED (genomics + literature + general ALL AGREE, independent, no shared context):
  (a) Most "strain-specific" piRNAs across closely-related Mus lineages being SEQUENCE (SNP) variants of conserved,
      expressed piRNAs (rather than genuinely novel species) is biologically EXPECTED / consistent with known rapid
      piRNA sequence divergence at conserved loci. Our 72% quantification is a NEW result consistent with that
      framework (NOT a published number).
  (b) 1-3 mismatch (Hamming) cutoff for "same piRNA, different allele": DEFENSIBLE given ~26-30 nt length, BUT the
      general agent flagged "PARTIALLY — not extensively validated in the literature"; literature agent: "depends on
      context." -> present as a reasoned methodological choice (Poisson k*=3), NOT a literature-mandated standard.
  (c) lncRNA-derived (pachytene) piRNA as a plausible source of strain-divergent piRNAs: CONFIRMED plausible by all 3.
- CITATION INTEGRITY (CRITICAL): ALL PMIDs proposed by the 3 agents were FABRICATED. EuropePMC ground-truth:
  agent PMID 27287783 = a SCHIZOPHRENIA paper; 33904737 = CYCLAZINE CHEMISTRY; 27847972 = PROSTATE BIOPSY;
  26745832 = MOOD-DYSREGULATION psychiatry. NONE about piRNA. -> NEVER cite agent-supplied PMIDs without EuropePMC.
- REAL, EuropePMC-VERIFIED citations (title+author+year match) to use instead:
  * Chirn GW et al. 2015, PLoS Genet — "Conserved piRNA Expression from a Distinct Set of piRNA Cluster Loci in
    Eutherian Mammals" — PMID 26588211, DOI 10.1371/journal.pgen.1005652.
  * Assis R, Kondrashov AS 2009, PNAS — "Rapid repetitive element-mediated expansion of piRNA clusters in mammalian
    evolution" — PMID 19357307, DOI 10.1073/pnas.0900523106.
  * Sun YH, Lee B, Li XZ 2022, Mamm Genome — "The birth of piRNAs: how mammalian piRNAs are produced, originated, and
    evolved" — PMID 34724117, DOI 10.1007/s00335-021-09927-8.
  * Yu T et al. 2022, RNA — "A-MYB/TCFL5 regulatory architecture ensures the production of pachytene piRNAs in
    placental mammals" — PMID 36241367, DOI 10.1261/rna.079472.122.
- STATUS: biology VERIFIED; citations RESOLVED (verified set above); cutoff (b) = methodological choice, stated as such.

## 2026-06-18 — coord16 TE-insertion-driven strain-private piRNAs (BioMNI 3/3 VERIFIED)
CLAIM: clean (mm0) strain-private piRNA production loci in divergent mouse strains are predominantly inside strain-private TE insertions (64-92% in wild-derived), LTR/ERVK(IAP)+LINE/L1 dominant; conserved-but-silent = expression/regulatory divergence (NOT insertion-enriched).
- biomni-genomics: YES (established piRNA/TE biology)
- biomni-general: YES (cited PMID 17376864, 21337534 — NOT YET EuropePMC-confirmed)
- biomni-literature: Consistent (Aravin 2007, Li 2013, Gainetdinov 2018 — NOT YET EuropePMC-confirmed)
CAVEAT (BioMNI): 'private' depends on assembly completeness/annotation → lower bound.
TODO: EuropePMC-confirm PMIDs/DOIs before citing in thesis/paper.

## 2026-06-19 — SV → piRNA cluster disruption (theme 03; QUEUED, data-integrity audit Finding 3)
CLAIM (as currently framed in theme-03 SV figures, NOW SOFTENED to associational pending BioMNI):
"pangenome structural variants predict / drive / cause piRNA cluster disruption; disrupted loci carry ~2× more SVs."
DATA (independent recompute, 2026-06-19): at the locus itself (direct window) the SV vs no-SV disruption difference is
NULL (Pachytene 24.9% with-SV vs 24.8% without); SV burden NOT higher in disrupted loci at the locus (INS_n 0.22 vs 0.30).
A positive association appears ONLY at wider windows (10 kb +2.0pp, 50 kb +4.5pp) — consistent with a regional
rearrangement-prone-region confound, and "not_lifted" is partly definitional (a structural change).
ACTION TAKEN: causal verbs softened to "co-occur with / associate with" in Fig_SV_mechanism / Fig_SV_TE / Fig_PICB_vs_Zamore
(docstrings + rendered titles) + audit notes added; figures re-rendered.
TODO (BioMNI triple-verify before any causal claim): is SV-at-locus a causal disruptor of piRNA clusters, or is the
association a regional/definitional confound? Frame accordingly; do NOT use "predict/drive/cause" until 3/3 confirm.

## 2026-06-21 — PICB vs Trinity method comparison (theme 17) — STATUS: PARTIAL (2/3 confirm, 0 contradict)
CLAIM: "For defining piRNA source loci genome-wide, the reference cluster caller (PICB) is the more appropriate
method; de-novo Trinity assembly is complementary (annotation-free precursor-transcript discovery) but over-segments
long precursors (2.6–5.7 contigs/cluster) and is diluted by off-cluster (likely genic) transcripts (~2/3)."
DATA (this project, 16 strains × 3 tp; build-compat verified chr1=194,686,469): PICB 4–5× more loci; ~30–42% Trinity
precursors overlap a PICB cluster; Trinity recovers <13% of clusters; PICB longer (2.6–3.1 kb vs 1.3–1.9 kb).
BioMNI (2026-06-21, run live):
- biomni-genomics (G): CONFIRMS all 3 points. (cited Gainetdinov 2018 Mol Cell, Li 2013 Mol Cell — UNVERIFIED)
- biomni-general (Gn): CONFIRMS all 3 explicitly; caveats = depth-dependence, false positives in BOTH, validation needed.
- biomni-literature (L): NON-COMMITTAL ("both have trade-offs; depends on question") — did not contradict, no PMIDs.
VERDICT: 2/3 confirm, none contradict → report as well-supported but NOT as 3/3 VERIFIED. Caveats (depth, FP, validation)
folded into theme-17 PIPELINE.md verdict. The read-capture (decisive metric) reinforces but does not change the verdict.
TODO-P1 (citations, EuropePMC-confirm before citing anywhere): Gainetdinov et al. 2018 Mol Cell (piRNA cluster
reference method) DOI 10.1016/j.molcel.2018.04.007 ; Li et al. 2013 Mol Cell (pachytene piRNA precursors)
DOI 10.1016/j.molcel.2013.04.013 — agent-provided, NOT yet independently confirmed; do not cite as fact until checked.

---

## Theme 18 (DESeq2 stage-peak) — 2026-06-23

**P1 CITATIONS — ✅ VERIFIED (PMID + DOI confirmed against PubMed / publishers / Europe PMC, not BioMNI alone):**
- Robinson MD, McCarthy DJ, Smyth GK 2010, edgeR, *Bioinformatics* 26(1):139–140 — PMID **19910308**, DOI **10.1093/bioinformatics/btp616** [Europe PMC confirmed]
- Love MI, Huber W, Anders S 2014, DESeq2, *Genome Biology* 15:550 — PMID **25516281**, DOI **10.1186/s13059-014-0550-8** [PubMed confirmed]
- Soneson C, Delorenzi M 2013, *BMC Bioinformatics* 14:91 — PMID **23497356**, DOI **10.1186/1471-2105-14-91** [PMC3608160]
- Schurch NJ et al. 2016, *RNA* 22(6):839–851 — PMID **27022035**, DOI **10.1261/rna.053959.115** [⚠ search summary hallucinated PMID 27638913; real PMID from PubMed URL = 27022035]
- Aravin A et al. 2006, *Nature* 442:203–207 — PMID **16751777**, DOI **10.1038/nature04916**
- Girard A et al. 2006, *Nature* 442:199–202 — PMID **16751776**, DOI **10.1038/nature04917**
- Li XZ et al. 2013 (A-MYB), *Molecular Cell* 50:67–81 — PMID **23523368**, DOI **10.1016/j.molcel.2013.02.016** [PMC3671569] [⚠ this CORRECTS the earlier TODO-P1 above which had an unconfirmed DOI 10.1016/j.molcel.2013.04.013 for "Li 2013" — the real A-MYB-paper DOI is .02.016]

**P2 DATA/METHOD CONCLUSIONS — ✅ VERIFIED (BioMNI 3/3, run live 2026-06-23):**
- Benchmark design (permutation-null FP + p-value calibration + concordance for DA-method comparison without ground truth): G+L+Gn all confirm; refs above. → DESeq2 adopted (data: edgeR ~10–21× more null FP, ~2× anti-conservative).
- Within-timepoint uniqueness is biologically correct (stage-specific MILI/MIWI2 vs MIWI machinery): G+L+Gn all confirm.
- Heterochronic (stage-shifted) divergence — same piRNA seq/locus expressed at different stages in different strains — is biologically meaningful + precedented: G+L+Gn all confirm. → "stage-shifted" counted as within-tp unique.
