# Theme 22 — graph-native cluster-PAV via `odgi inject` (vs HAL-liftover)

## Goal
Redo the piRNA-cluster PAV the reference-free way with `odgi inject`, and **keep both** alongside the existing
HAL-liftover cluster-PAV (theme 21 / `cluster_pav`). Liftover projects clusters onto GRCm39 (reference-biased); inject
places each strain's clusters on its **own graph path** (reference-free), so the two are independent methods for the
same question — *which strains have a piRNA cluster at this locus?* See `REVIEW_NOTES.md` for the design + tool checks.

## Method (code/)
1. `01_build_inject_bed.py` — map **364,924** native piCB clusters onto the fragmented strain paths
   (`STRAIN#0#chr#offset`; offset from path name, length from `vg paths -E`). **0.4% dropped** (1,471 fragment-boundary
   + 51 non-graph-chrom). Emits an inject BED in path-space + the cluster→strain path-groups.
2. `02_inject_overlap.slurm` — `odgi inject` (reuses theme-21 `graph.og`) → `graph_inj.og` (74 GB, 364,924 cluster-paths).
   *(The `odgi overlap` co-location step here was O(n)/path → ~300 h for 365k paths → killed.)*
3. `03_untangle_clusters.slurm` — `odgi untangle` to project clusters to GRCm39: **crashed** (an `sdsl` step-index
   assertion on the injected graph). Aborted.
4. `04_pav_cluster.slurm` — **`odgi pav -p`** (inject cluster-paths grouped by strain) at the theme-21 master loci →
   graph cluster-PAV = fraction of each locus's nodes covered by each strain's CLUSTER-paths. Fast (~10 min), no crash.
5. `05_compare_cluster_pav.py` — graph-inject cluster-PAV vs HAL-liftover cluster-PAV.

## Result
- The graph-inject and HAL-liftover cluster-PAVs **agree 99%** per strain×locus (median 4 strains each; classes
  graph 5608/23101/13481 vs liftover 5983/22866/13535 core/disp/private). Robust across coverage thresholds 0.05–0.2.
- The graph is slightly **more conservative**: **graph-only = 0** (never adds spurious clusters), lift-only ~200–300
  /strain (from the 0.4% fragment-boundary drops + the coverage threshold). Agreement is **99–100% for every strain**,
  including the divergent wild strains (SPRET/CAST/PWK).
- **Interpretation:** two *independent* methods (reference-free graph placement vs reference projection) **cross-validate**
  the cluster conservation classification → it is method-robust. (Contrast theme 21: graph SEQUENCE-PAV vs liftover
  CLUSTER-PAV measure *different* things → they disagree = silencing; here both measure cluster presence → they concur.)

## Non-reference clusters — caught + characterised + confounder-checked (steps 6–13)
The reference-free piece, resolved a simpler way than graph co-location: re-lift each strain's piCB clusters (tagged with
per-cluster IDs) to GRCm39 (`06_identify_nonreference.slurm`); a cluster that yields **no** lifted interval = non-reference
(sequence absent from GRCm39). **1,393 (0.4%)** caught; per strain 42–132, wild-derived toward the high end.
- `07_characterize_nonref.py` — **93%** overlap a TE (RepeatMasker), dominated by **LTR/ERVK + LINE-1** (young active
  elements) = the piRNA–TE arms-race signature.
- `08`/`12_multimapping.py` — expression. `clusters_fpm.bed` cols = [chrom,start,end,**allFPM,uniqFPM**,strand,tp]
  (per `build_arch_switch.py`). Non-ref clusters are **25% multimapping** (vs 0.2% reference — TE-driven); on all-primary
  FPM they look more expressed (12.4 vs 6.9), but on **UNIQUE reads they are comparable to reference (6.6 vs 6.4)** →
  genuinely expressed, NOT 'higher' (the all-primary signal was multimapping-inflated; corrected).
- `09`/`10_te_burden.py`/`11_confounding.py` — RIGOROUS cross-strain evolution test + confounders. Per-strain non-ref
  TE-insertion burden from the deconstructed VCF (`10`; SPRET 76k ≫ C57BL_6NJ 2.3k = clean positive control). **Non-ref
  cluster count tracks the DIRECT TE-insertion burden: Spearman ρ=0.61, p=0.012**; unique-read non-ref piRNA share ρ=0.51,
  p=0.044; wild>classical p=0.001 / 0.021. CONFOUNDERS (`11`): NOT total-output (partial r=0.50 ✓); multimapping CORRECTED
  (survives on unique reads ✓); lift-artifact NOT divergence-differential (`13_flank_lift.slurm`: clean-insertion fraction
  does not fall with divergence, ρ=+0.45 ✓). **CAVEAT** — the effect is wild-vs-classical driven (within-classical n.s.;
  the 4 wild strains are phylogenetically non-independent), so treat as a group contrast, not a continuous gradient.
- Biology BioMNI-triple-verified (`biomni_verify_nonref.py`; all 3 agents 'Established': Aravin 2007, Girard 2010,
  Gainetdinov 2015, Frazer 2011).
- `Fig_nonreference_clusters.{pdf,svg,png}` — the confounder-checked 4-panel figure.

**Co-location + TE age (steps 14–16):**
- `14_colocation.sh`/`15_colocation_aggregate.py` — reference-free co-location by all-vs-all strain→strain halLiftover
  (240 pairs; each non-ref cluster lifted into every other strain). **66% strain-private** (insertion in no other
  strain), **34% shared**; sharing **follows phylogeny (62% within the wild/classical group** = inherited by descent);
  a ~117 tail present in 14–15 strains but absent from GRCm39 (C57BL/6J-lineage deletions / reference gaps).
- `16_te_age.py` — **TE-age test** (RepeatMasker % divergence from consensus; lower = younger). Non-reference cluster
  TEs are **younger** than reference/conserved: median **14.3% vs 16.4%, p=7×10⁻⁶⁸**; holds WITHIN the driver families
  (L1 12.4/13.9, ERVL-MaLR 20.0/20.7, ERV1 11.2/12.8; all p<10⁻¹¹) = insertion age, not family composition (ERVK
  already youngest, n.s.). → conserved = OLD shared TEs; non-reference = YOUNG strain-private insertions.
- `Fig_nonreference_colocation_age.{pdf,svg,png}` — co-location + TE-age, 4-panel.

## Scope / limitation (honest)
Cluster-PAV compared at the **lifted-cluster master loci** (theme-21's 42,384). The non-reference clusters above are
caught by re-liftover (non-lifting = absent), characterised by TE/expression/evolution, and confounder-checked — but
their reference-FREE **co-location** (are the 1,393 shared between strains or strain-private?) still needs graph-space
`odgi overlap`/`untangle`, which hit scale/robustness limits at 365k injected paths on the 74 GB graph. Now tractable
(the subset is 1,393, not 365k) via a restricted-search `odgi overlap` if pursued — the remaining open piece.

## Outputs
`graph_inj.og` (74 GB, gitignored), `graph_cluster_pav_matrix.tsv`, `cluster_pav_comparison.tsv`,
`Fig_cluster_pav_graph_vs_liftover.{pdf,svg,png}`. Non-reference: `data/nonref/{strain}.nonref.bed` + `.summary`,
`nonref_te_summary.csv`, `nonref_te_families.csv`, `nonref_expression_summary.csv`, `te_insertion_burden.json`,
`multimapping_test.csv`, `evolution_test.csv`, `confounding.csv`, `flank_genuine_insertion.csv`, `colocation.csv`,
`te_age_test.csv`, `te_age_byfamily.csv`, `Fig_nonreference_clusters.{pdf,svg,png}`,
`Fig_nonreference_colocation_age.{pdf,svg,png}`. Reuses theme-21 `graph.og`, `liftover_pav_matrix.tsv`,
`master_loci_pav_fixed.bed`. odgi/vg/halLiftover in `cactus_v2.9.3.sif`.
