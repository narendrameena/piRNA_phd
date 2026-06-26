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

## Scope / limitation (honest)
Compared at the **lifted-cluster master loci** (theme-21's 42,384). A reference-FREE locus definition (to additionally
place **non-reference** clusters — the graph's unique potential) needs graph-space co-location (`odgi overlap`/`untangle`),
which hit scale/robustness limits at 365k injected paths on the 74 GB graph. `odgi pav` (index-light, the theme-21
workhorse) is the tractable route. The reference-free locus definition is the open piece (per-chromosome batching would
make `overlap`/`untangle` tractable if pursued).

## Outputs
`graph_inj.og` (74 GB, gitignored), `graph_cluster_pav_matrix.tsv`, `cluster_pav_comparison.tsv`,
`Fig_cluster_pav_graph_vs_liftover.{pdf,svg,png}`. Reuses theme-21 `graph.og`, `liftover_pav_matrix.tsv`,
`master_loci_pav_fixed.bed`. odgi/vg in `cactus_v2.9.3.sif`.
