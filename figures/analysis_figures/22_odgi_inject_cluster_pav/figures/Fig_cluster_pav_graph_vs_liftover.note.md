# Fig_cluster_pav_graph_vs_liftover

**Graph-inject vs HAL-liftover cluster-PAV: two independent methods cross-validate the piRNA-cluster conservation classification (99% agreement)**

- **Shows:** A confusion of HAL-liftover cluster class × graph-inject cluster class; B per-strain agreement (graph vs liftover); C graph class counts vs the coverage threshold (with the liftover-core reference); D the interpretation.
- **Result:** the graph-inject cluster-PAV (clusters placed on each strain's own graph path via `odgi inject`, then `odgi pav` node-coverage) and the HAL-liftover cluster-PAV **agree 99%** per strain×locus (median 4 strains each; classes 5608/23101/13481 vs 5983/22866/13535). Agreement is **99–100% for every strain, incl. divergent wild strains**. The graph is slightly more conservative (**graph-only = 0** — never adds spurious clusters; lift-only ~200–300/strain from the 0.4% fragment-boundary drops + coverage threshold). Robust across thresholds 0.05–0.2. → the cluster conservation classification is **method-robust** (cross-validated by a reference-free placement vs a reference projection). Contrasts theme 21 (graph SEQUENCE vs liftover CLUSTER = silencing); here both measure cluster presence → they concur.
- **Why trustworthy:** 364,924 native clusters injected onto the graph (0.4% dropped); cluster-presence = `odgi pav` node-coverage of each strain's injected cluster-paths at the same 42,384 master loci as the liftover PAV; two thresholds + per-strain breakdown reported.
- **How:** `code/01_build_inject_bed.py` → `02_inject_overlap.slurm` (inject) → `04_pav_cluster.slurm` (`odgi pav -p`) → `05_compare_cluster_pav.py` → `Fig_cluster_pav_graph_vs_liftover.py`. (`odgi overlap`/`untangle` for reference-free loci hit scale limits — see PIPELINE.)
- **Data:** `data/graph_cluster_pav_matrix.tsv`, `cluster_pav_comparison.tsv`; reuses theme-21 `liftover_pav_matrix.tsv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
