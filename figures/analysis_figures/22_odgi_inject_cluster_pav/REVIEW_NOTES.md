# Theme 22 — graph-native cluster-PAV via `odgi inject` (REVIEW / DESIGN NOTES)

LOCAL design notes for the reference-free cluster-PAV. Decisions, verified tool interfaces, and the plan, captured
before the heavy build. (Keep both methods: the existing HAL-liftover cluster-PAV **and** this odgi-inject one.)

## Goal & rationale
Redo the piRNA-cluster PAV the **reference-FREE** way with `odgi inject`, alongside the existing HAL-liftover
cluster-PAV (`cluster_pav/`). The two stay side by side:
- **HAL-liftover** projects each strain's clusters onto GRCm39 → cluster presence in one linear frame. Reference-biased:
  divergent strains fragment under the lift (SPRET 18 k native clusters → 145 k lifted intervals) and non-reference /
  TE-insertion clusters cannot lift at all → mis-called private.
- **odgi inject** places each strain's clusters directly on **its own graph path** → cluster loci defined in graph
  space, reference-free: captures non-reference clusters, no lift fragmentation, divergent strains handled natively.

This is the cluster (regulatory) analog of theme 21's graph sequence-PAV — done on the graph instead of via liftover.

## Verified tool interfaces (odgi v0.9.0 in cactus_v2.9.3.sif — checked, NOT assumed)
- `odgi inject -i graph.og -b targets.bed -o out.og` — the BED is **over PATH SPACE of the graph**: col1 = a graph
  **path name**, col2/col3 = range **in that path's own coordinates**, col4 = annotation name. Each record becomes a
  NEW path in the output graph. ⇒ I must supply **strain fragment paths + LOCAL coords**, not native chr coords.
- `odgi overlap -i out.og -R cluster_paths.txt [-s cluster_paths.txt]` — "find the paths touched by given input paths."
  Run on the injected cluster-paths (restrict the search space to cluster-paths with `-s`) ⇒ for each cluster-path,
  which OTHER strains' cluster-paths overlap it in graph space = **cross-strain cluster sharing**.
- `odgi flatten -b out.bed` — path↔node linearization (fallback to compute node-set overlap directly if `overlap` is
  too slow at scale).

## Key obstacle — fragmented strain paths (the reason native coords don't inject directly)
Graph strain paths are `STRAIN#0#<chr>#<offset>` (many contig fragments per chromosome; the 4th PanSN field = the
fragment's START position in that chromosome). GRCm39 is a single path/chrom (`GRCm39#0#<chr>`). Native cluster coords
(`cluster_pav/{strain}.clusters.bed`, bare chrom, chr coordinates) must be mapped to (fragment_path, local coords):
- fragment chr-span = `[offset, offset+length)`; **offset** from the path name, **length** from `vg paths -E` (gbz).
- a native cluster `(chr, cs, ce)` → the fragment with `offset ≤ cs` and `ce ≤ offset+length`; local = `(cs−offset, ce−offset)`.
- clusters spanning a fragment boundary / assembly gap → drop or split (must quantify the drop rate; flag it).

## Pipeline (planned)
1. `vg paths -E` (gbz) → per-fragment lengths → fragment map per strain.  **[RUNNING — /tmp/graph_path_lengths.txt]**
2. Build inject BED: native clusters → `(STRAIN#0#chr#offset, local_start, local_end, "CLU|strain|clusterid")`. Report
   the fraction of clusters successfully mapped (vs boundary/gap drops).
3. `odgi inject -i graph.og -b inject.bed -o graph_inj.og` — SLURM big-mem (loads the 33 GB og). ~380 k cluster-paths.
4. `odgi overlap -R cluster_paths.txt -s cluster_paths.txt` → co-located cluster-path pairs (record the strain of each).
5. Union-find the cluster-paths by overlap → **graph cluster loci** → distinct strains per locus → core / dispensable /
   private (the reference-free cluster-PAV) → `data/graph_inject_cluster_pav.tsv`.
6. **Compare, keep BOTH:** graph-inject cluster-PAV vs HAL-liftover cluster-PAV (`cluster_PAV_catalogue`). Expected:
   the graph **unifies** divergent-strain clusters that liftover shatters, and **recovers** non-reference clusters
   liftover drops → fewer spurious "private", more correctly-shared loci. Figure panel + note; commit.

## Risks / open questions
- **Scale:** ~380 k cluster intervals (≈24 k/strain × 16) on a 259 M-node graph — inject is heavy → SLURM big-mem;
  `odgi overlap` on 380 k paths may need per-chromosome batching.
- **Fragment-boundary clusters:** quantify drop rate at step 2 (expect small; flag if not).
- **Non-reference clusters:** confirm they inject (they should — they sit on the strain's own path, which is in the graph).
- **"Keep both" framing:** graph-inject = reference-free cluster loci; liftover = the projected-to-GRCm39 view. Neither
  is discarded — they answer the same question (which strains express a cluster here) in two coordinate frames, and the
  disagreement is itself the reference-bias readout (mirrors theme 21's sequence-PAV finding).
