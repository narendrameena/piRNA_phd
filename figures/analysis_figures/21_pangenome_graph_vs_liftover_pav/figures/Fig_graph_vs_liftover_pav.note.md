# Fig_graph_vs_liftover_pav

**Pangenome graph vs reference liftover: strain-variable piRNA clusters are mostly EPIGENETICALLY SILENCED, not genetically lost**

- **Shows:** A confusion (HAL-liftover cluster-PAV class × pangenome-graph sequence-PAV class); B strain×locus event breakdown (concordant / silencing / genetic loss) across PAV-ratio thresholds 0.5–0.95; C fraction of each liftover class that is sequence-shared in the graph; D the genetic-vs-regulatory model.
- **Result:** **99 %** of liftover-"private" piRNA clusters are sequence-**SHARED** in the graph (graph-core); **52–57 %** of strain×locus events are **SILENCING** (sequence present, no cluster) vs only **3–8 %** genetic loss; median strains with the **sequence** (graph) = **16** vs with a **cluster** (liftover) = **4**. So cross-strain piRNA-cluster variation is overwhelmingly **epigenetic** — the graph SEPARATES silencing from loss, which a single linear reference (liftover) conflates. Conclusion robust across all thresholds.
- **Why trustworthy:** 42,384 master loci = union of halLiftover-projected piCB clusters (GRCm39 frame). Graph PAV = `odgi pav -S` on the minigraph-cactus graph (paths imported via `vg convert -W` → P-lines → `odgi build`; 6,442 paths, 259 M nodes); per-strain PAV ratio = fraction of a locus's graph nodes the strain traverses. Threshold-robust (private-rescue = 99 % at every threshold 0.5–0.95).
- **How:** `code/02_odgi_build_pav.slurm` + `09_rebuild_plines_pav.slurm` (graph PAV), `03_compare_graph_vs_liftover.py`, `Fig_graph_vs_liftover_pav.py`.
- **Data:** `data/graph_pav_matrix.tsv`, `liftover_pav_matrix.tsv`, `graph_vs_liftover_comparison.tsv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
