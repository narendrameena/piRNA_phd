# Fig_picb_combined_clusters

**PICB clusters per strain (combined-replicate run)**

- **Shows:** Grouped bars: #PICB clusters per strain × timepoint, 16 strains, canonical order.
- **How:** Fig_picb_combined_clusters.py. PICB cluster calling (R PICB) on STAR sRNA BAMs (unmasked strain genome, piRNA params); final `clusters` sheet. Replicates pooled before PICB.
- **Data:** SourceData_PICB_cluster_counts.csv (replicate=='combined')
- **Provenance:** combined-replicate PICB run.

Full raw→figure pipeline: [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.

- **Error bars (replicates):** ±SD across the 3 sequencing replicates per strain × timepoint; the combined/mean bar value ≈ the replicate mean. See [PIPELINE.md](../PIPELINE.md).
