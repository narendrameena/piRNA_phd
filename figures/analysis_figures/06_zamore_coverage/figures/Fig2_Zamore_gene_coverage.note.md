# Fig2_Zamore_gene_coverage

**Zamore coverage: sRNA coverage of Zamore pachytene piRNA genes**

- **Shows:** sRNA coverage of Zamore pachytene piRNA genes
- **How:** generate_P12_5_P20_5_figures.py. Coverage = fraction of each Zamore locus outer span overlapped by PICB cluster intervals (bedtools-style interval intersect; mm10 (= GRCm38) annotation, NO liftover).
- **Data:** P12_5_P20_5_zamore_coverage_per_gene.csv + claude_biomni_figures inputs
- **Provenance:** C57BL/6 P12.5/P20.5 (external data).

Full raw→figure pipeline: [`PIPELINE.md`](../PIPELINE.md). Originals under `analysis/claude_biomni_analysis/`.
