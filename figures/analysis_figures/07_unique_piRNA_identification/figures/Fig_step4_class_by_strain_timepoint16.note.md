# Fig_step4_class_by_strain_timepoint16

- **Shows:** the 5-class unique-piRNA classification **annotated by strain AND timepoint together** — 3 panels (E16.5 / P12.5 / P20.5), each = 16 strains (canonical order; wild = red labels + shaded block) × stacked class composition (expressed-elsewhere · SNP-variant · low-quality · unique:conserved-but-silent · unique:strain-private). Above each bar: **% genuinely unique** and the **total candidate n** for that strain×timepoint.
- **Code:** `code/Fig_step4_class_by_strain_timepoint16.py` · **Data:** `data/SourceData_step4_class_by_strain_timepoint16.csv`
- **Pipeline:** see [`PIPELINE.md`](../PIPELINE.md). Source = `unique16/final_classified_clean.csv.gz` (klass5).

- **Relation to siblings:** [[Fig_step4_classification16]] Panel B is the class composition **pooled over timepoints**; `Fig_step4_class_strain_timepoint16` stacks **timepoint** and facets by **class** (the other marginal). This figure is the missing view — stack **class**, facet by **timepoint**.
- **Reading it:** the strain-private (purple) and conserved-but-silent (blue) genuinely-unique fractions are largest at **E16.5** (prepachytene) and in **wild-derived** strains; by P20.5 the candidates are dominated by SNP-variant/low-quality (not genuinely unique). Same data + caveats as [[Fig_step4_classification16]].
