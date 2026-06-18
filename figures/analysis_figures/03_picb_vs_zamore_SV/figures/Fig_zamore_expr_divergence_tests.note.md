# Fig_zamore_expr_divergence_tests
- **Shows:** A = expression concentration (Lorenz, Gini 0.90; ~25 loci → 90% of output); B = conserved master set (top-25 loci fire in every strain, ρ̄=0.88); C = divergence ↔ expression variability, split wild/classical (classical ρ=0.40 p=1e-9; wild ρ=0.15); D = wild slightly less concentrated (Gini, Mann-Whitney p=0.042).
- **Method:** quantitative locus FPM = Σ overlapping pangenome PICB clusters; Spearman/Mann-Whitney/Kruskal-Wallis. TE-class→divergence effect tested and REJECTED (KW p=0.38; not claimed).
- **Code:** `code/zamore_expression_tests.py`, `code/Fig_zamore_expr_divergence_tests.py` · **Data:** `data/zamore_locus_expression_P20.5.csv`, `data/zamore_expression_tests.txt`
- **Pipeline:** [`PIPELINE.md`](../PIPELINE.md)
