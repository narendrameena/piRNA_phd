# Fig_filter_order

**Effect of applying the 27/30 nt filter BEFORE vs AFTER DESeq2 → filter-before is more sensitive**

- **Shows:** A overall stage-peak candidates per stage, Order A (filter BEFORE DESeq2) vs Order B (DESeq2 on full 25–32 then subset); B per-strain (E16.5).
- **Result:** **filtering BEFORE recovers ~6–7 % more** stage-peak candidates (E16.5 3,777 vs 3,551; P12.5 9,717 vs 9,259; P20.5 2,473 vs 2,371), Jaccard 0.85–0.93, BEFORE ≥ AFTER across strains. The two orders differ only in size factors, dispersion trend and the **BH denominator** (27/30 subset vs full 25–32) — filtering first reduces the multiple-testing burden and focuses normalization/dispersion on the relevant length class.
- **Why trustworthy:** identical DESeq2 one-vs-rest + ≥2-read presence/absence; only the filter position changes. Order A = the production candidate set.
- **How:** `analysis/.../deseq_filter_order_test.R` (SLURM `run_filter_order.sh`) → `deseq16_lenfilt/{tp}.candidates_orderA_before.csv.gz` + `*.filterorder_perstrain.csv`; `code/Fig_filter_order.py`.
- **Data:** `data/SourceData_Fig_filter_order.csv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
