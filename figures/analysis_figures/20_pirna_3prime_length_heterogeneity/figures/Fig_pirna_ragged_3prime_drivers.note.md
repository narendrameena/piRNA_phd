# Fig_pirna_ragged_3prime_drivers

**Ragged 3′ ends are a BIOLOGICAL SIGNATURE of piRNA processing — abundant, strain-universal, secondary/ping-pong-linked, and pachytene-cluster-derived**

- **Shows:** A read- vs sequence-weighted raggedness (+ abundance); B per-strain raggedness across all 16 strains; C E16.5 ping-pong vs non-ping-pong raggedness; D TE-derived vs genic raggedness across development; E pachytene-cluster vs non-cluster (and high- vs low-expression cluster) raggedness; F the BioMNI-verified biology model.
- **Result:** ragged isoforms are **not noise** — (A) they **dominate the read pool** (59→86 % read-weighted vs 39→68 % by sequence) and are **more abundant** (median 189 vs 108 reads, P20.5); (B) **universal** across 16 strains (rise in **8/8**; wild ≈ classical); (C) **secondary (ping-pong) piRNAs are far raggeder, 66 % vs 39 %** (10-nt overlap z=7.2, p=3e-22; MILI↔MIWI2); (D) origin **flips** — TE-derived raggeder early (43 vs 37 %), genic raggeder in pachytene (70 vs 60 %); (E) **90 %** of pachytene piRNAs are **cluster-derived** (cluster 67 % vs non-cluster 57 %, p=9e-14; high-expression clusters raggedest, 73 %).
- **Why trustworthy:** exact-prefix raggedness (see Fig 1 validation). Each panel statistically tested (Fisher exact). Ping-pong = 10-nt 5′-5′ overlap z-score from cand_self16 coords; cluster overlap on self-genome coords vs 20.5dpp `clusters_fpm.bed`; TE/orientation from per-cand annotation. Biology **BioMNI 3/3-verified** (MIWI binds 30 nt · PNLDC1 trims · HENMT1 caps · A-MYB drives · Pnldc1-KO → longer untrimmed).
- **How:** `code/Fig_pirna_ragged_3prime_drivers.py`; computes `code/compute_{expression_ragged,perstrain_ragged,te_orientation,pingpong_e16,cluster_origin}.py`.
- **Data:** `data/SourceData_{expression_ragged,perstrain_ragged,pingpong_e16,cluster_origin}.csv` (+ TE/orientation recomputed from clean_2read × sense_antisense).

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
