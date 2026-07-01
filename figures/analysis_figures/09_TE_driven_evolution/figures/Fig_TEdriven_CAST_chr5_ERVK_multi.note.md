# Fig_TEdriven_CAST_chr5_ERVK_multi — CAST/EiJ TE-driven piRNA cluster (ERVK provirus, chr5:21.46 Mb)

**Locus:** CAST/EiJ-private ERVK provirus (RLTR20B4 + RMER19A + MurERV4_19-int), chr5:21,457,300–21,460,300.
**Shows:** a single TE insertion present in only CAST (1/16 strains) that created a piRNA cluster — 14,312
 P20.5 primary piRNAs, ~98% 1U, ~100% ANTISENSE-to-ERVK (silencing). MULTI: Panel A = standard FPM bar chart (CAST expresses; 15 strains TE-absent ○); Panel B = per-timepoint coverage + antisense-silencing second bar; Panel C = base resolution.

**How generated** — shared FASTQ→cutadapt→STAR→PICB→fetch_primary pipeline
([`../11_locus_catalogue/PIPELINE.md`](../11_locus_catalogue/PIPELINE.md) §1–8) + the coordinate-based
TE-driven test (private insertion ∩ PICB cluster, real-TE + high-1U + high-antisense filter; see this
theme's `PIPELINE.md` and `code/`).
- **Code:** `code/make_te_driven_locus_multi.py` + `code/pav_clusters.py`
- **Command:** `make_te_driven_locus_multi.py CAST_EiJ 5 21457300 21460300 "ERVK provirus (RLTR20B4 + RMER19A + MurERV4_19-int)" 15 Fig_TEdriven_CAST_chr5_ERVK_multi`
- **Data:** `data/TE_driven_COORDINATE_CAST_EiJ.csv` (coordinate-verified TE-driven candidates)
- **[finding · ⏳ pending BioMNI]:** "TE insertion → antisense-1U silencing cluster" is a project conclusion (VERIFICATION_QUEUE), not yet triple-confirmed.  **Formats:** PDF+SVG+PNG.
