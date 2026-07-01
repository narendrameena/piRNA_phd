# Theme 23 — non-reference piRNA-cluster locus atlas (per-locus figures, with the odgi graph)

## Goal
A visual atlas: one detailed per-locus figure (make_pav_locus style) for each piRNA cluster **genetically absent from
GRCm39** — the theme-22 non-reference findings (graph-confirmed absent, theme-22 steps 20/21) — with the pangenome
GRAPH (odgi inject + pav) supplying the cross-strain evidence that liftover cannot.

## The figure (code/make_nonref_locus.py)
Three panels (panels B/C reuse the strain-anchored `pav_clusters` drawing blocks; only panel A is graph-native):
- **(A) Pangenome × timepoint** — per-strain PICB-cluster FPM (log), each bar **split by HEIGHT: solid = + strand /
  pale = − strand** (piRNA strand architecture, from each strain's own sRNA via the co-location coords). Below it a
  **GRAPH-COVERAGE STRIP** (odgi pav): each strain's *continuous* coverage (full ≈ 1.0 · partial/segmental, **%
  labelled** · absent), GRCm39 = 0 — the clever odgi angle (liftover gives only binary present/absent). Tally:
  present / partial / genetically-absent / present-but-silent → the **genetic-loss vs regulatory-silencing** read-out
  (sequence carried but unexpressed = present-but-silent).
- **(B)** the carrier strain's per-timepoint sRNA coverage + TE + gene tracks.
- **(C)** base resolution (5′ arrow red = antisense-to-TE = silencing).

## Loci (code/gen_nonref_loci_batch.py)
All **19 mapped** shared-subset non-reference loci (those in theme-22 `graph_check_pav3.tsv`, so the graph-PAV panel A
works). **Featured** (curated names): `WSB_chr4_L1` (major L1, 213 FPM, present 10/16, partial 5), `FVB_chr1_ERVK`
(C57BL/6J-lineage ERVK, segmental ~50–70% across many strains), `SPRET_chr6_SINEB4` (SINE/B4, full pan-strain yet
absent from GRCm39). **+16 more** incl. a `C57BL_6NJ chrX L1` (a 6N-vs-6J substrain difference) and a non-TE cluster.

## Data dependencies (reads; not duplicated here)
Theme 22 `data/`: `graph_check_pav3.tsv` (graph PAV), `colo/*` (co-location lifts → per-strain coords + FPM),
`shared_subset_loci.csv` (loci). Plus `analysis/.../unique_pirna/cluster_pav/*.clusters_fpm.bed`, the
`results/STAR_srna_strain_wise` sRNA BAMs, RepeatMasker, and GFF (via `pav_clusters`).

## Run
`python make_nonref_locus.py <strain> <chrom> <start> <end> <graph_pav_name> <title> <outbase>` then the figure lands
in `figures/`. Batch: `python gen_nonref_loci_batch.py`. Use the env python that carries pysam (biomni_e1).
