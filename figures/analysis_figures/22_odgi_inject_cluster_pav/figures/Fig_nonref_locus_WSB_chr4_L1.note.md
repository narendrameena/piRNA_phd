# Fig_nonref_locus_WSB_chr4_L1 — a major piRNA cluster genetically absent from GRCm39

A per-locus view (in the style of `make_pav_locus`) of the headline theme-22 non-reference finding, with the pangenome
GRAPH as the cross-strain evidence (odgi inject + pav).

**The locus:** WSB/EiJ chr4:38,687,900–38,698,300 — a **LINE/L1-driven piRNA cluster peaking at P20.5 (213 FPM all-read;
247 unique-read)**, the highest-expression cluster in the 'present-in-most-but-not-GRCm39' subset.

- **(A) GRAPH PAV (odgi inject + pav — the odgi analysis):** the cluster's sequence is present in **10/16 strains**
  (graph coverage ≥0.5) but **GRCm39 = 0** — genetically absent from the reference path (with the validated controls of
  step 20/21). A major piRNA cluster the single reference entirely lacks.
- **(B)** WSB per-timepoint sRNA coverage (height = expression, colour = timepoint) over TE + gene tracks — the L1, the
  P20.5 expression peak, architecture, and per-tp TE/silencing summary.
- **(C)** base resolution — individual piRNAs at their true coordinates; 5′ arrow red = antisense-to-TE (silencing).

Built with `make_nonref_locus.py`, which reuses the strain-anchored `pav_clusters` drawing blocks for B/C (they need no
GRCm39 anchor) and makes panel A graph-native. Re-runnable for any non-reference locus present in `graph_check_pav3.tsv`:
`make_nonref_locus.py <strain> <chrom> <start> <end> <graph_pav_name> <title> <outbase>` (run with the biomni_e1 python,
which has pysam). Run with the env's python that carries pysam.
