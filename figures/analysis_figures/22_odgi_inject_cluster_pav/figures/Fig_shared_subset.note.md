# Fig_shared_subset — the 'present-in-most-but-not-GRCm39' subset (GRAPH-CONFIRMED genuinely absent)

The non-reference clusters shared with >=14 other strains (present in 15-16/16 strains, flagged non-reference because
they did not halLiftover to GRCm39). Characterised on three axes:

- **(A) Collapse** (`17`): 117 strain-entries -> **22 distinct loci** (co-location graph), one present in all 16 strains.
- **(B) TE profile** (`17`): distinct from the bulk non-reference set — **SINE/B4-enriched** (vs ERVK/L1) and **younger**
  (median RepeatMasker divergence 6.7% vs 14.3% bulk non-ref, 16.4% conserved).
- **(C/D) Why absent from GRCm39 — settled by the GRAPH** (`19`/`20`/`21`; `odgi inject` + `odgi pav` on graph_inj.og,
  the actual pangenome MSA — no halLiftover, no minimap2). GRCm39's coverage of each cluster's nodes, with controls:
  - **GRCm39-frame controls** (regions in GRCm39's own coords) -> **1.000**; **C57BL_6NJ lifted controls** (known-present)
    -> **1.000** => the GRCm39 group is correct.
  - the **19 shared-subset loci -> 0.000** => GRCm39 genuinely does NOT traverse them = **genuinely absent from the
    reference path** (present in the strains and in C57BL_6NJ, missing from C57BL/6J/GRCm39).

**Correction:** an earlier minimap2 pass (`18`) called 15/22 'present in GRCm39 (technical lift-artifact)'. The graph
shows that was the **TE-sequence-matches-elsewhere confound** — the TE *sequence* has a copy in GRCm39, but the
*insertion at this locus* is not on GRCm39's path. So these are **NOT liftover artifacts**: they are genuine
**C57BL/6J-lineage absences** — piRNA clusters most strains carry but the single reference lacks, including the
high-expression WSB chr4 L1 (247 FPM). The whole non-reference catch is genuine; the graph-native inject PAV is the
reliable arbiter (it is the MSA, not a heuristic). A single reference under-represents real piRNA clusters.
