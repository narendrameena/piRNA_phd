# Fig_shared_subset — the 'present-in-most-but-not-GRCm39' subset (mostly a liftover artifact)

The non-reference clusters shared with >=14 other strains (present in 15-16/16 strains, yet flagged non-reference
because they did not halLiftover to GRCm39). Characterised on three axes:

- **(A) Collapse** (`17`): 117 strain-entries -> **22 distinct loci** (co-location graph), one present in all 16 strains.
- **(B) TE profile** (`17`): distinct from the bulk non-reference set — **SINE/B4-enriched** (vs ERVK/L1) and **younger**
  (median RepeatMasker divergence 6.7% vs 14.3% bulk non-ref, 16.4% conserved).
- **(C/D) Why absent from GRCm39** (`18`; minimap2 of the unique flanks + the TE-rich cluster body to GRCm39 + N-check):
  - **15/22 are PRESENT in GRCm39** (cluster qcov ~1.0 at the homologous locus) — the TE-rich cluster body simply
    **failed to halLiftover** = FALSE non-reference (a liftover limitation for repetitive clusters). This includes the
    high-expression WSB chr4 L1 (247 FPM) — it IS in the reference; the earlier 'reference misses it' was wrong.
  - **6/22 are genuine C57BL/6J-lineage DELETIONS** — flanks lift, cluster body absent, bracket clean (0% N): the panel
    kept a cluster the reference lineage lost (e.g. FVB chr1 19kb ERVK 7 FPM; PWK chr1 9kb L1 6 FPM).
  - 1 divergent/absent, **0 assembly gaps** (GRCm39 is gap-free here).

**Implication:** the 'shared-but-not-GRCm39' subset (~8% of the non-reference entries) is **mostly a liftover artifact**
for TE-rich clusters that are actually present in GRCm39; the graph-native odgi-inject PAV (theme 22) sidesteps it.
The strain-PRIVATE bulk (66%, genuine young post-divergence insertions) is unaffected — this only refines the small
shared tail. A real minority (6 loci) are genuine C57BL/6J-lineage losses.
