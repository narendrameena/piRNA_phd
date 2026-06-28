# Theme 24 — strain-private TE-insertion piRNA-cluster locus atlas

## Goal
The MIRROR of theme 23: per-locus figures (same style) for the top **strain-PRIVATE** non-reference piRNA clusters —
present in ONE strain only, absent from every other strain AND GRCm39. These are the young, strain-specific TE
insertions that **founded** a brand-new piRNA cluster (the leading edge of piRNA-cluster evolution).

## The figure
Same 3-panel figure as theme 23, via theme-23 `code/make_nonref_locus.py` (a documented cross-theme code dependency;
the generator routes `Fig_private_locus_*` output here). Panel A is the mirror of the shared subset: a single carrier
with FPM bars + graph coverage ≈ 1.0; **every other strain and GRCm39 flat at 0** ("present in 1/16").

## Loci (code/gen_private_loci_batch.py)
The 8 top INJECTED autosomal/X strain-private clusters by peak FPM, each **graph-confirmed carrier-only** (carrier ≈
1.0, max-other = 0, GRCm39 = 0): PWK chr2 (24.7 kb, 3440 FPM, LTR/ERVK), SPRET chr2/8/13, WSB chr17/8 — 3.4–24.7 kb,
86–3440 FPM. Mitochondrial "private" hits were excluded (mito is maternally conserved — a graph artifact, not a real
private locus).

## Causality (the honest key point)
The TE insertion **FOUNDS** the locus (origin): the piRNAs ARE the TE sequence and are antisense-to-TE (= silencing).
For a PRIVATE locus the insertion is the ONLY strain-specific difference, so it accounts for the cluster existing in
that one strain and nowhere else — the **clean causal case**. The insertion does NOT *set the expression level*: in the
SHARED subset (theme 23) the identical insertion is "present-but-silent" in some strains (a regulatory layer, shown in
panel A's tally). So neither atlas claims the insertion *determines expression* — only that it provides the locus.

## Limitation + a cautionary check (`code/25_dropped_locus_graph_pav.py`)
Only **50 / 865** autosomal strain-private loci are in the graph — the inject drops clusters that span fragment
boundaries OR fall in fragment GAPS (regions the minigraph-cactus graph doesn't represent), which biases toward the
LARGEST. **These dropped loci are NOT graph-confirmable and can be FALSE positives.** Worked example — the top one, the
"82-kb / 31,188-FPM WSB chr7 private cluster": it sits in an 846-kb WSB chr7 graph gap (3.01–3.86 Mb), so step 25 found
0 fragments → no graph PAV. A direct minimap2 of its 82 kb to GRCm39 returns **qcov = 1.00, mapq = 60 at GRCm39
chr7:6.57 Mb** — i.e. the sequence is fully in the reference at a structurally-shifted coordinate (likely a major
CONSERVED pachytene cluster). Its "private" call was a **halLiftover artifact** (divergent + graph-gapped region), not a
real private cluster. **Lesson:** the 8 graph-CONFIRMED loci above are the trustworthy set; any dropped/large "private"
locus must be minimap2-verified against GRCm39 before it is believed — graph confirmation is what filters these out.

## Data dependencies
Theme 22 `data/`: `private_te_loci.csv`, `graph_check_pav3.tsv` (with the appended private rows), `colo/*`. Plus
`unique_pirna/cluster_pav/*.clusters_fpm.bed`, `results/STAR_srna_strain_wise` BAMs, RepeatMasker, GFF. Theme 23
`code/make_nonref_locus.py` (shared generator).

## Run
`python code/gen_private_loci_batch.py` (biomni_e1 python with pysam). Locus selection: theme-22 `data/private_te_loci.csv`.
