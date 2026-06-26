# Theme 21 — piRNA-cluster PAV: pangenome GRAPH (odgi) vs HAL-LIFTOVER

## Goal
Implement the pangenome graph-based piRNA-cluster presence/absence (PAV) classification (after the navin
`pirna_cluster_pangenome_comparison.sh` concept) on the project's built minigraph-cactus graph, and **compare it
against the existing HAL-liftover PAV** (`cluster_pav/`).

## Biology (grounded)
piCB piRNA clusters are called per strain **from small-RNA reads**. To ask *"is a cluster conserved or strain-variable
across the 16 strains?"* the two methods differ fundamentally:
- **HAL-liftover** projects every strain's clusters onto **one** reference (GRCm39) via the cactus HAL. A cluster in
  non-reference / TE-insertion sequence **cannot lift** → it is mis-called private/absent (reference bias), and
  divergent strains' clusters **fragment** heavily under lift (e.g. SPRET 18 k native clusters → 145 k lifted intervals).
- **Pangenome graph (`odgi pav`)** carries every strain's sequence as nodes. At a GRCm39 cluster locus, the per-strain
  PAV ratio = fraction of the locus's graph nodes a strain traverses = **genetic sequence presence**.
- The liftover PAV measures **cluster (regulatory) presence**; the graph PAV measures **sequence (genetic) presence**.
  Their disagreement separates **genetic loss/gain (TE insertion)** from **epigenetic silencing** (sequence present,
  no cluster) — the central axis of piRNA-cluster evolution.

## Inputs (all verified; three coordinate conventions reconciled)
- Graph: `results/pangenome/output/mouse_17strain_pangenome.gfa.gz` (final; for `odgi build`) / `.gbz` (vg) /
  `.full.hal`. Samples = **GRCm39 + 16 strains** (17). Paths: GRCm39 = `GRCm39#0#<chr>` (one path/chrom); strains =
  `STRAIN#0#<chr>#<offset>` (**fragmented** contig paths). Tools `odgi` v0.9.0 + `vg` v1.61.0 in `cactus_v2.9.3.sif`.
- Precomputed graph products (reused, no recompute): `mouse_17strain_pangenome.snarls`, `...vcf.gz`
  (vg deconstruct vs GRCm39, 16 genotypes) — for SV/TE step.
- Clusters (piCB, from sRNA): `cluster_pav/{strain}.clusters.bed` (native, bare chrom) +
  `{strain}.GRCm39.merged.bed` (lifted) + `bytp/{strain}.{tp}.clusters.bed` (per-tp; **all 3 tps incl. fetal E16.5**).
- Liftover PAV baseline (existing): `cluster_pav/cluster_PAV_catalogue.csv.gz`, `locus_genome_pav.tsv`.
- Bulk/whole-testis RNA (readthrough, stranded RF): `results/STAR_rna_strain_wise/{strain}-{tp}.{rep}/` (`STRAIN#1#chr<N>`).
- Gene models: `resources/annotation/{strain}_v3.3.gff3` (78 files). Genome FASTA: `results/pangenome/prepared/{strain}.fa`.
- RepeatMasker `.out` **precomputed** (21): `resources/repeatMasker/{strain}_*.out` (no RM run). `gffread` (conda) +
  system `blastn` for pseudogene-fragment step. Big-mem SLURM: TEST (1 TB), NXFL/2004/2204 (385 GB).

## Key corrections found (verified, NOT assumed)
1. The concept script's `odgi depth --path-name --path-group` flags **do not exist** in odgi v0.9.0 → use the dedicated
   **`odgi pav`** subcommand (`-S` group-by-sample, `-M` matrix).
2. Anchoring `odgi pav` on **GRCm39** = sequence-PAV (genetic), the biologically right complement to the liftover
   cluster-PAV; anchoring on fragmented strain paths is reserved for non-lifting (graph-only) clusters.
3. Snarls + deconstructed VCF already exist → reuse for the SV/TE step.

## Steps (code/)
1. `01_prep_master_loci.sh` — 42,384 master loci (GRCm39 frame, union of lifted clusters); **liftover cluster-PAV**
   (locus×16) → **core 5,983 / dispensable 22,866 / private 13,535**; odgi-pav BED (`GRCm39#0#chrom`).  [DONE]
2. `09_rebuild_plines_pav.slurm` — graph PAV. **CRITICAL FIX:** `odgi build` imports GFAv1 **P-lines only**, but the
   cactus GFA stores paths as **W-lines** → the first build (`02_…`) gave a graph with **0 paths** (pav impossible).
   Fix = `vg convert -W` (gbz → P-line GFA) → `odgi build` (6,442 paths) → `odgi pav -S -M` at master loci → graph
   sequence-PAV matrix (validated on a 2 Mb chunk first).  [DONE]
3. `03_compare_graph_vs_liftover.py` — graph PAV ratio ≥0.5 = sequence present; per-locus graph vs liftover class;
   disagreements: liftover-absent & graph-present = **silencing**; both absent = **loss**.  [DONE]
4. Downstream: `05_sv_te_variable_loci.sh` (SV/TE from precomputed VCF + TE consensus blastn); `04_devclass_genecontext.sh`
   (dev-class from per-tp lifted clusters incl. E16.5 + gene context vs GRCm39); `07_pgf_pachytene.sh` (antisense
   pseudogene fragments); `06_readthrough.slurm` (piC-DoG readthrough from whole-testis bulk RNA).  [DONE]
5. Figures `Fig_graph_vs_liftover_pav.py` + `Fig_pangenome_cluster_drivers.py` + `.note.md`.  [DONE]

## Key corrections found (verified, NOT assumed)
4. `odgi build` reads **P-lines only**; cactus GFA = **W-lines** → must `vg convert -W` first, else 0 paths.

## Results
- **Graph vs liftover (Fig 1):** 99 % of liftover-"private" clusters are sequence-**SHARED** in the graph; 52–57 % of
  strain×locus events = **silencing** (sequence+, cluster−) vs only 3–8 % genetic loss; median strains: graph-sequence
  16 vs liftover-cluster 4. → strain-variable piRNA clusters are mostly **epigenetic**, not genetic; the graph separates
  silencing from loss (liftover conflates). Robust across PAV thresholds 0.5–0.95.
- **Drivers (Fig 2):** genetic-novelty clusters = young **L1MdTf/L1MdA + IAP** insertions; pachytene → intergenic
  (A-MYB), pre-pachytene/hybrid → 3′UTR-enriched; **1,285 antisense** pseudogene fragments; **~65 %** piC-DoG readthrough.

## Status
COMPLETE — graph PAV (P-line-fixed), comparison, downstream biology, both figures + notes done. Ready to commit.
