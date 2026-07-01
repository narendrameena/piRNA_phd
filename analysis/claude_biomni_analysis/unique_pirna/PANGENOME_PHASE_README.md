# Scale-16 pangenome phase — cross-strain expression test (the 16-strain Step-4)

**Goal.** For every strain-specific candidate piRNA from the 16-strain edgeR DA (`edger16/{tp}.strain_specific_DA.csv.gz`), decide whether it is **genuinely unique by expression** or merely looks strain-specific — the 16-strain equivalent of the pilot Step 4 — **without** 16×15 pairwise genome mapping. We split into the same 4 classes used in the pilot:

| class | meaning | unique? |
|---|---|---|
| `expressed elsewhere (exact)` | exact sequence is in another strain's expressed pool | NO |
| `unique: conserved-but-silent` | novel sequence; orthologous locus EXISTS in other strains but is not expressed there | YES (expression) |
| `unique: strain-private locus` | novel sequence; orthologous locus exists in NO other strain | YES (locus gain) |
| `SNP-variant (1–3 mm)` *(refinement)* | a ≤3-mm allele of the candidate is expressed at the orthologous locus in another strain | NO |

**Why this is "via the pangenome, not pairwise."** The expression question is answered at the **sequence** level (exact match against each strain's expressed pool — build-independent, unambiguous). The locus-homology question (does an orthologous locus exist elsewhere?) is answered by projecting each candidate's production locus into the **GRCm39 common frame through the minigraph-cactus pangenome HAL** (`halLiftover`, the same SIF + graph that produced the cluster PAV), then lifting back out to each strain — 16 + 16 liftovers total, replacing 240 pairwise genome alignments.

## Inputs (dependencies)
- `edger16/{tp}.strain_specific_DA.csv.gz` — 16-strain strain-specific candidates (job **3302624**, chained after matrices **3302621**).
- `unique16/pools/{strain}.pool.txt.gz` — per-strain expressed pools, 24–32 nt, ≥2 reads (job **3302627**, running now).
- pangenome: `cactus_v2.9.3.sif` + `results/pangenome/output/mouse_17strain_pangenome.full.hal` (HAL genome names = strain names + `GRCm39`; chroms = Ensembl `1..19,X,Y,MT`).
- per-strain STAR indices `results/indexs/{strain}` (PanSN chroms `{strain}#1#chrN`).

## Run order
```
# 0. (running) pools           : sbatch run_pools16.sh                 -> unique16/pools/{strain}.pool.txt.gz   [job 3302627]
# 0. (chained) matrices+edgeR  : 3302621 -> 3302624                    -> edger16/{tp}.strain_specific_DA.csv.gz

# 1. exact cross-strain expression (needs pools + DA). Big-mem node (loads 16 pools as a bitmask, ~30-50 GB):
sbatch --mem=240G --partition=2004,TEST --wrap \
  '/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python classify_unique16.py'
#    -> unique16/expr_exact_classified.csv.gz  +  per-strain/tp  unique16/{X}.{tp}.novel.fasta

# 2a. map NOVEL candidates to own genome (validated sRNA piRNA STAR params):
sbatch cand_loci16.sh            # array 0-15 -> unique16/loci/{X}.cand_loci.ens.bed
# 2b. lift candidate loci  X -> GRCm39  (through the pangenome HAL):
sbatch --dependency=afterok:<2a> lift_cand16.sh   # array 0-15 -> unique16/loci/{X}.cand_GRCm39.bed
# 2c. build the union, then lift GRCm39 -> each strain Y for cross-strain presence:
cat unique16/loci/*.cand_GRCm39.bed > unique16/loci/all.cand_GRCm39.bed
sbatch lift_presence16.sh        # array 0-15 -> unique16/loci/present_in_{Y}.bed

# 3. final 4-way classification:
/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python classify_unique16_locus.py
#    -> unique16/final_classified.csv.gz   (klass + homolog_strains, per strain x timepoint)
```

## Design decisions & caveats (reviewer-facing)
1. **Expression = sequence in the strain's expressed pool** (pilot definition). Pool floor = ≥2 reads (non-singleton), pooled over a strain's 9 collapse files; built from collapse (not the min.total≥15 union matrix) so single-strain sequences are not lost.
2. **Exact (0-mm) expression is build-independent** and needs no genome/pangenome — done first, removes the bulk of non-unique candidates unambiguously.
3. **Locus homology via halLiftover** (coordinate orthology through the cactus graph). Same caveat as the cluster PAV: halLiftover fragments long intervals, but candidate piRNA loci are short (≤32 bp + the bamtobed interval), so fragmentation is minimal; presence = "any lifted row for this candidate id in strain Y". Chrom names reconciled PanSN→Ensembl exactly as the cluster PAV (`sed s/^{X}#1#chr//`).
4. **TE-multimapping:** a repetitive candidate maps to many loci; presence is the UNION over all its loci (homolog present if *any* locus exists elsewhere) — conservative against over-calling strain-private.
5. **SNP-variant (1–3 mm) refinement is a documented add-on, not yet wired in.** First pass folds "novel sequence + homolog present elsewhere" into `conserved-but-silent`. To split out true SNP-variants: at strain Y's orthologous locus, reverse-lift GRCm39→Y, fetch Y's genomic allele (needs the per-strain unmasked FASTA), compare ≤3 mm to the candidate (Poisson cutoff k*=3 from Step 5) and check it is in Y's pool. This slightly inflates `conserved-but-silent` at the expense of `SNP-variant` until added — flagged so it is not over-interpreted.
6. **DRY-RUN before the heavy run:** verify `cand_loci16.sh` STAR maps a small novel FASTA cleanly (read names with `|` survive bamtobed) and that `lift_cand16.sh` produces non-empty `cand_GRCm39.bed` on one strain, before launching all 16 — the STAR params/paths are validated but the candidate-FASTA-as-reads path is new for this project.

## Verification status of the biology this rests on (VERIFICATION_QUEUE.md)
- "Uniqueness = expression, not DNA presence" and the TE-driven-locus premise: BioMNI-verified (genomics + general).
- M3–M6 drained 2026-06-11: M3/M4/M5 VERIFIED; M6 verified-with-correction (mapping rationale = exact genomic origin + repeat multimapping).
- All agent citations remain PROVISIONAL pending EuropePMC.
