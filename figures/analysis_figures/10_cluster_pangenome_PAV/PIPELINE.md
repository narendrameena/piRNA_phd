# 10 — piRNA-cluster pangenome (presence/absence) + the two evolutionary mechanisms

**What these figures are.** The open piRNA-cluster pangenome across 16 strains (core/accessory/private),
expression concentration, Zamore stage conservation, and the decisive test splitting strain-restricted
clusters into **insertion-driven** vs **divergence-driven**. The two-mechanism split was **BioMNI
triple-verified** (memory `project_te_driven_finding`, 2026-06-15).

---

## STEP-BY-STEP (tool · version · parameters · result)

**S1 · clusters → pangenome frame.** PICB (R 4.2.3) combined clusters per strain projected to GRCm39 via
**cactus v2.9.3 halLiftover** (`build_picb_pangenome_fpm.py` → `picb_pangenome_clusters.tsv`, 2.89 M rows;
own + GRCm39 coords + FPM). **No pairwise liftover.**

**S2 · cluster PAV (core/accessory/private).** **bedtools 2.31.1** merge projected clusters → per-GRCm39
cluster, which strains carry/express it → openness/rarefaction. **Result #:** open pangenome (`Fig_pirna_pangenome16`).

**S3 · expression concentration.** **Python 3.11.15** Lorenz curve / Gini across clusters. → `Fig_pirna_lorenz16`.

**S4 · Zamore stage conservation.** pachytene-cluster (Zamore, **mm10→mm39 UCSC liftOver**) conservation
across strains. → `Fig_zamore_stage_conservation`.

**S5 · two-mechanism test [finding, BioMNI-verified].** **Python/scipy** —
**insertion-driven** (TE inserts → NEW locus, absent elsewhere) = wild-specific, #loci ∝ insertion burden
**ρ ≈ 0.94**; **divergence-driven** (conserved locus, present-but-silent, activated in some strains) = broad,
**~83 %** of strain-restricted-cluster loci. → `Fig_te_driven_pangenome16`, `Fig_divergence_pangenome16`.

**S6 · figures.** matplotlib (Python 3.11.15). → **6** figures.

## TOOLS
| Tool | Version | What/why | Key params |
|---|---|---|---|
| PICB | R 4.2.3 | combined clusters | `clusters` sheet |
| cactus / halLiftover | v2.9.3 | project clusters → GRCm39 (pangenome) | strain→GRCm39 |
| bedtools | 2.31.1 | cluster PAV merge/overlap | `merge`/`intersect` |
| UCSC liftOver | 2023 | Zamore mm10→mm39 (only UCSC use) | `mm10ToMm39.over.chain.gz` |
| Python / scipy | 3.11.15 | Gini, ρ, two-mechanism split, figures | Spearman; matplotlib |

## INPUTS  `unique_pirna/cluster_pav/picb_pangenome_clusters.tsv`; `cluster_PAV_catalogue.csv.gz`; Zamore loci → `data/SourceData_pirna_pangenome16.csv`.
## OUTPUTS (`figures/`, 6)  `Fig_pangenome_pav`, `Fig_pirna_pangenome16`, `Fig_pirna_lorenz16`, `Fig_zamore_stage_conservation`, `Fig_te_driven_pangenome16`, `Fig_divergence_pangenome16`.

## DOUBLE-VERIFICATION
- Coordinate/pangenome-based (NOT sequence-containment) — the proxy was replaced after the 2026-06-15 confound check.
- Two-mechanism split **BioMNI triple-verified**; cross-ref themes 09 (creation) and 13 (divergence).
- mm10→mm39 liftover enforced for Zamore loci.
