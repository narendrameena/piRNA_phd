#!/usr/bin/env python3
"""QC + statistics for the pangenome insertion-driven piRNA finding. Tests whether 'wild-derived strains gain
more insertion-driven piRNA source loci' is biological (tracks private-insertion burden) or a confound of
(i) sRNA library depth, (ii) total PICB-cluster count, or (iii) reference-bias inflation of private insertions.
Outputs per-strain QC table + stats."""
import os, sys, numpy as np, pandas as pd, pysam
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER, WILD
from scipy.stats import spearmanr, mannwhitneyu
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
PG = f"{U}/pangenome_te"; CP = f"{U}/cluster_pav"; RB = f"{ROOT}/results/STAR_srna_strain_wise"; SD = f"{ROOT}/analysis/claude_biomni_analysis/source_data"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
def depth(X):
    t = 0
    for tp in TPS:
        for r in (1, 2, 3):
            b = f"{RB}/{X}/{X}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
            if os.path.exists(b):
                try: t += sum(s.mapped for s in pysam.AlignmentFile(b, "rb").get_index_statistics())
                except Exception: pass
    return t
def nclusters(X):
    cl = set()
    for ln in open(f"{CP}/{X}.clusters_fpm.bed"):
        f = ln.split("\t"); cl.add((f[0], f[1], f[2]))
    return len(cl)
def n_ins150(X):
    return sum(1 for ln in open(f"{PG}/ins16/{X}.private_insertions.fasta") if ln[0] == ">" and int(ln.strip().rsplit("_", 1)[-1]) >= 150)
def n_loci(X):
    return pd.read_csv(f"{PG}/TE_driven_COORDINATE16_{X}.csv").drop_duplicates(["chrom", "start"]).shape[0]
rows = [(X, X in WILD, n_ins150(X), n_loci(X), nclusters(X), depth(X)) for X in ORDER]
df = pd.DataFrame(rows, columns=["strain", "wild", "private_ins_ge150bp", "insertion_driven_loci", "total_PICB_clusters", "srna_depth"])
df["frac_clusters_ins_driven"] = df.insertion_driven_loci / df.total_PICB_clusters
df["loci_per_1k_insertions"] = df.insertion_driven_loci / df.private_ins_ge150bp * 1000
df.to_csv(f"{SD}/Fig_te_driven_pangenome16_QC_sourcedata.csv", index=False)
print(df.to_string(index=False)); print()
def mw(col):
    w = df[df.wild][col]; c = df[~df.wild][col]; u, p = mannwhitneyu(w, c, alternative="greater")
    return w.median(), c.median(), w.min(), w.max(), c.min(), c.max(), u, p
print("=== Mann-Whitney U (wild > classical, one-sided; n=4 wild vs 12 classical) ===")
for col in ["insertion_driven_loci", "frac_clusters_ins_driven", "loci_per_1k_insertions", "private_ins_ge150bp", "srna_depth", "total_PICB_clusters"]:
    wm, cm, wmn, wmx, cmn, cmx, u, p = mw(col)
    print(f"  {col:28s}: wild med={wm:.4g} [{wmn:.3g}-{wmx:.3g}]  classical med={cm:.4g} [{cmn:.3g}-{cmx:.3g}]  U={u:.0f}  p={p:.4g}")
print("\n=== Spearman correlations (what does loci count track?) ===")
for x in ["private_ins_ge150bp", "srna_depth", "total_PICB_clusters"]:
    r, p = spearmanr(df[x], df.insertion_driven_loci); print(f"  loci ~ {x:24s}: rho={r:.3f}  p={p:.3g}")
r, p = spearmanr(df.private_ins_ge150bp, df.frac_clusters_ins_driven); print(f"  frac_clusters ~ private_ins  : rho={r:.3f}  p={p:.3g}")
r, p = spearmanr(df.srna_depth, df.private_ins_ge150bp); print(f"  depth ~ private_ins (confound?): rho={r:.3f}  p={p:.3g}")
print("\nwrote QC source data")
