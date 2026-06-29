#!/usr/bin/env python3
"""Confound fix step 1: classify insertion-driven loci by cluster expression breadth (clusters_at, authoritative).
creation = new strain-private locus (breadth<=3); propagation = private insertion LANDED in a conserved cluster
(breadth>=10). Per strain: estimate counts (sampled fraction x total) + test if wild-specificity holds for CREATION."""
import sys; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"); sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
import pav_clusters as pc, pandas as pd, numpy as np
from strain_order import STRAIN_ORDER, WILD
from scipy.stats import mannwhitneyu, spearmanr
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"; SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
P = pc._pang(); rows = []
for X in ORDER:
    d = pd.read_csv(f"{PG}/TE_driven_COORDINATE16_{X}.csv").drop_duplicates(["chrom", "start"]); total = len(d)
    samp = d.sample(min(150, len(d)), random_state=11); br = []
    for _, r in samp.iterrows():
        row = P[(P.strain == X) & (P.own_chrom == str(r.chrom)) & (P.own_start >= int(r.start) - 3000) & (P.own_start <= int(r.end) + 3000)]
        if len(row): br.append(pc.clusters_at(row.iloc[0].g39_chrom, int(row.start.min()), int(row.end.max())).strain.nunique())
    br = np.array(br); fc = (br <= 3).mean(); fp = (br >= 10).mean()
    rows.append((X, X in WILD, total, len(br), round(fc, 3), round(fp, 3), int(round(total * fc)), int(round(total * fp))))
df = pd.DataFrame(rows, columns=["strain", "wild", "total_insertion_loci", "n_sampled", "frac_creation", "frac_propagation", "est_creation", "est_propagation"])
# private-insertion burden (>=150bp) for the creation-vs-insertion-burden test
def n_ins150(X): return sum(1 for ln in open(f"{PG}/ins16/{X}.private_insertions.fasta") if ln[0] == ">" and int(ln.strip().rsplit("_", 1)[-1]) >= 150)
df["private_ins"] = [n_ins150(X) for X in df.strain]
df.to_csv(f"{SD}/Fig_insertion_creation_vs_propagation_sourcedata.csv", index=False)
print(df.to_string(index=False))
def mw(col): w = df[df.wild][col]; c = df[~df.wild][col]; return w.median(), c.median(), mannwhitneyu(w, c, alternative="greater").pvalue
for col in ["est_creation", "est_propagation", "frac_creation"]:
    wm, cm, p = mw(col); print(f"{col}: wild med={wm:.0f} classical med={cm:.0f} MW P(wild>cl)={p:.3g}")
print("Spearman est_creation ~ private_ins:", spearmanr(df.private_ins, df.est_creation))
print("overall mean fractions: creation %.0f%%  propagation %.0f%%" % (df.frac_creation.mean() * 100, df.frac_propagation.mean() * 100))
