#!/usr/bin/env python3
"""Developmental program of strain-specific TE-driven piRNA loci. Both INSERTION-driven and DIVERGENCE-driven loci
are enriched in the PREPACHYTENE program (E16.5, fetal TE-silencing / de-novo-methylation pathway) vs the all-
cluster baseline — i.e. strain-specific piRNA innovation acts on the fetal TE-silencing machinery, not the
pachytene (A-MYB) program. Dominant-timepoint distribution + Fisher E16.5-enrichment test."""
import sys, pandas as pd, numpy as np
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"); sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER
import pav_clusters as pc
from scipy.stats import fisher_exact
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"; CP = f"{U}/cluster_pav"; SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]; TPL = {"16.5dpc": "E16.5\n(prepachytene)", "12.5dpp": "P12.5", "20.5dpp": "P20.5\n(pachytene)"}; TPC = {"16.5dpc": "#0072B2", "12.5dpp": "#009E73", "20.5dpp": "#D55E00"}
ins_tp = []
for X in ORDER:
    d = pd.read_csv(f"{PG}/TE_driven_COORDINATE16_{X}.csv"); ins_tp += list(d.loc[d.groupby(["chrom", "start"]).FPM.idxmax()].tp)
P = pc._pang(); base_tp = list(P.loc[P.groupby(["strain", "own_chrom", "own_start"]).all_primary_FPM.idxmax()].tp)
dv = pd.read_csv(f"{PG}/divergence_candidates.tsv", sep="\t", dtype={"g39_chrom": str})
sup = pd.read_csv(f"{CP}/locus_genome_pav_divergence.tsv", sep="\t", dtype={"g39_chrom": str})
gn = sup[sup.present].groupby(["g39_chrom", "g39_start"]).size().rename("gn").reset_index()
dv = dv.merge(gn, on=["g39_chrom", "g39_start"]); dv = dv[dv.gn >= 10]
dv_tp = []
for _, r in dv.sample(min(150, len(dv)), random_state=1).iterrows():
    sub = pc.clusters_at(r.g39_chrom, int(r.g39_start), int(r.g39_end)); cc = sub[sub.strain == r.carrier]
    if len(cc): dv_tp.append(cc.loc[cc.all_primary_FPM.idxmax()].tp)
def cnt(lst): s = pd.Series(lst).value_counts(); return np.array([int(s.get(t, 0)) for t in TPS])
ins_c, base_c, dv_c = cnt(ins_tp), cnt(base_tp), cnt(dv_tp)
def fish(c): return fisher_exact([[c[0], c.sum() - c[0]], [base_c[0], base_c.sum() - base_c[0]]], alternative="greater")
ori, pi = fish(ins_c); ordv, pdv = fish(dv_c)
sets = {"insertion-driven": ins_c, "divergence-driven": dv_c, "baseline (all clusters)": base_c}
pd.DataFrame({k: v for k, v in sets.items()}, index=TPS).to_csv(f"{SD}/Fig_developmental_program16_sourcedata.csv")
print("E16.5 enrichment vs baseline: insertion OR=%.2f P=%.2g | divergence OR=%.2f P=%.2g" % (ori, pi, ordv, pdv))
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8})
fig, ax = plt.subplots(1, 1, figsize=(7.2, 4.6), dpi=300)
labels = list(sets); x = np.arange(len(labels)); w = 0.26
for j, t in enumerate(TPS):
    vals = [sets[k][j] / sets[k].sum() * 100 for k in labels]
    ax.bar(x + (j - 1) * w, vals, w, color=TPC[t], label=TPL[t].replace("\n", " "), edgecolor="white")
    for xi, v in zip(x + (j - 1) * w, vals): ax.text(xi, v + 1, f"{v:.0f}", ha="center", fontsize=6.5, color=TPC[t], fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([l.replace("-driven", "-\ndriven").replace(" (all clusters)", "\n(all clusters)") for l in labels], fontsize=8)
ax.set_ylabel("% of loci (dominant timepoint)", fontsize=9); ax.set_ylim(0, 98); ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=8)
ax.legend(fontsize=7.2, frameon=False, ncol=1, loc="upper right")
def star(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."
ax.text(0, 70, f"E16.5 vs baseline\nP={pi:.0e} ({star(pi)})", ha="center", fontsize=6.5, color="#0072B2")
ax.text(1, 88, f"E16.5 vs baseline\nP={pdv:.0e} ({star(pdv)})", ha="center", fontsize=6.5, color="#0072B2")
ax.set_title("Strain-specific TE-driven piRNA loci are PREPACHYTENE (E16.5, fetal TE-silencing)-biased —\nboth mechanisms enriched for E16.5 vs the all-cluster baseline, not the pachytene program", fontsize=8.6, fontweight="bold", loc="left", pad=10)
fig.tight_layout()
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_developmental_program16.{e}", bbox_inches="tight")
print("wrote Fig_developmental_program16.{png,pdf,svg} + source data")
