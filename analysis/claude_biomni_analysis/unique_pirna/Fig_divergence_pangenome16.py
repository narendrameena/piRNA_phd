#!/usr/bin/env python3
"""Divergence-driven piRNA loci across the 16-strain pangenome — numbers + QC. Among the strongest strain-
restricted-CLUSTER loci (clusters_at breadth <=3, top 600 by FPM), most are GENOME-CONSERVED (present in many
strains) yet expressed in only 1-3 -> present-but-silent = divergence/regulatory, NOT genetic loss. Unlike the
insertion-driven set (wild-specific), divergence carriers span wild AND classical strains.
  (a) genome-presence (n strains) distribution -> peak at 16 = fully conserved;
  (b) expression breadth (how many strains express) -> mostly 1;
  (c) expression level (maxFPM) -> weak (vs strong insertion loci);
  (d) carriers per strain (wild vs classical) -> broad, not wild-specific;
  (e,f) QC: per-strain divergence-locus count vs sRNA depth + total PICB clusters (confounds)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"; CP = f"{U}/cluster_pav"
SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; CW, CC = "#C0392B", "#4393C3"
cand = pd.read_csv(f"{PG}/divergence_candidates.tsv", sep="\t", dtype={"g39_chrom": str})
sup = pd.read_csv(f"{CP}/locus_genome_pav_divergence.tsv", sep="\t", dtype={"g39_chrom": str})
gn = sup[sup.present].groupby(["g39_chrom", "g39_start"]).size().rename("genome_n").reset_index()
cand = cand.merge(gn, on=["g39_chrom", "g39_start"], how="left").fillna({"genome_n": 0})
cand["genome_n"] = cand.genome_n.astype(int)
div = cand[cand.genome_n >= 10].copy()                       # divergence-driven = genome-conserved + strain-restricted
div.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/10_cluster_pangenome_PAV/data/source_data/Fig_divergence_pangenome16_sourcedata.csv", index=False)
qc = pd.read_csv(f"{SD}/Fig_te_driven_pangenome16_QC_sourcedata.csv")   # per-strain depth + total clusters
percar = div.carrier.value_counts().reindex(ORDER).fillna(0).astype(int)
m = qc.set_index("strain"); m["div_loci"] = percar
rho_dep, p_dep = spearmanr(m.srna_depth, m.div_loci); rho_cl, p_cl = spearmanr(m.total_PICB_clusters, m.div_loci)
print(f"strain-restricted cluster loci={len(cand)} | genome-conserved (divergence)={len(div)} ({100*len(div)/len(cand):.0f}%)")
print(f"divergence carriers: wild={percar[[s for s in ORDER if s in WILD]].sum()} classical={percar[[s for s in ORDER if s not in WILD]].sum()}")
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8})
fig, ax = plt.subplots(2, 3, figsize=(13.5, 8.2), dpi=300)
def nostpr(a): a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=7)
# a: genome-presence distribution (all 600 strain-restricted cluster loci)
a = ax[0, 0]; a.hist(cand.genome_n, bins=np.arange(0.5, 17.5, 1), color="#7a8a99", edgecolor="white")
a.axvspan(9.5, 16.5, color=CC, alpha=0.10); a.set_xlabel("strains with the locus in genome (n/16)", fontsize=8.5); a.set_ylabel("strain-restricted cluster loci", fontsize=8.5); nostpr(a)
a.set_title(f"a   Most strain-restricted loci are genome-CONSERVED\n{len(div)}/{len(cand)} ({100*len(div)/len(cand):.0f}%) present in ≥10 strains → divergence", fontsize=8.4, fontweight="bold", loc="left")
# b: expression breadth
b = ax[0, 1]; vc = div.breadth.value_counts().reindex([1, 2, 3]).fillna(0); b.bar([1, 2, 3], vc.values, color="#6a3d9a", edgecolor="white", width=0.6)
b.set_xticks([1, 2, 3]); b.set_xlabel("strains EXPRESSING the locus (n)", fontsize=8.5); b.set_ylabel("divergence loci", fontsize=8.5); nostpr(b)
b.set_title("b   Expression is strain-restricted\n(mostly a single strain)", fontsize=8.4, fontweight="bold", loc="left")
# c: expression level (weak)
c = ax[0, 2]; c.hist(np.log10(div.maxFPM), bins=20, color="#B8860B", edgecolor="white")
c.set_xlabel("max PICB FPM (log10)", fontsize=8.5); c.set_ylabel("divergence loci", fontsize=8.5); nostpr(c)
c.set_title("c   Divergence loci are WEAKLY expressed\n(present-but-near-silent)", fontsize=8.4, fontweight="bold", loc="left")
# d: carriers per strain (wild vs classical)
d = ax[1, 0]; cols = [CW if s in WILD else CC for s in ORDER]
d.bar(range(len(ORDER)), percar.values, color=cols, edgecolor="white", linewidth=0.4)
d.set_xticks(range(len(ORDER))); d.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=6.3)
for t, s in zip(d.get_xticklabels(), ORDER):
    if s in WILD: t.set_color(CW); t.set_fontweight("bold")
d.set_ylabel("divergence loci (carrier strain)", fontsize=8.5); d.set_ylim(0, percar.max() * 1.22); nostpr(d)
d.set_title("d   Divergence carriers span wild AND classical\n(NOT wild-specific, unlike insertion-driven)", fontsize=8.4, fontweight="bold", loc="left")
d.legend(handles=[Patch(facecolor=CW, label="wild"), Patch(facecolor=CC, label="classical")], fontsize=6.8, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.98))
# e: QC depth
e = ax[1, 1]; e.scatter(m.srna_depth / 1e9, m.div_loci, c=cols, s=34, edgecolor="white", linewidth=0.5, zorder=3)
nostpr(e); e.set_xlabel("sRNA library depth (billion mapped reads)", fontsize=8.5); e.set_ylabel("divergence loci per strain", fontsize=8.5)
e.set_title(f"e   QC: vs sRNA depth\nSpearman ρ={rho_dep:.2f}, P={p_dep:.2f}", fontsize=8.4, fontweight="bold", loc="left")
# f: QC total clusters
fa = ax[1, 2]; fa.scatter(m.total_PICB_clusters / 1000, m.div_loci, c=cols, s=34, edgecolor="white", linewidth=0.5, zorder=3)
nostpr(fa); fa.set_xlabel("total PICB clusters per strain (×1,000)", fontsize=8.5); fa.set_ylabel("divergence loci per strain", fontsize=8.5)
fa.set_title(f"f   QC: vs total cluster count\nSpearman ρ={rho_cl:.2f}, P={p_cl:.2f}", fontsize=8.4, fontweight="bold", loc="left")
fig.suptitle("Divergence-driven piRNA loci across the 16-strain pangenome — genome-conserved but strain-restricted expression (present-but-silent), spanning wild and classical strains", fontsize=10.5, fontweight="bold", y=1.005)
fig.tight_layout()
for ext in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_divergence_pangenome16.{ext}", bbox_inches="tight")
print(f"wrote Fig_divergence_pangenome16.{{png,pdf,svg}} (rho_depth={rho_dep:.2f}, rho_clusters={rho_cl:.2f})")
