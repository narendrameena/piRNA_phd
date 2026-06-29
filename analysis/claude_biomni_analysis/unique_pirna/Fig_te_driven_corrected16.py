#!/usr/bin/env python3
"""CORRECTED insertion-driven finding (confound-fix): split insertion-driven loci into CREATION (new strain-private
locus, cluster breadth<=3) vs PROPAGATION (private insertion LANDED in a conserved cluster, breadth>=10). The chr17
hotspot etc. are propagation. (a) per-strain creation vs propagation; (b) does CREATION still track private-insertion
burden? (c) is CREATION wild-specific?  -> separates the true new-locus signal from insertions hitting conserved clusters."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
df = pd.read_csv(f"{SD}/Fig_insertion_creation_vs_propagation_sourcedata.csv")
CW, CC, C_CRE, C_PROP = "#C0392B", "#4393C3", "#6a3d9a", "#c9b3e0"
cols = [CW if w else CC for w in df.wild]; x = np.arange(len(df)); xl = [s.replace("_", "/") for s in df.strain]
w, c = df[df.wild], df[~df.wild]
def mwp(col): return mannwhitneyu(w[col], c[col], alternative="greater").pvalue
p_cre, p_prop = mwp("est_creation"), mwp("est_propagation")
rho_cre, pc_cre = spearmanr(df.private_ins, df.est_creation); rho_prop, pc_prop = spearmanr(df.private_ins, df.est_propagation)
def star(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8})
fig, ax = plt.subplots(1, 3, figsize=(14, 4.7), dpi=300)
def xt(a): a.set_xticks(x); a.set_xticklabels(xl, rotation=90, fontsize=6.3); [t.set_color(CW) or t.set_fontweight("bold") for t, wd in zip(a.get_xticklabels(), df.wild) if wd]
def ns(a): a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=7)
# a: per-strain stacked creation + propagation
a = ax[0]; a.bar(x, df.est_creation, color=C_CRE, label="creation (new locus, breadth ≤3)")
a.bar(x, df.est_propagation, bottom=df.est_creation, color=C_PROP, label="propagation (insertion in conserved cluster, ≥10)")
xt(a); ns(a); a.set_ylabel("insertion-driven loci (n)", fontsize=8.5)
a.legend(fontsize=6.6, frameon=False, loc="upper left")
a.set_title(f"a   Insertion-driven loci split: creation vs propagation\noverall ≈{df.frac_creation.mean()*100:.0f}% creation / {df.frac_propagation.mean()*100:.0f}% propagation", fontsize=8.4, fontweight="bold", loc="left")
# b: creation vs insertion burden
b = ax[1]; b.scatter(df.private_ins, df.est_creation, c=cols, s=34, edgecolor="white", linewidth=0.5, zorder=3)
for _, r in df.iterrows(): b.annotate(r.strain.replace("_", "/"), (r.private_ins, r.est_creation), fontsize=5, color=(CW if r.wild else "#666"), xytext=(3, 1), textcoords="offset points")
b.set_xscale("log"); b.set_yscale("log"); ns(b); b.set_xlabel("private insertions ≥150 bp (log)", fontsize=8.5); b.set_ylabel("CREATION loci (new strain-private, log)", fontsize=8.5)
b.set_title(f"b   True new loci (creation) still track insertion burden\nSpearman ρ = {rho_cre:.2f}, P = {pc_cre:.1e}", fontsize=8.4, fontweight="bold", loc="left")
# c: wild vs classical, creation & propagation (grouped)
cc = ax[2]; ww = 0.36
cc.bar(x - ww / 2, df.est_creation, ww, color=C_CRE, label="creation"); cc.bar(x + ww / 2, df.est_propagation, ww, color=C_PROP, label="propagation")
xt(cc); ns(cc); cc.set_yscale("log"); cc.set_ylabel("loci (n, log)", fontsize=8.5)
cc.legend(fontsize=6.8, frameon=False, loc="upper left")
cc.set_title(f"c   Wild-specificity holds for BOTH\ncreation P={p_cre:.1e} ({star(p_cre)}) · propagation P={p_prop:.1e} ({star(p_prop)})", fontsize=8.4, fontweight="bold", loc="left")
fig.suptitle("CORRECTED insertion-driven finding — most 'insertion-driven' loci are PROPAGATION (private insertion landed in a conserved cluster); a minority are CREATION (new locus), but creation still tracks insertion burden and stays wild-enriched", fontsize=9.4, fontweight="bold", y=1.02)
fig.tight_layout()
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_te_driven_corrected16.{e}", bbox_inches="tight")
print(f"wrote Fig_te_driven_corrected16 | creation: wild>cl P={p_cre:.2g}, ρ(burden)={rho_cre:.2f} | propagation: P={p_prop:.2g} | overall {df.frac_creation.mean()*100:.0f}% creation")
