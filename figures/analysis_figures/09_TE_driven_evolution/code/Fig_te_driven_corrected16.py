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
def _smart_labels(ax, xs, ys, labels, colors, fs=5.6, off=0.085):
    "greedy leader-arrow labels (log-log): place each label in the candidate direction farthest from points + placed labels"
    lx = np.log10(np.asarray(xs, float)); ly = np.log10(np.asarray(ys, float))
    xr = (lx.max() - lx.min()) or 1.0; yr = (ly.max() - ly.min()) or 1.0
    nx = lx / xr; ny = ly / yr
    dirs = [(np.cos(t), np.sin(t)) for t in np.linspace(0, 2*np.pi, 16, endpoint=False)]
    cand = [(d[0]*r, d[1]*r) for r in (off, off*1.7, off*2.5) for d in dirs]
    placed = []; pos = {}
    for i in sorted(range(len(nx)), key=lambda i: ny[i]):
        best = None; bs = -1
        for dx, dy in cand:
            cx, cy = nx[i] + dx, ny[i] + dy
            dmin = min([(cx-nx[j])**2 + (cy-ny[j])**2 for j in range(len(nx))] +
                       [(cx-p[0])**2 + (cy-p[1])**2 for p in placed] + [99])
            if dmin > bs: bs = dmin; best = (cx, cy)
        placed.append(best); pos[i] = best
    for i in range(len(nx)):
        cx, cy = pos[i]
        ax.annotate(labels[i], xy=(xs[i], ys[i]), xytext=(10**(cx*xr), 10**(cy*yr)),
                    fontsize=fs, color=colors[i], ha="center", va="center",
                    fontweight=("bold" if colors[i] == CW else "normal"),
                    arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.45, shrinkA=0, shrinkB=2.5), zorder=5)
# a: per-strain stacked creation + propagation
a = ax[0]; a.bar(x, df.est_creation, color=C_CRE, label="creation (new locus, breadth ≤3)")
a.bar(x, df.est_propagation, bottom=df.est_creation, color=C_PROP, label="propagation (insertion in conserved cluster, ≥10)")
xt(a); ns(a); a.set_ylabel("insertion-driven loci (n)", fontsize=8.5)
a.legend(fontsize=6.6, frameon=False, loc="upper left")
a.set_title(f"a   Insertion-driven loci split: creation vs propagation\noverall ≈{df.frac_creation.mean()*100:.0f}% creation / {df.frac_propagation.mean()*100:.0f}% propagation", fontsize=8.4, fontweight="bold", loc="left")
# b: creation vs insertion burden
b = ax[1]; b.scatter(df.private_ins, df.est_creation, c=cols, s=34, edgecolor="white", linewidth=0.5, zorder=3)
b.set_xscale("log"); b.set_yscale("log")
_smart_labels(b, df.private_ins.values, df.est_creation.values, [s.replace("_", "/") for s in df.strain], [CW if w else "#444" for w in df.wild])
b.margins(0.18)   # room for the offset leader-arrow labels
ns(b); b.set_xlabel("private insertions ≥150 bp (log)", fontsize=8.5); b.set_ylabel("CREATION loci (new strain-private, log)", fontsize=8.5)
b.set_title(f"b   True new loci (creation) still track insertion burden\nSpearman ρ = {rho_cre:.2f}, P = {pc_cre:.1e}", fontsize=8.4, fontweight="bold", loc="left")
# c: wild vs classical, creation & propagation (grouped)
cc = ax[2]; ww = 0.36
cc.bar(x - ww / 2, df.est_creation, ww, color=C_CRE, label="creation"); cc.bar(x + ww / 2, df.est_propagation, ww, color=C_PROP, label="propagation")
xt(cc); ns(cc); cc.set_yscale("log"); cc.set_ylabel("loci (n, log)", fontsize=8.5)
cc.legend(fontsize=6.8, frameon=False, loc="upper left")
cc.set_title(f"c   Wild-specificity holds for BOTH\ncreation P={p_cre:.1e} ({star(p_cre)}) · propagation P={p_prop:.1e} ({star(p_prop)})", fontsize=8.4, fontweight="bold", loc="left")
fig.suptitle("CORRECTED insertion-driven finding — most 'insertion-driven' loci are PROPAGATION (private insertion landed in a conserved cluster); a minority are CREATION (new locus), but creation still tracks insertion burden and stays wild-enriched", fontsize=9.4, fontweight="bold", y=1.02)
fig.tight_layout()
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/09_TE_driven_evolution/data/source_data"; _os.makedirs(_SD,exist_ok=True)
df.to_csv(f"{_SD}/SourceData_Fig_te_driven_corrected16.csv",index=False)   # per-strain creation vs propagation counts + private-insertion burden (all 3 panels)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_te_driven_corrected16.{e}", bbox_inches="tight")
print(f"wrote Fig_te_driven_corrected16 | creation: wild>cl P={p_cre:.2g}, ρ(burden)={rho_cre:.2f} | propagation: P={p_prop:.2g} | overall {df.frac_creation.mean()*100:.0f}% creation")
