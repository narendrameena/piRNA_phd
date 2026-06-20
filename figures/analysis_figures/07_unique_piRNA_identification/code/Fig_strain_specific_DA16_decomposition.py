#!/usr/bin/env python3
"""Filter decomposition of the strain-specific piRNA call (companion to Fig_strain_specific_DA16), 16 strains x 3
timepoints, each filter SEPARATELY:
  A  edgeR DA only        (FDR<0.05 & logFC>0; strain X vs mean of other 15)
  B  presence/absence only (present >=1 read in >=2/3 reps of X, absent in all other 15 strains)
  C  strain-specific = A & B (the final set, = Fig_strain_specific_DA16)
Finding: DA alone flags ~100k-685k per strain (most are merely strain-ENRICHED, not specific); presence/absence
alone is essentially identical to the intersection -> presence/absence is the specificity-defining filter and DA
is nearly redundant given it. Sources: edger16/da_only_counts.csv (parsed from edger16 logs: DA-only + strain-
specific); edger16/<tp>.presence_only_counts.csv (presence_only_counts.R, filterByExpr + presence, no glmQLFit)."""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U  = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPO = ["E16.5", "P12.5", "P20.5"]
TPC = {"E16.5": "#0072B2", "P12.5": "#009E73", "P20.5": "#D55E00"}
TPLAB = {"E16.5": "E16.5 (prepachytene)", "P12.5": "P12.5", "P20.5": "P20.5 (pachytene)"}
da = pd.read_csv(f"{U}/da_only_counts.csv"); da["tp"] = da.timepoint.map(TPMAP)
po = pd.concat([pd.read_csv(f"{U}/{f}.presence_only_counts.csv") for f in TPMAP]); po["tp"] = po.timepoint.map(TPMAP)
cov = pd.concat([pd.read_csv(f"{U}/{f}.coverage_probe.csv") for f in TPMAP]); cov["tp"] = cov.timepoint.map(TPMAP)   # n_presence_2read2 = ADOPTED ≥2-read presence/absence-only
sp2 = pd.concat([pd.read_csv(f"{U}/{f}.strain_specific_DA_2read.csv.gz") for f in TPMAP])                            # ADOPTED ≥2-read intersection (DA ∩ ≥2-read)
sp2c = sp2.groupby(["strain","timepoint"]).size().reset_index(name="n"); sp2c["tp"] = sp2c.timepoint.map(TPMAP)
CANON = [s for s in STRAIN_ORDER if s in set(da.strain)]
WPOS = [i for i, s in enumerate(CANON) if s in WILD]
def pv(df, val): return df.pivot(index="strain", columns="tp", values=val).reindex(index=CANON, columns=TPO).fillna(0)
mats = {"A  edgeR DA only (FDR<0.05 & logFC>0) — strain-ENRICHED, not specific": pv(da, "da_only"),
        "B  presence/absence only (≥2-read absence — ADOPTED; present ≥2/3 reps in X, absent <2 reads in all 15 others)": pv(cov, "n_presence_2read2"),
        "C  strain-specific = DA ∩ ≥2-read presence/absence (ADOPTED final set = Fig_strain_specific_DA16)": pv(sp2c, "n")}
def klab(v):
    if v >= 1e6: return f"{v/1e6:.1f}M"
    if v >= 1e4: return f"{v/1e3:.0f}k"
    if v >= 1e3: return f"{v/1e3:.1f}k"
    return f"{int(v)}" if v > 0 else ""
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, axes = plt.subplots(3, 1, figsize=(13, 11.2), dpi=300)
x = np.arange(len(CANON)); w = 0.27
for ax, (title, mat) in zip(axes, mats.items()):
    ymax = mat.values.max()
    if WPOS: ax.axvspan(min(WPOS) - 0.5, max(WPOS) + 0.5, color="#C0392B", alpha=0.06, zorder=0)
    for j, t in enumerate(TPO):
        xs = x + (j - 1) * w; vals = mat[t].values.astype(float)
        ax.bar(xs, vals, w, color=TPC[t], edgecolor="white", linewidth=0.3, label=TPLAB[t] if ax is axes[0] else None, zorder=3)
        for xi, v in zip(xs, vals):
            if v > 0: ax.text(xi, v * 1.18, klab(v), ha="center", va="bottom", fontsize=4.2, rotation=90, color=TPC[t])
    ax.set_yscale("log"); ax.set_ylim(1, ymax * 4.5)
    ax.set_xticks(x); labs = ax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7.5)
    for lab, s in zip(labs, CANON): lab.set_color("#C0392B" if s in WILD else "#333")
    ax.set_ylabel("# piRNAs (log)", fontsize=8); ax.set_title(title, fontsize=9.2, fontweight="bold", loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
axes[2].set_xticklabels([])   # bottom-panel strain labels carried by the classical/wild companion below
from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=TPC[t], label=TPLAB[t]) for t in TPO], loc="upper center",
           bbox_to_anchor=(0.5, 0.995), ncol=3, frameon=False, fontsize=8, columnspacing=1.8)
fig.suptitle("Filter decomposition of strain-specific piRNA calls — DA-only vs presence/absence-only vs their intersection (16 strains × 3 timepoints)", fontsize=10.5, fontweight="bold", y=1.035)
fig.text(0.5, -0.02,
    "edgeR DA alone flags ~100k–685k piRNAs per strain (mostly just strain-ENRICHED, not specific) · presence/absence alone ≈ the intersection -> presence/absence is the specificity-defining filter, DA is nearly redundant given it · "
    "x = canonical order, red = wild-derived (they dominate the specific set) · y = log",
    ha="center", fontsize=6, color="#666")
fig.tight_layout(rect=[0, 0, 1, 0.965])
# classical(blue)/wild(orange) companion: total strain-specific (intersection) piRNAs per strain
fig.subplots_adjust(bottom=0.20)
_tot=list(mats.values())[2].sum(1).reindex(CANON).fillna(0).values
_cax=add_classical_wild_companion(fig, axes[2], CANON, _tot, gap=0.10, height_frac=0.20, ylabel="total\n(log)")
_cax.set_xticks(x); _cax.set_xticklabels([s.replace("_","/") for s in CANON], rotation=45, ha="right", fontsize=6.5)
for lab,s in zip(_cax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — total strain-specific piRNAs per strain (intersection)", fontsize=7.5, fontweight="bold", loc="left")
for e in ("pdf", "svg", "png"): fig.savefig(f"{U}/Fig_strain_specific_DA16_decomposition.{e}", bbox_inches="tight")
out = da[["strain", "timepoint", "da_only", "strain_specific"]].merge(po[["strain", "timepoint", "presence_only"]], on=["strain", "timepoint"])
out.to_csv(f"{SD}/Fig_strain_specific_DA16_decomposition.csv", index=False)
print("wrote Fig_strain_specific_DA16_decomposition.{png,pdf,svg} + source data")
print("presence-only vs intersection identical?", (out.presence_only >= out.strain_specific).all(), "| max diff", int((out.presence_only - out.strain_specific).abs().max()))
