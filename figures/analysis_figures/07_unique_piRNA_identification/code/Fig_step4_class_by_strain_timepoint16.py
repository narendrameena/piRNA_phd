#!/usr/bin/env python3
"""Unique-piRNA CLASS composition per strain, FACETED BY TIMEPOINT — the classification annotated by strain
AND timepoint together. 3 panels (E16.5 / P12.5 / P20.5); each = 16 strains (canonical order; wild = red labels
+ shaded block) × stacked 5-class composition; % above each bar = genuinely-unique fraction; n above = total
candidates that strain×timepoint. Complements Fig_step4_classification16 (Panel B is pooled over timepoints)
and Fig_step4_class_strain_timepoint16 (which stacks timepoint and facets by class — the other marginal)."""
import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")   # ADOPTED ≥2-read absence
CANON = [s for s in STRAIN_ORDER if s in set(d.strain)]
WPOS = [i for i, s in enumerate(CANON) if s in WILD]
KL  = ["expressed elsewhere (exact)", "SNP-variant (1-3mm)", "low-quality: no mm0 own-genome locus", "unique: conserved-but-silent", "unique: strain-private locus"]
LAB = ["expressed-elsewhere (not unique)", "SNP-variant (allelic — not unique)", "low-quality (no mm0 own locus)", "unique: conserved-but-silent", "unique: strain-private locus (clean)"]
COL = ["#9e9e9e", "#E69F00", "#cdb892", "#0072B2", "#7a3b9a"]   # same 5-class palette as Fig_step4_classification16
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPO = ["E16.5", "P12.5", "P20.5"]
d["tp"] = d.timepoint.map(TPMAP)
x = np.arange(len(CANON))
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, axes = plt.subplots(3, 1, figsize=(13, 13.2), dpi=300)
src = []
for ax, tp in zip(axes, TPO):
    sub = d[d.tp == tp]
    ct = pd.crosstab(sub.strain, sub.klass5).reindex(index=CANON, columns=KL).fillna(0)
    n = ct.sum(1); prop = ct.div(n.replace(0, 1), axis=0); bot = np.zeros(len(CANON))
    if WPOS: ax.axvspan(min(WPOS) - 0.5, max(WPOS) + 0.5, color="#C0392B", alpha=0.06, zorder=0)   # wild-derived block
    for i, (k, lab) in enumerate(zip(KL, LAB)):
        ax.bar(x, prop[k].values, 0.74, bottom=bot, color=COL[i], label=lab if ax is axes[0] else None, zorder=3); bot += prop[k].values
    gu = ct[["unique: conserved-but-silent", "unique: strain-private locus"]].sum(1) / n.replace(0, 1) * 100
    ax.set_xticks(x); labs = ax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7.5)
    for lab, s in zip(labs, CANON): lab.set_color("#C0392B" if s in WILD else "#333")
    ax.set_ylim(0, 1); ax.set_ylabel("class composition\n(proportion)", fontsize=8)
    ax.set_title(f"{tp}  —  unique-piRNA class composition per strain  (n above = total candidates;  % = genuinely unique by expression)", fontsize=9.4, fontweight="bold", loc="left", pad=12)
    for xi, s in zip(x, CANON):
        ax.text(xi, 1.012, f"{gu[s]:.0f}%", ha="center", va="bottom", fontsize=5.0, color="#333")
        ax.text(xi, 1.07, f"{int(n[s]):,}", ha="center", va="bottom", fontsize=4.2, color="#999")
    ax.spines[["top", "right"]].set_visible(False)
    src.append(ct.assign(tp=tp))
axes[0].legend(fontsize=7, frameon=False, ncol=5, loc="lower left", bbox_to_anchor=(0, 1.22), columnspacing=1.0, handlelength=1.2)
fig.suptitle("Unique-piRNA classification by strain AND timepoint — 5-class pangenome-syntenic composition per strain at each developmental window", fontsize=11, fontweight="bold", y=0.999)
pd.concat(src).to_csv(f"{PG}/SourceData_step4_class_by_strain_timepoint16.csv")
fig.tight_layout(rect=[0, 0, 1, 0.985])
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_step4_class_by_strain_timepoint16.{e}", bbox_inches="tight")
print("wrote Fig_step4_class_by_strain_timepoint16.{png,pdf,svg} + source data")
print(d.groupby("tp").size())
