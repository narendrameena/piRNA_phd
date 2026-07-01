#!/usr/bin/env python3
"""Developmental timepoint origin for EACH class × strain (the full 3-way breakdown; panels C/D of
Fig_step4_classification16 show only the two 2-way marginals: strain×tp and class×tp). One facet per klass5
class; 16 strains (canonical order) × 3 timepoints, stacked proportion; text = n per strain."""
import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
d = pd.read_csv(f"{U}/unique16/final_classified_clean.csv.gz")
CANON = [s for s in STRAIN_ORDER if s in set(d.strain)]
KL = ["expressed elsewhere (exact)", "SNP-variant (1-3mm)", "low-quality: no mm0 own-genome locus", "unique: conserved-but-silent", "unique: strain-private locus"]
KLAB = {"expressed elsewhere (exact)": "expressed-elsewhere (not unique)", "SNP-variant (1-3mm)": "SNP-variant (allelic, not unique)", "low-quality: no mm0 own-genome locus": "low-quality (no mm0 own locus)", "unique: conserved-but-silent": "unique: conserved-but-silent", "unique: strain-private locus": "unique: strain-private locus"}
TPC = {"16.5dpc": "#0072B2", "12.5dpp": "#009E73", "20.5dpp": "#D55E00"}; TPL = {"16.5dpc": "E16.5 (prepachytene)", "12.5dpp": "P12.5", "20.5dpp": "P20.5 (pachytene)"}
x = np.arange(len(CANON))
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, axes = plt.subplots(len(KL), 1, figsize=(13, 15.5), dpi=300)
src = []
for ax, k in zip(axes, KL):
    sub = d[d.klass5 == k]
    ct = pd.crosstab(sub.strain, sub.timepoint).reindex(index=CANON).fillna(0)
    for t in TPC:
        if t not in ct.columns: ct[t] = 0.0
    n = ct[list(TPC)].sum(1); prop = ct[list(TPC)].div(n.replace(0, 1), axis=0); bot = np.zeros(len(CANON))
    for t in TPC:
        ax.bar(x, prop[t].values, 0.74, bottom=bot, color=TPC[t], label=TPL[t]); bot += prop[t].values
    ax.set_xticks(x); ax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7); ax.set_ylim(0, 1)
    ax.set_ylabel("timepoint\n(proportion)", fontsize=8)
    ax.set_title(f"{KLAB[k]}  —  per-strain timepoint composition (text = n;  class total {int(n.sum()):,})", fontsize=9.4, fontweight="bold", loc="left", pad=10)
    for xi, s in zip(x, CANON): ax.text(xi, 1.012, f"{int(n[s]):,}", ha="center", va="bottom", fontsize=4.6, color="#777")
    ax.spines[["top", "right"]].set_visible(False)
    src.append(ct[list(TPC)].assign(klass=KLAB[k]))
axes[0].legend(fontsize=7.5, frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.18))
fig.suptitle("Developmental timepoint origin for EACH class × strain — 16-strain pangenome unique-piRNA classification (3-way breakdown)", fontsize=11, fontweight="bold", y=0.999)
pd.concat(src).to_csv(f"{PG}/SourceData_step4_class_strain_timepoint16.csv")
fig.tight_layout(rect=[0, 0, 1, 0.99])
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_step4_class_strain_timepoint16.{e}", bbox_inches="tight")
print("wrote Fig_step4_class_strain_timepoint16.{png,pdf,svg} + source data")
