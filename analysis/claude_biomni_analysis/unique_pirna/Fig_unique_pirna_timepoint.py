#!/usr/bin/env python3
"""Which timepoint produces each strain's UNIQUE (strain-private-locus) piRNAs? From unique16/final_classified.csv.gz
(klass == 'unique: strain-private locus'), count strain-private piRNA detections per (strain, timepoint). Two panels,
strains in canonical order (wild in red): (A) absolute counts, stacked by timepoint (log y because wild strains
dominate); (B) timepoint COMPOSITION (% per strain). Biology: prepachytene (E16.5/16.5dpc) vs pachytene
(P20.5/20.5dpp) windows differ in piRNA programme, so the timepoint mix reveals when a strain's private piRNAs are
born. Descriptive; framing queued for BioMNI."""
import sys, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
SD = f"{ROOT}/analysis/claude_biomni_analysis/source_data"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}
TPS = [("E16.5", "#4393C3"), ("P12.5", "#FDB863"), ("P20.5", "#B2182B")]   # developmental blue->orange->red
d = pd.read_csv(f"{U}/unique16/final_classified_clean.csv.gz")   # mm0-clean strain-private (klass5); excludes mm1-3/error 'low-quality' reads
g = d[d.klass5 == "unique: strain-private locus"].copy()
g["tp"] = g.timepoint.map(TPMAP)
ct = g.groupby(["strain", "tp"]).size().unstack(fill_value=0).reindex(index=CANON, columns=[t for t, _ in TPS]).fillna(0)
os.makedirs(SD, exist_ok=True); ct.to_csv(f"{SD}/SourceData_unique_pirna_timepoint.csv")
frac = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0) * 100
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, (axA, axB, axC) = plt.subplots(3, 1, figsize=(11, 11.0), dpi=300, gridspec_kw={"hspace": 0.62})
x = np.arange(len(CANON))
bottom = np.zeros(len(CANON))
for tp, col in TPS:
    axA.bar(x, ct[tp].values, bottom=bottom, width=0.8, color=col, edgecolor="white", linewidth=0.3, label=tp)
    bottom = bottom + ct[tp].values
axA.set_yscale("log"); axA.set_ylim(1, max(ct.sum(axis=1).max() * 1.6, 10))
axA.set_xticks(x); axA.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7.5)
for t, X in zip(axA.get_xticklabels(), CANON):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.set_ylabel("# strain-private piRNAs\n(stacked, log)", fontsize=8.5)
axA.set_title("A  Strain-private piRNA count by timepoint (wild strains dominate — open pan-piRNA-ome)", fontsize=9.4, fontweight="bold", loc="left")
axA.legend(handles=[Patch(facecolor=c, label=t) for t, c in TPS], ncol=3, fontsize=8, frameon=False, loc="upper left")
bottom = np.zeros(len(CANON))
for tp, col in TPS:
    axB.bar(x, frac[tp].values, bottom=bottom, width=0.8, color=col, edgecolor="white", linewidth=0.3)
    bottom = bottom + np.nan_to_num(frac[tp].values)
axB.set_ylim(0, 100); axB.set_xticks(x); axB.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7.5)
for t, X in zip(axB.get_xticklabels(), CANON):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axB.set_ylabel("timepoint composition (%)", fontsize=8.5)
axB.set_title("B  Timepoint composition per strain (which developmental window is responsible)", fontsize=9.4, fontweight="bold", loc="left")
# Panel C — EXACT (linear) counts, total labelled
bottom = np.zeros(len(CANON))
for tp, col in TPS:
    axC.bar(x, ct[tp].values, bottom=bottom, width=0.8, color=col, edgecolor="white", linewidth=0.3)
    bottom = bottom + ct[tp].values
tot = ct.sum(axis=1).values
for xi, t in zip(x, tot):
    if t > 0: axC.text(xi, t + tot.max() * 0.012, f"{int(t):,}", ha="center", va="bottom", fontsize=6.2)
axC.set_ylim(0, tot.max() * 1.12)
axC.set_xticks(x); axC.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7.5)
for t, X in zip(axC.get_xticklabels(), CANON):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axC.set_ylabel("# strain-private piRNAs\n(EXACT, linear)", fontsize=8.5)
axC.set_title("C  Strain-private piRNA count by timepoint — EXACT counts (linear scale; total labelled per strain)", fontsize=9.4, fontweight="bold", loc="left")
fig.suptitle("Which timepoint produces each strain's UNIQUE (strain-private-locus) piRNAs — 16 strains × 3 developmental windows",
             fontsize=11, fontweight="bold", y=0.98)
for e in ("pdf", "svg", "png"): fig.savefig(f"{U}/pangenome_te/Fig_unique_pirna_timepoint.{e}", bbox_inches="tight")
print("wrote Fig_unique_pirna_timepoint + source data")
print(ct.astype(int).to_string())
