#!/usr/bin/env python3
"""piRNA genic-FEATURE overlap for ALL piRNAs (companion to the unique-piRNA coarse figure), in parallel across
16 strains x 3 timepoints. Source = thesis output analysis/sRNA_deseq/genric_regions/list2_count.csv (all piRNA-seq
expression summed per sample x genic feature). list2 reports CDS / exon / 5'UTR / 3'UTR / intron, but these are NOT
all mutually exclusive (exon = CDS u 5'UTR u 3'UTR u non-coding-exon). To keep the stacked bar a true partition we
use the mutually-exclusive set {CDS, 5'UTR, 3'UTR, intron} and DROP the redundant 'exon' superset; fractions are
within the gene body. Replicates summed per strain x timepoint. Three timepoint panels, strains in canonical order
(wild in red). This answers: within genes, do all piRNAs sit in coding vs UTR vs intronic sequence, and does that
shift by strain/timepoint? Descriptive; interpretation queued for BioMNI."""
import sys, re, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
SRC = f"{ROOT}/analysis/sRNA_deseq/genric_regions/list2_count.csv"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
SD = f"{ROOT}/analysis/claude_biomni_analysis/source_data"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
TPS = ["E16.5", "P12.5", "P20.5"]
# mutually-exclusive gene-body partition (exon dropped as the CDS+UTR superset); CB-safe sequential
CLASSES = [("CDS", "CDS (coding)", "#1A9850"),
           ("five_prime_UTR", "5′UTR", "#91CF60"),
           ("three_prime_UTR", "3′UTR", "#FEE08B"),
           ("intron", "intron", "#998EC3")]
df = pd.read_csv(SRC)                                   # cols: V8 (sample), V9 (feature), summed_V7 (count)
def parse(v):
    strain = v.split("/")[0]
    m = re.search(r"(16\.5dpc|12\.5dpp|20\.5dpp)", v)
    tp = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}.get(m.group(1)) if m else None
    return strain, tp
df["strain"], df["tp"] = zip(*df.V8.map(parse))
g = df.groupby(["strain", "tp", "V9"])["summed_V7"].sum().reset_index()
piv = g.pivot_table(index=["strain", "tp"], columns="V9", values="summed_V7", fill_value=0.0)
for cls, _, _ in CLASSES:
    if cls not in piv.columns: piv[cls] = 0.0
sub = piv[[c for c, _, _ in CLASSES]]
frac = sub.div(sub.sum(axis=1), axis=0)
# ---- source data ----
os.makedirs(SD, exist_ok=True)
out = frac.copy(); out.columns = [n for _, n, _ in CLASSES]
out.reset_index().to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/data/source_data/SourceData_pirna_genic_features_all.csv", index=False)
print("intron fraction by strain (mean over tp):")
for X in CANON:
    vals = [out.loc[(X, tp), "intron"] for tp in TPS if (X, tp) in out.index]
    if vals: print(f"  {'*' if X in WILD else ' '} {X:14s} {np.mean(vals):.3f}")
# ---- figure ----
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), dpi=300, sharey=True)
x = np.arange(len(CANON))
for ax, tp in zip(axes, TPS):
    bottom = np.zeros(len(CANON))
    for cls, name, col in CLASSES:
        h = np.array([frac.loc[(X, tp), cls] if (X, tp) in frac.index else np.nan for X in CANON])
        ax.bar(x, h, bottom=bottom, width=0.82, color=col, edgecolor="white", linewidth=0.3)
        bottom = bottom + np.nan_to_num(h)
    ax.set_title(tp, fontsize=11, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=90, fontsize=7)
    for t, X in zip(ax.get_xticklabels(), CANON):
        if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
    ax.set_ylim(0, 1); ax.set_xlim(-0.7, len(CANON) - 0.3)
    ax.axvspan(len(CANON) - 3.5, len(CANON) - 0.5, color="#C0392B", alpha=0.05, zorder=0)
axes[0].set_ylabel("fraction of all-piRNA signal within gene body", fontsize=10)
leg = [Patch(facecolor=col, label=name) for _, name, col in CLASSES]
fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=4, fontsize=10, frameon=False,
           title="gene-body feature class of all piRNAs (mutually exclusive; 'exon' superset excluded)", title_fontsize=9.5)
fig.suptitle("Genic-FEATURE overlap of ALL piRNAs — 16 strains × 3 timepoints, in parallel\n"
             "all piRNA-seq signal within genes split into CDS / 5′UTR / 3′UTR / intron (thesis genic-region method, list2); wild strains in red",
             fontsize=12, fontweight="bold", y=1.02, linespacing=1.5)
fig.tight_layout(rect=[0, 0.04, 1, 0.98])
for e in ("pdf", "svg", "png"):
    fig.savefig(f"{U}/pangenome_te/Fig_pirna_genic_features_all.{e}", bbox_inches="tight")
print("wrote Fig_pirna_genic_features_all.{png,pdf,svg} + SourceData_pirna_genic_features_all.csv")
