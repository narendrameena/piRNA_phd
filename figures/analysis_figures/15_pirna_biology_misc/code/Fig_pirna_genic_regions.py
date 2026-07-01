#!/usr/bin/env python3
"""piRNA genomic-REGION overlap, in parallel across 16 strains x 3 timepoints (the thesis genic-region analysis,
restricted to UNIQUE / strain-private piRNAs). Each sample's unique-piRNA signal partitions into 3 mutually-exclusive
region classes from featureCounts: protein-coding GENE / lncRNA / INTERGENIC. Source = thesis output
analysis/sRNA_deseq/genric_regions/uniq_piRNA_list1_count_v3.3.csv (already summed per sample x class). Replicates are
summed per strain x timepoint, then converted to a fraction. Layout = three timepoint panels (E16.5/P12.5/P20.5),
strains side-by-side in the canonical thesis order (wild strains in red). The biological question this answers:
do wild strains (CAST/SPRET/PWK), which carry the large open pan-piRNA-ome, route MORE of their strain-private
piRNAs through intergenic (candidate novel-locus) space vs genic/lncRNA? Purely descriptive here; interpretation
is queued for BioMNI verification."""
import sys, re, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
SRC = f"{ROOT}/analysis/sRNA_deseq/genric_regions/uniq_piRNA_list1_count_v3.3.csv"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
SD = f"{ROOT}/analysis/claude_biomni_analysis/source_data"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]          # C57BL_6 = external public data, absent here
TPS = ["E16.5", "P12.5", "P20.5"]
CLASSES = [("gene_gene", "protein-coding gene", "#4575B4"),  # CB-safe: blue / orange / purple
           ("lnc_RNA_transcript", "lncRNA", "#FDAE61"),
           ("intergenic", "intergenic", "#7B3294")]
df = pd.read_csv(SRC)                                        # cols: V8 (sample), V9 (class), summed_V7 (count)
def parse(v):
    strain = v.split("/")[0]
    m = re.search(r"(E16\.5|P12\.5|P20\.5)", v)
    return strain, (m.group(1) if m else None)
df["strain"], df["tp"] = zip(*df.V8.map(parse))
g = df.groupby(["strain", "tp", "V9"])["summed_V7"].sum().reset_index()        # sum replicates
piv = g.pivot_table(index=["strain", "tp"], columns="V9", values="summed_V7", fill_value=0.0)
for cls, _, _ in CLASSES:
    if cls not in piv.columns: piv[cls] = 0.0
frac = piv[[c for c, _, _ in CLASSES]].div(piv[[c for c, _, _ in CLASSES]].sum(axis=1), axis=0)
# ---- source data ----
os.makedirs(SD, exist_ok=True)
out = frac.copy(); out.columns = [n for _, n, _ in CLASSES]
out.reset_index().to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/data/source_data/SourceData_pirna_genic_regions.csv", index=False)
inter = out["intergenic"]
print("intergenic fraction by strain (mean over tp):")
for X in CANON:
    vals = [out.loc[(X, tp), "intergenic"] for tp in TPS if (X, tp) in out.index]
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
    ax.axvspan(len(CANON) - 3.5, len(CANON) - 0.5, color="#C0392B", alpha=0.05, zorder=0)  # wild block
axes[0].set_ylabel("fraction of unique-piRNA signal", fontsize=10)
leg = [Patch(facecolor=col, label=name) for _, name, col in CLASSES]
fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=10, frameon=False,
           title="genomic-region class of strain-private piRNAs (mutually exclusive)", title_fontsize=9.5)
fig.suptitle("Genomic-region overlap of UNIQUE (strain-private) piRNAs — 16 strains × 3 timepoints, in parallel\n"
             "each sample's strain-private piRNA signal split into protein-coding gene / lncRNA / intergenic (thesis genic-region method); wild strains in red",
             fontsize=12, fontweight="bold", y=1.02, linespacing=1.5)
fig.tight_layout(rect=[0, 0.04, 1, 0.98])
for e in ("pdf", "svg", "png"):
    fig.savefig(f"{U}/pangenome_te/Fig_pirna_genic_regions.{e}", bbox_inches="tight")
print("wrote Fig_pirna_genic_regions.{png,pdf,svg} + SourceData_pirna_genic_regions.csv")
