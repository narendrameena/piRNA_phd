#!/usr/bin/env python3
"""piRNA genomic-region overlap for ALL thesis region-lists, in parallel across 16 strains x 3 timepoints. The
thesis genic-region analysis (analysis/sRNA_deseq/genric_regions) produced several featureCounts summaries that
differ in annotation granularity; this renders each as a 16-strain x 3-timepoint stacked-bar figure. Category =
the leading token of V9 (before ':'), so it works for both the coarse list (gene/lncRNA/intergenic) and the
gene-body lists (CDS/exon/UTR/intron). For the CDS-less gene-body lists (list3/list4) we derive a mutually-
exclusive partition exon_other = max(0, exon - 5'UTR - 3'UTR) so the stack does not double-count (exon is the
UTR/CDS superset). list1-uniq and list2-all already have dedicated figures (Fig_pirna_genic_regions /
Fig_pirna_genic_features_all); this adds the missing companions. Replicates summed per strain x timepoint, then
fraction. Source data saved per figure. Descriptive; biology queued for BioMNI."""
import sys, re, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
GR = f"{ROOT}/analysis/sRNA_deseq/genric_regions"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
SD = f"{ROOT}/analysis/claude_biomni_analysis/source_data"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
TPS = ["E16.5", "P12.5", "P20.5"]
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}
GENIC = [("CDS_or_exon", "CDS / coding exon", "#1A9850"), ("five_prime_UTR", "5′UTR", "#91CF60"),
         ("three_prime_UTR", "3′UTR", "#FEE08B"), ("intron", "intron", "#998EC3")]
# (out_name, file, scope-label, ylab, cats, derive_exon_other)
LISTS = [
    ("Fig_pirna_genic_features_uniq", "uniq_piRNA_list2_count.csv", "UNIQUE piRNAs — gene-body feature (list2)",
     "fraction of unique-piRNA signal within gene body",
     [("CDS", "CDS (coding)", "#1A9850")] + GENIC[1:], False),
    ("Fig_pirna_genic_list3_uniq", "uniq_piRNA_list3_count.csv", "UNIQUE piRNAs — gene-body feature (list3, CDS-less)",
     "fraction of unique-piRNA signal within gene body", GENIC, True),
    ("Fig_pirna_genic_list3_all", "list3_count.csv", "ALL piRNAs — gene-body feature (list3, CDS-less)",
     "fraction of all-piRNA signal within gene body", GENIC, True),
    ("Fig_pirna_genic_list4_uniq", "uniq_piRNA_list4_count.csv",
     "UNIQUE piRNAs — transcript-level gene-body feature (list4, per-strain Ensembl annotation)",
     "fraction of unique-piRNA signal within gene body", GENIC, True),
]
def render(out, fn, scope, ylab, cats, derive):
    df = pd.read_csv(f"{GR}/{fn}")
    df["strain"] = df.V8.str.split("/").str[0]
    df["tp"] = df.V8.str.extract(r"(16\.5dpc|12\.5dpp|20\.5dpp)")[0].map(TPMAP)
    df["cat"] = df.V9.str.split(":").str[0]
    piv = df.groupby(["strain", "tp", "cat"])["summed_V7"].sum().reset_index() \
            .pivot_table(index=["strain", "tp"], columns="cat", values="summed_V7", fill_value=0.0)
    if derive:                                  # exon_other = exon - UTRs (clamp >=0) -> mutually exclusive
        for c in ("exon", "five_prime_UTR", "three_prime_UTR"):
            if c not in piv.columns: piv[c] = 0.0
        piv["CDS_or_exon"] = (piv["exon"] - piv["five_prime_UTR"] - piv["three_prime_UTR"]).clip(lower=0)
    for c, _, _ in cats:
        if c not in piv.columns: piv[c] = 0.0
    sub = piv[[c for c, _, _ in cats]]; frac = sub.div(sub.sum(axis=1), axis=0)
    os.makedirs(SD, exist_ok=True)
    o = frac.copy(); o.columns = [n for _, n, _ in cats]
    o.reset_index().to_csv(f"{SD}/SourceData_{out}.csv", index=False)
    plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), dpi=300, sharey=True)
    x = np.arange(len(CANON))
    for ax, tp in zip(axes, TPS):
        bottom = np.zeros(len(CANON))
        for c, name, col in cats:
            h = np.array([frac.loc[(X, tp), c] if (X, tp) in frac.index else np.nan for X in CANON])
            ax.bar(x, h, bottom=bottom, width=0.82, color=col, edgecolor="white", linewidth=0.3)
            bottom = bottom + np.nan_to_num(h)
        ax.set_title(tp, fontsize=11, fontweight="bold")
        ax.set_xticks(x); ax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=90, fontsize=7)
        for t, X in zip(ax.get_xticklabels(), CANON):
            if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
        ax.set_ylim(0, 1); ax.set_xlim(-0.7, len(CANON) - 0.3)
        ax.axvspan(len(CANON) - 3.5, len(CANON) - 0.5, color="#C0392B", alpha=0.05, zorder=0)
    axes[0].set_ylabel(ylab, fontsize=10)
    leg = [Patch(facecolor=col, label=name) for _, name, col in cats]
    fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=len(cats), fontsize=10,
               frameon=False, title="region class (mutually exclusive)", title_fontsize=9.5)
    fig.suptitle(f"Genomic-region overlap — {scope}, 16 strains × 3 timepoints, in parallel\n"
                 "thesis genic-region method; wild strains in red", fontsize=12, fontweight="bold", y=1.02, linespacing=1.5)
    fig.tight_layout(rect=[0, 0.04, 1, 0.98])
    for e in ("pdf", "svg", "png"): fig.savefig(f"{U}/pangenome_te/{out}.{e}", bbox_inches="tight")
    print(f"wrote {out} + SourceData_{out}.csv")
if __name__ == "__main__":
    for cfg in LISTS: render(*cfg)
