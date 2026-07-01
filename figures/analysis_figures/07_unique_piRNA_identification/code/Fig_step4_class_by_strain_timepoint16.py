#!/usr/bin/env python3
"""Unique-piRNA classification by strain AND timepoint — TWO panels, each strains × timepoints.
A (NUMBERS): total candidates per strain at each timepoint (grouped bars, log; colour = timepoint) — the magnitude,
which spans 0-95k. B (PERCENTAGE): 5-class composition per strain at each timepoint (grouped STACKED bars; within
each strain group the three sub-bars are E16.5 / P12.5 / P20.5 left->right). Wild strains red. klass5 ≥2-read."""
import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")   # ADOPTED ≥2-read absence
CANON = [s for s in STRAIN_ORDER if s in set(d.strain)]; WPOS = [i for i, s in enumerate(CANON) if s in WILD]
KL  = ["expressed elsewhere (exact)", "SNP-variant (1-3mm)", "low-quality: no mm0 own-genome locus", "unique: conserved-but-silent", "unique: strain-private locus"]
LAB = ["expressed-elsewhere (not unique)", "SNP-variant (allelic — not unique)", "low-quality (no mm0 own locus)", "unique: conserved-but-silent", "unique: strain-private locus (clean)"]
COL = ["#9e9e9e", "#E69F00", "#cdb892", "#0072B2", "#7a3b9a"]; WHITE_SEG = {"#9e9e9e", "#E69F00", "#0072B2", "#7a3b9a"}
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPO = ["E16.5", "P12.5", "P20.5"]; TPCOL = {"E16.5": "#4393C3", "P12.5": "#E8852B", "P20.5": "#B2182B"}
d["tp"] = d.timepoint.map(TPMAP)
x = np.arange(len(CANON)); nb = len(TPO); bw = 0.82 / nb
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(15, 12.6), dpi=300)
gs = fig.add_gridspec(3, 1, height_ratios=[1, 0.11, 1], hspace=0.30)
axN = fig.add_subplot(gs[0]); axL = fig.add_subplot(gs[1]); axP = fig.add_subplot(gs[2], sharex=axN); axL.axis("off")   # axL = dedicated legend band between the panels
# ---- Panel A: NUMBERS (total candidates per strain x timepoint; grouped bars, log) ----
for sp_ax in (axN, axP):
    if WPOS: sp_ax.axvspan(min(WPOS) - 0.5, max(WPOS) + 0.5, color="#C0392B", alpha=0.06, zorder=0)
tot = {}
for j, tp in enumerate(TPO):
    vals = [int(((d.strain == s) & (d.tp == tp)).sum()) for s in CANON]; tot[tp] = vals
    axN.bar(x + (j - 1) * bw, vals, bw, color=TPCOL[tp], edgecolor="white", linewidth=0.3, label=tp, zorder=3)
    for xi, v in zip(x + (j - 1) * bw, vals):
        if v > 0: axN.text(xi, v * 1.18, f"{v//1000}k" if v >= 1000 else f"{v}", ha="center", va="bottom", fontsize=4.6, rotation=90, color=TPCOL[tp], fontweight="bold")
        else: axN.text(xi, 1.2, "0", ha="center", va="bottom", fontsize=6, color=TPCOL[tp], fontweight="bold")   # explicit 0 (no bar on log scale) — e.g. FVB/NJ & NOD/ShiLtJ at P20.5: libraries exist but 0 candidates survived DA+absence
axN.set_yscale("log"); axN.set_ylim(1, max(max(v) for v in tot.values()) * 3); axN.set_ylabel("total candidates (log)", fontsize=9.5)
axN.legend(fontsize=8.5, frameon=False, ncol=3, loc="upper left", title="timepoint", title_fontsize=8.5)
axN.set_title("A  NUMBERS — total strain-specific candidates per strain at each timepoint", fontsize=10.5, fontweight="bold", loc="left")
axN.tick_params(labelbottom=False)   # strain labels only once (rotated) on Panel B — avoids the overlapping horizontal row here
axN.spines[["top", "right"]].set_visible(False)
# ---- Panel B: PERCENTAGE (5-class composition per strain x timepoint; grouped stacked bars) ----
for j, tp in enumerate(TPO):
    sub = d[d.tp == tp]; ct = pd.crosstab(sub.strain, sub.klass5).reindex(index=CANON, columns=KL).fillna(0)
    prop = ct.div(ct.sum(1).replace(0, 1), axis=0); bot = np.zeros(len(CANON))
    for i, k in enumerate(KL):
        axP.bar(x + (j - 1) * bw, prop[k].values, bw, bottom=bot, color=COL[i], edgecolor="white", linewidth=0.2, label=LAB[i] if j == 0 else None, zorder=3)
        for xi, v, b in zip(x + (j - 1) * bw, prop[k].values, bot):
            if v >= 0.10: axP.text(xi, b + v / 2, f"{v*100:.0f}", ha="center", va="center", fontsize=3.9, fontweight="bold", color="white" if COL[i] in WHITE_SEG else "#333", zorder=4)
        bot += prop[k].values
axP.set_ylim(0, 1); axP.set_ylabel("class composition (proportion)", fontsize=9.5)
axP.set_title("B  PERCENTAGE — 5-class composition per strain at each timepoint  (within each strain: E16.5 · P12.5 · P20.5, left→right)", fontsize=10.5, fontweight="bold", loc="left")
axL.legend(handles=[Patch(facecolor=COL[i], label=LAB[i]) for i in range(len(KL))], fontsize=8.2, frameon=False, ncol=5, loc="center", columnspacing=1.6, handlelength=1.5, title="5-class composition  (Panel B segments)", title_fontsize=8.8)
axP.set_xticks(x); axP.set_xticklabels([])   # strain labels carried by the classical/wild companion below
axP.spines[["top", "right"]].set_visible(False)
fig.suptitle("Unique-piRNA classification by strain AND timepoint — NUMBERS (A) and PERCENTAGE composition (B), 16 strains × 3 windows", fontsize=11.5, fontweight="bold", y=0.995)
pd.DataFrame({tp: tot[tp] for tp in TPO}, index=CANON).to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data/SourceData_step4_class_by_strain_timepoint16_numbers.csv")
fig.text(0.5, 0.006, "“0” above a bar = that strain × timepoint had 0 strain-specific candidates surviving the edgeR DA + ≥2-read absence filter (the sRNA libraries DO exist; e.g. FVB/NJ and NOD/ShiLtJ at P20.5 — the signal collapses by pachytene).", ha="center", fontsize=7.5, color="#666")
fig.tight_layout(rect=[0, 0.02, 1, 0.97])
# classical(blue)/wild(orange) companion: total candidates per strain (subspecies colour scheme)
fig.subplots_adjust(bottom=0.16)
_tot=d.strain.value_counts().reindex(CANON).fillna(0).values
_cax=add_classical_wild_companion(fig, axP, CANON, _tot, gap=0.085, height_frac=0.16, ylabel="total\n(log)")
_cax.set_xticks(x); _cax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=6.5)
for lab, s in zip(_cax.get_xticklabels(), CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — total candidates per strain", fontsize=7.5, fontweight="bold", loc="left")
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_step4_class_by_strain_timepoint16.{e}", bbox_inches="tight")
print("wrote Fig_step4_class_by_strain_timepoint16 (2-panel: numbers + percentage)")
