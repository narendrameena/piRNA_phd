#!/usr/bin/env python3
"""THEME 07 (additional, non-destructive) -- STAGE-PEAK-LENGTH unique-piRNA subset (user-requested 2026-06-23).

From the genuinely-unique piRNAs (conserved-but-silent + strain-private), keep ONLY the user-specified
stage-characteristic EXACT length(s) per developmental timepoint:
   E16.5 (16.5dpc) -> 27 nt          (pre-pachytene peak)
   P12.5 (12.5dpp) -> 27 nt AND 30 nt (transition: both populations)
   P20.5 (20.5dpp) -> 30 nt          (pachytene peak)
The FULL unique set is untouched; this is a separate subset + figure.

BIOLOGICAL CAVEAT (documented, BioMNI triple-verified 2026-06-23): piRNA length is set by the PIWI footprint
(MILI/MIWI2 ~24-28 nt pre-pachytene; MIWI ~29-32 nt pachytene), each spanning ~3 nt -- so a length WINDOW
(26-28 / 29-31) is biologically more accurate than single exact lengths, and 7 strains actually peak at 28 nt
(MIWI2). The user chose EXACT 27/30 (modal length) with this caveat in front of them: this is the cleanest
single-length stage-characteristic set, not the full pre-pachytene/pachytene populations. Panel A shows the
full distribution so the exact cut and its caveat are transparent.

Inputs : unique16/final_classified_clean_2read.csv.gz   (full genuinely-unique, for Panel A distribution)
         unique16/unique_stage_peak_27_30.csv.gz         (the exact-27/30 subset, make_stage_peak_unique.py)
Outputs: figures/Fig_unique_stage_peak_lengths.{pdf,svg,png} + data/source_data/SourceData_Fig_unique_stage_peak_lengths*.csv
"""
import sys; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion, color_wild_labels
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
T = f"{ROOT}/figures/analysis_figures/07_unique_piRNA_identification"
FULL = sys.argv[1] if len(sys.argv) > 1 else f"{U}/unique16/final_classified_clean_2read.csv.gz"
SUB  = sys.argv[2] if len(sys.argv) > 2 else f"{U}/unique16/unique_stage_peak_27_30.csv.gz"
GU = ["unique: conserved-but-silent", "unique: strain-private locus"]
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]; TPN = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}
KEEP = {"16.5dpc": {27}, "12.5dpp": {27, 30}, "20.5dpp": {30}}
TPCOL = {"16.5dpc": "#4393C3", "12.5dpp": "#E8852B", "20.5dpp": "#B2182B"}
CCOL = {"unique: conserved-but-silent": "#0072B2", "unique: strain-private locus": "#7a3b9a"}
CLAB = {"unique: conserved-but-silent": "conserved-but-silent", "unique: strain-private locus": "strain-private locus"}

full = pd.read_csv(FULL, usecols=["sequence", "strain", "timepoint", "klass5"])
full = full[full.klass5.isin(GU)].copy(); full["L"] = full.sequence.str.len()
sub = pd.read_csv(SUB)
order = [s for s in STRAIN_ORDER if s in set(sub.strain)]; WPOS = [i for i, s in enumerate(order) if s in WILD]

fig = plt.figure(figsize=(15, 13), dpi=300)
gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.25, 0.95], hspace=0.62, wspace=0.28,
                      left=0.07, right=0.985, top=0.91, bottom=0.08)
# ---- Panel A: length distribution per tp, KEPT lengths highlighted (rationale + breadth caveat) ----
LENS = list(range(25, 33))
for j, tp in enumerate(TPS):
    ax = fig.add_subplot(gs[0, j])
    cnt = full[full.timepoint == tp].L.value_counts().reindex(LENS, fill_value=0)
    cols = [TPCOL[tp] if L in KEEP[tp] else "#d9d9d9" for L in LENS]
    ax.bar(LENS, cnt.values, color=cols, edgecolor=["#222" if L in KEEP[tp] else "none" for L in LENS], linewidth=1.0, zorder=2)
    for L in KEEP[tp]:
        ax.text(L, cnt[L] + cnt.max() * 0.02, f"{int(cnt[L]):,}", ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=TPCOL[tp])
    ax.set_title(f"{TPN[tp]}  (keep {'/'.join(str(x) for x in sorted(KEEP[tp]))} nt)", fontsize=10, fontweight="bold")
    ax.set_xlabel("piRNA length (nt)", fontsize=9); ax.set_xticks(LENS)
    if j == 0: ax.set_ylabel("genuinely-unique piRNAs", fontsize=9.5)
    ax.tick_params(labelsize=8); ax.spines[["top", "right"]].set_visible(False)
fig.text(0.07, 0.925, "A", fontsize=14, fontweight="bold")
fig.text(0.5, 0.925, "Length distribution of genuinely-unique piRNAs per stage — coloured bars = kept (exact peak); grey = excluded",
         ha="center", fontsize=10.5, fontweight="bold")
# ---- Panel B: per-strain stage-peak unique counts, stacked by tp (canonical order, wild band, companion) ----
axB = fig.add_subplot(gs[1, :]); x = np.arange(len(order)); bw = 0.26
piv = sub.groupby(["strain", "timepoint"]).size().unstack(fill_value=0).reindex(index=order, columns=TPS, fill_value=0)
if WPOS: axB.axvspan(min(WPOS) - 0.5, max(WPOS) + 0.5, color="#C0392B", alpha=0.06, zorder=0)
for k, tp in enumerate(TPS):
    axB.bar(x + (k - 1) * bw, piv[tp].values, bw, color=TPCOL[tp], edgecolor="white", linewidth=0.3, label=TPN[tp], zorder=2)
axB.set_yscale("log"); axB.set_ylim(bottom=1)   # 125x range across strains -> log so classical strains stay readable
axB.set_xticks(x); axB.set_xticklabels([s.replace("_", "/") for s in order], rotation=45, ha="right", fontsize=8)
color_wild_labels(axB, order); axB.set_ylabel("stage-peak unique piRNAs\n(27/30 nt subset, log)", fontsize=9.5)
axB.legend(title="timepoint", fontsize=8.5, title_fontsize=8.5, frameon=False, loc="upper left", ncol=3)
axB.set_title("B  Stage-peak (exact 27/30 nt) unique piRNAs per strain, by stage (log — wild-derived strains carry ~10–100× more)", fontsize=10.5, fontweight="bold", loc="left")
axB.spines[["top", "right"]].set_visible(False)
if WPOS: axB.text(np.mean(WPOS), axB.get_ylim()[1] * 0.45, "wild-derived", ha="center", fontsize=8.5, fontweight="bold", color="#C0392B")
add_classical_wild_companion(fig, axB, order, piv.sum(axis=1).values, ylabel="total\n(log)")
# ---- Panel C: class composition (CBS vs strain-private) within the subset, per tp ----
axC = fig.add_subplot(gs[2, 0])
comp = sub.groupby(["timepoint", "klass5"]).size().unstack(fill_value=0).reindex(index=TPS, columns=GU, fill_value=0)
xb = np.arange(len(TPS)); bottom = np.zeros(len(TPS))
for k in GU:
    axC.bar(xb, comp[k].values, bottom=bottom, color=CCOL[k], edgecolor="white", linewidth=0.4, label=CLAB[k], zorder=2)
    bottom += comp[k].values
axC.set_xticks(xb); axC.set_xticklabels([TPN[t] for t in TPS], fontsize=9); axC.set_ylabel("unique piRNAs", fontsize=9.5)
axC.legend(fontsize=7.6, frameon=False, loc="upper right")
axC.set_title("C  Class composition of subset", fontsize=10.5, fontweight="bold", loc="left")
axC.spines[["top", "right"]].set_visible(False)
# ---- Panel C2: % kept of genuinely-unique per tp (what the exact cut retains) ----
axC2 = fig.add_subplot(gs[2, 1])
kept = sub.groupby("timepoint").size().reindex(TPS, fill_value=0)
totg = full.groupby("timepoint").size().reindex(TPS, fill_value=0)
pct = (100 * kept / totg).values
axC2.bar(xb, pct, color=[TPCOL[t] for t in TPS], edgecolor="white", linewidth=0.4, zorder=2)
for xi, p, k, tg in zip(xb, pct, kept.values, totg.values):
    axC2.text(xi, p + 1.0, f"{p:.0f}%\n{k:,}/{tg:,}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
axC2.set_xticks(xb); axC2.set_xticklabels([TPN[t] for t in TPS], fontsize=9); axC2.set_ylabel("% of genuinely-unique kept", fontsize=9.5)
axC2.set_ylim(0, max(pct) * 1.35); axC2.set_title("D  Fraction retained by the exact cut", fontsize=10.5, fontweight="bold", loc="left")
axC2.spines[["top", "right"]].set_visible(False)
# ---- Panel C3: biological caveat text ----
axT = fig.add_subplot(gs[2, 2]); axT.axis("off")
axT.text(0, 1.0, "Biological note (BioMNI triple-verified)", fontsize=9.5, fontweight="bold", va="top")
axT.text(0, 0.86,
         "piRNA length = PIWI footprint:\n"
         "• MILI/MIWI2 (pre-pachytene) ~24–28 nt\n"
         "• MIWI (pachytene) ~29–32 nt\n\n"
         "Each footprint spans ~3 nt, so a window\n"
         "(26–28 / 29–31) captures the full peak;\n"
         "7 strains actually peak at 28 nt (MIWI2).\n\n"
         "This figure uses the user-specified EXACT\n"
         "27/30 nt (modal) cut: cleanest single-length\n"
         "stage-characteristic set, not the full\n"
         "pre-pachytene / pachytene populations.",
         fontsize=7.8, va="top", linespacing=1.45, color="#333")
n_sub, n_full = len(sub), len(full)
fig.suptitle(f"Stage-peak-length unique piRNAs (exact 27/30 nt) — {n_sub:,} of {n_full:,} genuinely-unique ({100*n_sub/n_full:.0f}%), 16 strains × 3 stages",
             fontsize=13, fontweight="bold", y=0.965)
fig.text(0.5, 0.018, "Subset of the genuinely-unique (conserved-but-silent + strain-private) within-tp set · exact stage lengths E16.5=27, P12.5=27&30, P20.5=30 nt · "
         "additional / non-destructive (full unique set unchanged) · biological caveat in panel.", ha="center", fontsize=7, color="#666")
out = f"{T}/figures/Fig_unique_stage_peak_lengths"
for e in ("pdf", "svg", "png"): fig.savefig(f"{out}.{e}", bbox_inches="tight")
# ---- source data ----
piv.assign(total=piv.sum(axis=1)).to_csv(f"{T}/data/source_data/SourceData_Fig_unique_stage_peak_lengths_perstrain.csv")
comp.assign(pct_of_unique=pct).to_csv(f"{T}/data/source_data/SourceData_Fig_unique_stage_peak_lengths_byclass.csv")
full.groupby(["timepoint", "L"]).size().unstack(fill_value=0).to_csv(f"{T}/data/source_data/SourceData_Fig_unique_stage_peak_lengths_distribution.csv")
print(f"wrote {out}.[pdf/svg/png]  (subset {n_sub:,} of {n_full:,} unique)")
print("per-tp kept:", {TPN[t]: int(kept[t]) for t in TPS})
