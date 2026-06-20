#!/usr/bin/env python3
"""TE families across development for the TWO genuinely-unique subcategories (klass5 ≥2-read): strain-private
(locus NEW) vs conserved-but-silent (locus SHARED). Per strain, each subcategory's piRNA loci (cand_self16 BAM)
intersect the per-strain RepeatMasker BED -> largest-overlap TE family, keeping timepoint (id=strain|tp|seq).
Layout: 3 columns (E16.5/P12.5/P20.5); row 0 = total piRNAs per strain at that window (grouped: strain-private vs
conserved-but-silent, log); row 1 = strain-private TE-family heatmap; row 2 = conserved-but-silent TE-family heatmap.
Reveals the prepachytene TE-silencing wave (E16.5 large, LTR/ERVK+LINE/L1) and how it differs between subcategories.
Per-locus by-class TE table cached. Descriptive; framing queued for BioMNI."""
import warnings; warnings.filterwarnings("ignore")
import sys, subprocess, tempfile, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
import pandas as pd, numpy as np, pysam
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPS = ["E16.5", "P12.5", "P20.5"]
SUB = [("strain-private", "unique: strain-private locus", "#7a3b9a"), ("conserved-but-silent", "unique: conserved-but-silent", "#0072B2")]
d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")
# total piRNAs per (strain, tp, klass) — top bars (from the classification directly)
d2 = d[d.klass5.isin([k for _, k, _ in SUB])].assign(tp=lambda r: r.timepoint.map(TPMAP), kl=lambda r: r.klass5.map({k: lab for lab, k, _ in SUB}))
ctall = {lab: d2[d2.kl == lab].groupby(["strain", "tp"]).size().unstack(fill_value=0).reindex(index=CANON, columns=TPS).fillna(0) for lab, _, _ in SUB}
cache = f"{PG}/SourceData_TE_timepoint_strain_byclass.csv"
if os.path.exists(cache):
    T = pd.read_csv(cache); print("loaded cached by-class per-locus TE table")
else:
    rows = []
    for X in CANON:
        g = d[(d.strain == X) & (d.klass5.isin([k for _, k, _ in SUB]))]
        cls = {X + "|" + r.timepoint + "|" + r.sequence: ("strain-private" if r.klass5.endswith("private locus") else "conserved-but-silent") for r in g.itertuples()}
        rm = f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
        if not cls or not os.path.exists(rm): continue
        bam = pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam", "rb")
        tb = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=PG); seen = set()
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.query_name not in cls or a.query_name in seen: continue
            seen.add(a.query_name); tb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
        tb.close(); bam.close()
        out = subprocess.run(f"sort -k1,1 -k2,2n {tb.name} | {BT} intersect -a - -b {rm} -wa -wb", shell=True, capture_output=True, text=True).stdout
        os.unlink(tb.name)
        best = {}
        for ln in out.splitlines():
            f = ln.split("\t")
            if len(f) < 8: continue
            ov = min(int(f[2]), int(f[6])) - max(int(f[1]), int(f[5])); fam = f[7].split("|")[-1] if "|" in f[7] else f[7]
            if f[3] not in best or ov > best[f[3]][0]: best[f[3]] = (ov, fam)
        agg = Counter()
        for cid, (ov, fam) in best.items(): agg[(cls[cid], TPMAP.get(cid.split("|")[1], cid.split("|")[1]), fam)] += 1
        for (klab, tp, fam), n in agg.items(): rows.append((klab, tp, fam, X, n))
        print(f"{X}: TE-annotated loci={len(best):,}")
    T = pd.DataFrame(rows, columns=["klass", "tp", "family", "strain", "count"]); T.to_csv(cache, index=False)
# top families pooled over both classes + all tp
topfam = T.groupby("family")["count"].sum().sort_values(ascending=False).head(14).index.tolist()
def mat(klab, tp):
    s = T[(T.klass == klab) & (T.tp == tp)]
    return s.pivot_table(index="family", columns="strain", values="count", aggfunc="sum").reindex(index=topfam, columns=CANON).fillna(0)
mats = {(klab, tp): mat(klab, tp) for lab_ in [0] for klab, _, _ in SUB for tp in TPS}
vmax = np.log10(max(m.values.max() for m in mats.values() if m.values.size) + 1)
ymax = max(ctall[lab].values.max() for lab, _, _ in SUB) * 1.6
plt.rcParams.update({"font.family": "Liberation Sans"})
cmap = LinearSegmentedColormap.from_list("vivBlue", ["#eaf3fb", "#9ecae8", "#3a8fd4", "#1565a8"]); cmap.set_bad("white")
fig = plt.figure(figsize=(18, 11.5), dpi=300)
gs = fig.add_gridspec(3, 3, height_ratios=[1, 2.4, 2.4], hspace=0.14, wspace=0.08)
heataxes = []
for j, tp in enumerate(TPS):
    axb = fig.add_subplot(gs[0, j]); x = np.arange(len(CANON)); ww = 0.4
    for s_i, (lab, _, col) in enumerate(SUB):
        axb.bar(x + (s_i - 0.5) * ww, ctall[lab][tp].values, ww, color=col, edgecolor="white", linewidth=0.3, label=lab if j == 0 else None)
    axb.set_yscale("log"); axb.set_ylim(1, ymax); axb.set_xlim(-0.5, len(CANON) - 0.5); axb.set_xticks([])
    axb.tick_params(labelsize=6); axb.spines[["top", "right"]].set_visible(False)
    axb.set_ylabel("# piRNAs (log)", fontsize=7) if j == 0 else axb.set_yticklabels([])
    axb.set_title(tp, fontsize=13, fontweight="bold", pad=4)
    if j == 0: axb.legend(fontsize=6.2, frameon=False, loc="upper left")
    for ri, (klab, _, _) in enumerate(SUB):
        axh = fig.add_subplot(gs[ri + 1, j]); heataxes.append(axh); Mt = mats[(klab, tp)]
        axh.imshow(np.ma.masked_where(Mt.values == 0, np.log10(Mt.values + 1)), aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
        axh.set_yticks(range(len(topfam))); axh.set_yticklabels(topfam if j == 0 else [], fontsize=7)
        if ri == len(SUB) - 1:
            axh.set_xticks(range(len(CANON))); axh.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=90, fontsize=6.5)
            for t, X in zip(axh.get_xticklabels(), CANON):
                if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
        else: axh.set_xticks([])
        axh.set_xticks(np.arange(-.5, len(CANON), 1), minor=True); axh.set_yticks(np.arange(-.5, len(topfam), 1), minor=True)
        axh.grid(which="minor", color="#e6e6e6", linewidth=0.5); axh.tick_params(which="minor", length=0)
        for i in range(len(topfam)):
            for k in range(len(CANON)):
                v = int(Mt.values[i, k])
                if v > 0: axh.text(k, i, f"{v:,}" if v < 1000 else f"{v//1000}k", ha="center", va="center", fontsize=3.6, color="white" if (np.log10(v+1)/vmax if vmax else 0) > 0.55 else "#222")
        if j == 0: axh.text(-0.42, 0.5, ["strain-private\n(locus NEW)", "conserved-but-silent\n(locus SHARED)"][ri], transform=axh.transAxes, rotation=90, va="center", ha="center", fontsize=8.5, fontweight="bold", color=SUB[ri][2])
im = heataxes[0].images[0]
fig.colorbar(im, ax=heataxes, fraction=0.012, pad=0.012).set_label("log10(# piRNAs +1)", fontsize=8)
fig.suptitle("TE families across development, split by genuinely-unique subcategory — strain-private (locus NEW) vs conserved-but-silent (locus SHARED), 16 strains × 3 windows\n"
             "top = how many piRNAs per strain (grouped by subcategory); heatmaps = which TE families drive them; wild strains (red) carry far more", fontsize=10.5, fontweight="bold", y=0.995)
# classical(blue)/wild(orange) companion below the bottom-left heatmap: total genuinely-unique piRNAs per strain
_tot=d[d.klass5.str.startswith("unique")].groupby("strain").size().reindex(CANON).fillna(0).values
_cax=add_classical_wild_companion(fig, heataxes[1], CANON, _tot, gap=0.085, height_frac=0.45, ylabel="total\nuniq (log)")
_cax.set_xticks(np.arange(len(CANON))); _cax.set_xticklabels([s.replace("_","/") for s in CANON], rotation=90, fontsize=6)
for lab,s in zip(_cax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — total genuinely-unique piRNAs per strain", fontsize=7, fontweight="bold", loc="left")
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_TE_timepoint_strain.{e}", bbox_inches="tight")
print("wrote Fig_TE_timepoint_strain (by subcategory)")
