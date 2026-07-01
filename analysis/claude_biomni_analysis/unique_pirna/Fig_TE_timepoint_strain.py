#!/usr/bin/env python3
"""At EACH timepoint: how MANY strain-private piRNAs each strain makes (top bars) AND which TE families drive them
(heatmaps) — one combined message. Per strain, strain-private-locus piRNA loci (cand_self16 BAM) intersect the
per-strain RepeatMasker BED; each locus -> largest-overlap TE family, keeping timepoint (id = strain|tp|seq).
Layout: 3 columns (E16.5 / P12.5 / P20.5); each column = top bar (TOTAL # strain-private piRNAs per strain at that
window, log) + heatmap (TE family x 16 strains, # piRNAs, log colour, white = none). Reveals the prepachytene
TE-silencing wave (E16.5 huge, LTR/ERVK+LINE/L1) collapsing by P20.5. Wild strains (red) dominate. Per-locus TE
table cached -> re-plots are instant. Descriptive; framing queued for BioMNI."""
import warnings; warnings.filterwarnings("ignore")
import sys, subprocess, tempfile, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np, pysam
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPS = ["E16.5", "P12.5", "P20.5"]
TPCOL = {"E16.5": "#4393C3", "P12.5": "#FDB863", "P20.5": "#B2182B"}
d = pd.read_csv(f"{U}/unique16/final_classified_clean.csv.gz")   # mm0-clean strain-private (klass5)
# total strain-private piRNAs per (strain, timepoint) — the "unique piRNA by timepoint" counts (top bars)
ctall = (d[d.klass5 == "unique: strain-private locus"].assign(tp=lambda r: r.timepoint.map(TPMAP))
         .groupby(["strain", "tp"]).size().unstack(fill_value=0).reindex(index=CANON, columns=TPS).fillna(0))
cache = f"{PG}/SourceData_TE_timepoint_strain.csv"
if os.path.exists(cache):
    T = pd.read_csv(cache); print("loaded cached per-locus TE table")
else:
    rows = []
    for X in CANON:
        g = d[(d.strain == X) & (d.klass == "unique: strain-private locus")].copy(); g["id"] = X + "|" + g.timepoint + "|" + g.sequence
        sp = set(g.id); rm = f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
        if not sp or not os.path.exists(rm): continue
        bam = pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam", "rb")
        tb = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=PG); seen = set()
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.query_name not in sp or a.query_name in seen: continue
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
        for locid, (ov, fam) in best.items():
            rows.append((X, TPMAP.get(locid.split("|")[1], locid.split("|")[1]), fam))
        print(f"{X}: TE-annotated loci={len(best):,}")
    T = pd.DataFrame(rows, columns=["strain", "tp", "family"]); T.to_csv(cache, index=False)
topfam = T.family.value_counts().head(14).index.tolist()
mats = {tp: T[T.tp == tp].groupby(["family", "strain"]).size().unstack(fill_value=0).reindex(index=topfam, columns=CANON).fillna(0) for tp in TPS}
vmax = np.log10(max(m.values.max() for m in mats.values()) + 1); ymax = ctall.values.max() * 1.5
plt.rcParams.update({"font.family": "Liberation Sans"})
fig = plt.figure(figsize=(18, 7.8), dpi=300)
gs = fig.add_gridspec(2, 3, height_ratios=[1, 3.4], hspace=0.07, wspace=0.08)
cmap = plt.get_cmap("magma").copy(); cmap.set_bad("white"); heataxes = []
for j, tp in enumerate(TPS):
    axb = fig.add_subplot(gs[0, j]); axh = fig.add_subplot(gs[1, j]); heataxes.append(axh)
    axb.bar(range(len(CANON)), ctall[tp].values, width=0.8, color=TPCOL[tp], edgecolor="white", linewidth=0.3)
    axb.set_yscale("log"); axb.set_ylim(1, ymax); axb.set_xlim(-0.7, len(CANON) - 0.3); axb.set_xticks([])
    axb.tick_params(labelsize=6); axb.spines[["top", "right"]].set_visible(False)
    axb.set_ylabel("# private\npiRNAs (log)", fontsize=7) if j == 0 else axb.set_yticklabels([])
    axb.set_title(tp, fontsize=13, fontweight="bold", pad=4)
    for xi, X in enumerate(CANON):
        v = int(ctall[tp].values[xi])
        if v >= 500: axb.text(xi, v * 1.15, f"{v // 1000}k" if v >= 1000 else f"{v}", ha="center", va="bottom", fontsize=5, color="#333")
    Mt = mats[tp]; im = axh.imshow(np.ma.masked_where(Mt.values == 0, np.log10(Mt.values + 1)), aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
    axh.set_xticks(range(len(CANON))); axh.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=90, fontsize=6.5)
    for t, X in zip(axh.get_xticklabels(), CANON):
        if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
    axh.set_yticks(range(len(topfam))); axh.set_yticklabels(topfam if j == 0 else [], fontsize=7)
    axh.set_xticks(np.arange(-.5, len(CANON), 1), minor=True); axh.set_yticks(np.arange(-.5, len(topfam), 1), minor=True)
    axh.grid(which="minor", color="#e6e6e6", linewidth=0.5); axh.tick_params(which="minor", length=0)
    for i in range(len(topfam)):
        for k in range(len(CANON)):
            v = int(Mt.values[i, k])
            if v > 0: axh.text(k, i, f"{v:,}" if v < 1000 else f"{v // 1000}k", ha="center", va="center", fontsize=3.8, color="white" if np.log10(v + 1) < vmax * 0.6 else "black")
fig.colorbar(im, ax=heataxes, fraction=0.012, pad=0.012).set_label("log10(# strain-private piRNAs +1)", fontsize=8)
fig.suptitle("Strain-private piRNAs across development — how MANY (top bars) and which TE families drive them (heatmaps), 16 strains × 3 windows\n"
             "prepachytene (E16.5) wave: wild strains make 1000s, dominated by LTR/ERVK + LINE/L1, collapsing by P20.5; wild strains (red) carry far more throughout",
             fontsize=11.5, fontweight="bold", y=1.01)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_TE_timepoint_strain.{e}", bbox_inches="tight")
print("wrote Fig_TE_timepoint_strain (with timepoint count bars)")
