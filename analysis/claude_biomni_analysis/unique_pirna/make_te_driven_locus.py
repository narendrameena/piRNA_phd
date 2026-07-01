#!/usr/bin/env python3
"""Figure for a VERIFIED coordinate-based TE-DRIVEN piRNA cluster: a strain-PRIVATE TE insertion that created an
antisense-1U silencing cluster (absent from every strain lacking the insertion). DATA ONLY: sRNA BAMs (PRIMARY
reads), RepeatMasker TE annotation (own genome), pangenome (private = singleton in the 16-strain VCF). NO liftover.
  (A) developmental activation: per-timepoint piRNA reads (log) with %1U and %antisense-to-TE;
  (B) the locus in the carrier strain: primary-read coverage per timepoint (plus ↑ green / minus ↓ purple) + the
      TE-insertion track (the provirus/element that drove it); strain-private;
  (C) nucleotide resolution at the silencing peak: piRNAs as coloured bases; 5'-U = 1U; RED 5' arrow = ANTISENSE-to-TE.
Usage: make_te_driven_locus.py <strain> <chrom> <start> <end> <TE_label> <n_absent> <outbase>"""
import sys, os, numpy as np
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from collections import Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, fetch_primary, te_at
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; PG = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb", "SINE/MIR": "#1b9e77"}
X, CH, PS, PE, TELAB, NABS, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6]), sys.argv[7]
nb = 220; N = max(1, PE - PS)
DAT = {tp: fetch_primary(X, CH, PS, PE, tp, nb) for tp in TPS}
tes = next((DAT[tp]["tes"] for tp in TPS if DAT[tp]["ntot"]) , te_at(X, CH, PS, PE))
def tst(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
CHOSEN = max(TPS, key=lambda tp: DAT[tp]["ntot"])
print(f"{X} chr{CH}:{PS:,}-{PE:,}  reads/tp=" + " ".join(f"{TPLAB[t]}:{DAT[t]['ntot']}" for t in TPS) + f"  peak={TPLAB[CHOSEN]}")
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(12.5, 11), dpi=300); gs = fig.add_gridspec(3, 1, height_ratios=[0.8, 2.1, 1.4], hspace=0.6); fig.subplots_adjust(top=0.9, bottom=0.06)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
# A: developmental activation (per-tp reads, log) + 1U + antisense-to-TE
xb = np.arange(3); tot = [DAT[t]["ntot"] for t in TPS]
axA.bar(xb, [max(v, 0.5) for v in tot], width=0.6, color=[TPCOL[t] for t in TPS], edgecolor="white")
axA.set_yscale("log"); axA.set_ylim(0.5, max(tot) * 3.5); axA.set_xticks(xb); axA.set_xticklabels([TPLAB[t] for t in TPS], fontsize=9)
for i, t in enumerate(TPS):
    d = DAT[t]; at = 100 * d["nat"] / max(1, d["nte"]); u = 100 * d["n1u"] / max(1, d["ntot"])
    axA.text(i, max(d["ntot"], 0.5) * 1.25, f"{d['ntot']:,}\n{u:.0f}% 1U\n{at:.0f}% AS→TE" + ("  ● peak" if t == CHOSEN else ""), ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=TPCOL[t], linespacing=1.3)
axA.set_ylabel("primary piRNA\nreads (log)", fontsize=8.5); axA.spines[["top", "right"]].set_visible(False)
axA.set_title(f"A  Developmental activation of the TE-driven cluster — ● = peak timepoint (nucleotide panel C)", fontsize=9.2, fontweight="bold", loc="left")
# B: carrier-strain per-timepoint coverage block — canonical multi-figure style (shared helper) + TE + gene tracks
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()   # 0–1 fractional x (like the locus figures); lims fixed before pc.rtext()
five = Counter((r[1] - 1 if r[2] else r[0]) for r in DAT[CHOSEN]["reads"]); pk = five.most_common(1)[0][0] if five else PS + N // 2; z0, z1 = pk - 30, pk + 50
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, DAT, tes, CH, PS, PE, None, 0.0, z0, z1, name=X, is_top=True, wild=(X in {"CAST_EiJ", "PWK_PhJ", "SPRET_EiJ", "WSB_EiJ"}), TECOL=TECOL, nb=nb)
_fams = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for _, _, _, f in tes))[:5]
axB.legend(handles=[Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#4a4a4a", label="+ strand = deeper shade"), Patch(facecolor="#cfcfcf", label="− strand = paler shade"), Patch(facecolor="#c9d6ea", edgecolor=pc.C_GENE, label="gene model")] + [Patch(facecolor=TECOL.get(f, "#ddd"), label=f) for f in _fams], fontsize=6.0, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), columnspacing=1.0, handlelength=1.0)
pc.pbadge(axB, "B", f"{X.replace('_','/')} chr{CH}:{PS:,}–{PE:,} — per-timepoint sRNA coverage (height = expression, colour = timepoint) above TE + gene tracks · {_ntot:,} primary piRNAs · {_p1u:.0f}% 1U", fs=7.0, y=1.05)
axB.text(0.012, 0.985, f"{X.replace('_','/')}-PRIVATE insertion (singleton in 16-strain pangenome VCF · absent from the other {NABS} genomes) → a piRNA cluster cannot exist where the TE is absent = TE-DRIVEN", transform=axB.transAxes, fontsize=5.6, color=pc.C_WILD, fontweight="bold", ha="left", va="center")
# C: nucleotide at the silencing peak (z0,z1 + zoom box already set in Panel B via draw_strain_block)
d = DAT[CHOSEN]; reads = d["reads"]
def anti_te(rs, re, isrev):
    s = tst((rs + re) // 2); return None if s is None else ((s == "-") != isrev)
zr = [r for r in reads if r[0] < z1 and r[1] > z0]; pr, mr = Counter(), Counter()
for rs, re, isrev, seq in zr: (mr if isrev else pr)[(rs, re, seq, isrev)] += 1
def draw(items, y0, dirn):
    y = y0
    for (rs, re, seq, isrev), cnt in items:
        for k, ch in enumerate(seq):
            xx = rs + k
            if z0 - 2 <= xx <= z1 + 2: axC.text(xx, y, ch, fontsize=4.6, ha="center", va="center", family="monospace", color="white", bbox=dict(boxstyle="square,pad=0.02", fc=NT.get(ch, "#999"), ec="none"))
        a = anti_te(rs, re, isrev); acol = "#C0392B" if a else ("#888" if a is not None else "#ccc"); fp = re - 1 if isrev else rs
        axC.annotate("", xy=(fp + (0.5 if not isrev else -0.5), y), xytext=(fp + (-2.6 if not isrev else 2.6), y), arrowprops=dict(arrowstyle="-|>", color=acol, lw=1.0))
        axC.text(z1 + 4, y, f"×{cnt}", fontsize=4.6, va="center", color="#666"); y += dirn
    return y
ytop = draw(pr.most_common(7), 1, 1); axC.axhline(0, color="#333", lw=0.7); ybot = draw(mr.most_common(7), -1, -1)
axC.set_xlim(z0 - 1, z1 + 8); axC.set_ylim(ybot - 1.6, ytop + 0.6)
for sp in ("top", "left", "right"): axC.spines[sp].set_visible(False)
axC.set_yticks([]); axC.spines["bottom"].set_position(("data", ybot - 1.0))
tk = np.linspace(z0, z1, 5).astype(int); axC.set_xticks(tk); axC.set_xticklabels([f"{t:,}" for t in tk], fontsize=6.5); axC.tick_params(axis="x", length=3)
axC.set_xlabel(f"{X.replace('_','/')} chr{CH} position (bp)", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
for x in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=((x - PS) / N, -1.30), coordsA=axB.transData, xyB=(x, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
axC.set_title(f"C  Base resolution at {TPLAB[CHOSEN]} silencing peak; 5′-U = 1U; 5′ arrow RED = ANTISENSE-to-TE (silencing), grey = sense-to-TE", fontsize=8.0, fontweight="bold", loc="left")
fig.suptitle(f"TE-DRIVEN piRNA cluster — {X.replace('_','/')} {TELAB} (coordinate-verified)", fontsize=12.5, fontweight="bold", y=0.965)
import os as _os, pandas as _pd; _os.makedirs("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/09_TE_driven_evolution/data/source_data",exist_ok=True); _pd.DataFrame([{"locus":f"{CH}:{PS}-{PE}","carrier":X,"TE":TELAB,"timepoint":TPLAB[_t],"reads_ntot":DAT[_t]["ntot"],"pct_1U":round(100*DAT[_t]["n1u"]/max(1,DAT[_t]["ntot"]),2),"pct_antisense_to_TE":round(100*DAT[_t]["nat"]/max(1,DAT[_t]["nte"]),2)} for _t in TPS]).to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/09_TE_driven_evolution/data/source_data/SourceData_{OUT}.csv",index=False)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
