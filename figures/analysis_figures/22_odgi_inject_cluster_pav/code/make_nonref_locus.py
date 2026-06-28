#!/usr/bin/env python3
"""THEME 22 — NON-REFERENCE locus figure: a piRNA cluster GENUINELY ABSENT from GRCm39, shown in the strain that carries
it, with the pangenome GRAPH (odgi inject + pav) as the cross-strain evidence.
  (A) GRAPH PAV (odgi): per-strain coverage of the cluster's nodes in graph_inj.og — the cluster's SEQUENCE is present
      in N strains but GRCm39 = 0 (genuinely absent from the reference path; controls validate the group, step 20/21);
  (B) the carrier strain's per-timepoint sRNA coverage (height = expression, colour = timepoint) + TE + gene tracks;
  (C) base resolution. Reuses the pav_clusters drawing blocks (strain-anchored) — only panel A is graph-native.
Usage: make_nonref_locus.py <strain> <chrom> <start> <end> <graph_pav_name> <title> <outbase>"""
import sys, os, numpy as np, pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
sys.path.insert(0, f"{ROOT}/analysis/claude_biomni_analysis"); sys.path.insert(0, f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, fetch_primary, te_at, dom_te_family
T22 = f"{ROOT}/figures/analysis_figures/22_odgi_inject_cluster_pav"; D = f"{T22}/data"; CP = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
STRAIN, CH, ps, pe, GPNAME, TITLE, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}; COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "SINE/B4": "#fb9a99", "Simple_repeat": "#bbbbbb"}
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; N = max(1, pe - ps); nb = 200
DAT = {}
for tp in TPS:
    d = fetch_primary(STRAIN, CH, ps, pe, tp, nb); DAT[tp] = d if (d and d["ntot"] > 0) else None
tes = next((DAT[tp]["tes"] for tp in TPS if DAT[tp]), te_at(STRAIN, CH, ps, pe)); domTE = dom_te_family(tes, ps, pe)
CHOSEN = max((tp for tp in TPS if DAT[tp]), key=lambda tp: DAT[tp]["ntot"], default=TPS[0])
fpm = pd.read_csv(f"{CP}/{STRAIN}.clusters_fpm.bed", sep="\t", header=None, names=["c", "s", "e", "allF", "uniqF", "st", "tp"], dtype={"c": str})
mm = fpm[(fpm.c == str(CH)) & (fpm.s < pe) & (fpm.e > ps)]; fpm_tp = pd.Series({tp: float(mm[mm.tp == tp].allF.sum()) for tp in TPS})
gp = pd.read_csv(f"{D}/graph_check_pav3.tsv", sep="\t"); grow = gp[gp.name == GPNAME]
cov = {s: (float(grow.iloc[0].get(s, 0.0)) if len(grow) else 0.0) for s in ORDER}; gr_cov = float(grow.iloc[0].get("GRCm39", 0.0)) if len(grow) else 0.0
npres = sum(v >= 0.5 for v in cov.values())
print(f"{TITLE}: {STRAIN} chr{CH}:{ps:,}-{pe:,} domTE={domTE} FPM_tp={dict(fpm_tp.round(1))} | graph: {npres}/16 strains>=0.5, GRCm39={gr_cov:.2f}")
def testr(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(12.5, 11.8), dpi=300); gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 2.3, 1.4], hspace=0.62); fig.subplots_adjust(top=0.9, bottom=0.06)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
# === A: GRAPH PAV (odgi) — per-strain coverage; GRCm39 absent ===
xs = np.arange(len(ORDER))
axA.bar(xs, [cov[s] for s in ORDER], width=0.74, color=["#C0392B" if s in WILD else "#4393C3" for s in ORDER], edgecolor="white")
axA.bar([len(ORDER) + 0.6], [gr_cov], width=0.74, color="#222", edgecolor="white")
axA.axhline(0.5, ls=(0, (3, 2)), color="#999", lw=0.8)
axA.text(len(ORDER) + 0.6, max(gr_cov, 0) + 0.04, "GRCm39\nABSENT" if gr_cov < 0.1 else "GRCm39", ha="center", va="bottom", fontsize=6.6, color="#C0392B" if gr_cov < 0.1 else "#222", fontweight="bold")
axA.annotate(f"carrier: {STRAIN.replace('_','/')}\n{fpm_tp.max():.0f} FPM", xy=(list(ORDER).index(STRAIN), cov[STRAIN]), xytext=(list(ORDER).index(STRAIN), 1.12), fontsize=6.6, ha="center", color="#1B7837", fontweight="bold", arrowprops=dict(arrowstyle="-", color="#1B7837", lw=0.7))
axA.set_xticks(list(xs) + [len(ORDER) + 0.6]); axA.set_xticklabels([s.replace("_", "/") for s in ORDER] + ["GRCm39"], rotation=90, fontsize=6.2)
[axA.get_xticklabels()[i].set_color("#C0392B") for i, s in enumerate(ORDER) if s in WILD]
axA.set_ylim(0, 1.25); axA.set_yticks([0, 0.5, 1.0]); axA.set_ylabel("graph coverage\n(odgi pav, fraction of nodes)", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
pc.pbadge(axA, "A", f"GRAPH PAV (odgi inject + pav, graph_inj.og) — the cluster's sequence is present in {npres}/16 strains but GENUINELY ABSENT from GRCm39 (=0); dashed = 0.5 presence threshold", fs=8.0, y=1.18)
# === B: strain sRNA coverage + TE + gene (reuse draw_strain_block) ===
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()
_dc = DAT[CHOSEN]; reads = _dc["reads"] if _dc else []
_f5 = Counter(r[0] for r in reads if not r[2]); _r5 = Counter(r[1] - 1 for r in reads if r[2]); _rn = lambda i: sum(_r5[j] for j in range(i - 2, i + 13))
_cd = [(min(_f5[i], _rn(i)), _f5[i] + _rn(i), i + 5) for i in _f5 if _rn(i) > 0]
pk = max(_cd)[2] if _cd else ((Counter((r[1] - 1 if r[2] else r[0]) for r in reads).most_common(1)[0][0]) if reads else ps + N // 2); z0, z1 = pk - 30, pk + 50
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, DAT, tes, CH, ps, pe, fpm_tp, 0.0, z0, z1, name=STRAIN, is_top=True, wild=(STRAIN in WILD), TECOL=TECOL, only_tp=CHOSEN)
_fams = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for _, _, _, f in tes))[:5]
_lh = [Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#6a3d9a", label="solid = + strand"), Patch(facecolor=pc.pale("#6a3d9a", 0.55), label="pale = − strand"), Patch(facecolor="#efefef", label="non-TE piRNA"), Patch(facecolor="#C0392B", label="antisense-to-TE = silencing"), Patch(facecolor="#cfcfcf", label="sense-to-TE"), Patch(facecolor="#c9d6ea", edgecolor=pc.C_GENE, label="gene model")] + [Patch(facecolor=pc.famcol(f), label=f) for f in _fams]
pc.color_legend(axB.legend(handles=_lh, fontsize=6.0, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), columnspacing=1.0, handlelength=1.0), _lh)
pc.pbadge(axB, "B", f"{STRAIN.replace('_','/')} chr{CH}:{ps:,}–{pe:,} — per-timepoint sRNA coverage (height = expression, colour = timepoint) above TE + gene tracks · {_ntot:,} primary piRNAs · domTE {domTE}", fs=7.4, y=1.07)
# === C: base resolution (reuse) ===
def anti_te(rs, re, isrev):
    st = testr((rs + re) // 2); return None if st is None else ((st == "-") != isrev)
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
axC.set_xlabel(f"{STRAIN.replace('_','/')} chr{CH} position (bp) — every base at its true genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
for x in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=((x - ps) / N, -1.30), coordsA=axB.transData, xyB=(x, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
pc.pbadge(axC, "C", f"Base resolution, {STRAIN.replace('_','/')} at {TPLAB[CHOSEN]} (top tp)   ·   5′ arrow RED = antisense-to-TE (silencing), grey = sense-to-TE", fs=7.6)
fig.suptitle(TITLE, fontsize=12, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{T22}/figures/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
