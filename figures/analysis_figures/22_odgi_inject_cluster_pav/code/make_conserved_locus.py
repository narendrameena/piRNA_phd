#!/usr/bin/env python3
"""THEME 22 — CONSERVED-cluster figure: WHY coordinate liftover fails (the textbook case for a pangenome). A major piRNA
cluster whose sequence is present (minimap2 qcov~1) in ALL strains + GRCm39 but at STRAIN-SPECIFIC coordinates, so
halLiftover mis-projects it and calls it 'non-reference'. minimap2 (sequence homology, not coordinate) is the arbiter
here because the region is also a pangenome-graph gap.
  (A) per-strain FPM at the homologous locus (expression) + a presence/COORDINATE strip (minimap2 qcov + each genome's
      aligned position, making the coordinate offset explicit); (B) carrier sRNA/TE/gene tracks; (C) base resolution.
Usage: make_conserved_locus.py <strain> <chrom> <start> <end> <conserved_tsv> <title> <outbase>"""
import sys, numpy as np, pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
sys.path.insert(0, f"{ROOT}/analysis/claude_biomni_analysis"); sys.path.insert(0, f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, fetch_primary, te_at, dom_te_family
T22 = f"{ROOT}/figures/analysis_figures/22_odgi_inject_cluster_pav"; CP = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
STRAIN, CH, ps, pe, CONS, TITLE, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7]
NT = {"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}; COMP = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
TECOL = {"LINE/L1":"#E69F00","LTR/ERVK":"#6a3d9a","LTR/ERVL-MaLR":"#b15928","SINE/Alu":"#33a02c","SINE/B2":"#1f78b4","LTR/ERVL":"#a6cee3","LTR/ERV1":"#cab2d6","SINE/B4":"#fb9a99","Simple_repeat":"#bbbbbb"}
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; N = max(1, pe - ps); nb = 200
_fc = {}
def fpm_at(X, c, s, e):
    if X not in _fc: _fc[X] = pd.read_csv(f"{CP}/{X}.clusters_fpm.bed", sep="\t", header=None, names=["c","s","e","aF","uF","st","tp"], dtype={"c":str})
    m = _fc[X][(_fc[X].c==str(c)) & (_fc[X].s<e) & (_fc[X].e>s)]; return {tp: float(m[m.tp==tp].aF.sum()) for tp in TPS}
cons = pd.read_csv(CONS, sep="\t", header=None, names=["strain","qlen","qcov","tc","ts","te","mapq"]).set_index("strain")
FPM = pd.DataFrame(0.0, index=ORDER, columns=TPS); QCOV = {}; COORD = {}
for Y in ORDER + ["GRCm39"]:
    if Y in cons.index and pd.notna(cons.loc[Y, "qcov"]):
        r = cons.loc[Y]; QCOV[Y] = float(r.qcov); _tc = str(r.tc).split(".")[0]; COORD[Y] = (_tc, int(r.ts))
        if Y in ORDER:
            for tp, v in fpm_at(Y, _tc, int(r.ts), int(r.te)).items(): FPM.loc[Y, tp] = v
    else: QCOV[Y] = 0.0; COORD[Y] = None
npres = sum(QCOV[Y] >= 0.5 for Y in ORDER + ["GRCm39"])
DAT = {tp: (d if (d := fetch_primary(STRAIN, CH, ps, pe, tp, nb)) and d["ntot"] > 0 else None) for tp in TPS}
tes = next((DAT[tp]["tes"] for tp in TPS if DAT[tp]), te_at(STRAIN, CH, ps, pe)); domTE = dom_te_family(tes, ps, pe)
CHOSEN = max((tp for tp in TPS if DAT[tp]), key=lambda tp: DAT[tp]["ntot"], default=TPS[0])
print(f"{TITLE}: {STRAIN} chr{CH}:{ps:,}-{pe:,} domTE={domTE} | minimap2 present (qcov>=0.5) in {npres}/17; GRCm39 qcov={QCOV.get('GRCm39',0):.2f} @ {COORD.get('GRCm39')}")
def testr(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig = plt.figure(figsize=(12.6, 12.6), dpi=300); gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 0.42, 2.15, 1.3], hspace=0.55); fig.subplots_adjust(top=0.9, bottom=0.06)
axA = fig.add_subplot(gs[0]); axAc = fig.add_subplot(gs[1], sharex=axA); axB = fig.add_subplot(gs[2]); axC = fig.add_subplot(gs[3])
G = ORDER + ["GRCm39"]; x = np.arange(len(ORDER)); gx = len(ORDER) + 0.4; bw = 0.27
# A: per-strain FPM (log, by tp) at the homologous locus
for j, tp in enumerate(TPS):
    h = FPM[tp].values
    axA.bar(x + (j - 1) * bw, np.maximum(h, 1e-3), width=bw, color=TPCOL[tp], edgecolor="white", linewidth=0.2, label=TPLAB[tp])
    for xi in range(len(ORDER)):
        if h[xi] > 0: axA.text(xi + (j - 1) * bw, h[xi] * 1.2, (f"{h[xi]:.0f}" if h[xi] >= 1 else f"{h[xi]:.1f}"), ha="center", va="bottom", fontsize=4.2, rotation=90, color=TPCOL[tp], fontweight="bold")
axA.bar([gx], [1e-3], width=bw, color="#ddd")
axA.set_yscale("log"); axA.set_ylim(0.1, max(FPM.values.max(), 1) * 24); axA.set_ylabel("PICB cluster\nFPM (log)", fontsize=8)
axA.spines[["top","right"]].set_visible(False); plt.setp(axA.get_xticklabels(), visible=False)
axA.legend(fontsize=6, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.004, 1.02), title="timepoint", title_fontsize=6.5)
pc.pbadge(axA, "A", "Per-strain FPM at the HOMOLOGOUS locus (expression) — the cluster is expressed across the panel; below: it is PRESENT everywhere but coordinate-shifted", fs=7.6, y=1.18)
# A-cov: minimap2 presence (qcov) + the aligned COORDINATE (the offset that breaks liftover)
for xi, Y in enumerate(ORDER):
    axAc.bar(xi, QCOV.get(Y, 0), width=0.74, color=("#C0392B" if Y in WILD else "#4393C3"), edgecolor="white", linewidth=0.3)
    if COORD.get(Y): axAc.text(xi, 1.06, f"{COORD[Y][0]}:{COORD[Y][1]/1e6:.1f}", ha="center", va="bottom", fontsize=4.0, rotation=90, color="#555")
axAc.bar(gx, QCOV.get("GRCm39", 0), width=0.74, color="#222", edgecolor="white")
if COORD.get("GRCm39"): axAc.text(gx, 1.06, f"{COORD['GRCm39'][0]}:{COORD['GRCm39'][1]/1e6:.1f}", ha="center", va="bottom", fontsize=4.4, rotation=90, color="#C0392B", fontweight="bold")
axAc.axhline(0.5, ls=(0,(2,2)), color="#bbb", lw=0.6); axAc.set_ylim(0, 1.5); axAc.set_yticks([0,1]); axAc.tick_params(labelsize=6); axAc.set_ylabel("minimap2\nqcov", fontsize=6.6); axAc.spines[["top","right"]].set_visible(False)
axAc.set_xticks(list(x) + [gx]); axAc.set_xticklabels([s.replace("_","/") for s in ORDER] + ["GRCm39"], rotation=90, fontsize=7.0)
for t, Y in zip(axAc.get_xticklabels()[:len(ORDER)], ORDER):
    if Y in WILD: t.set_color("#C0392B")
axAc.get_xticklabels()[-1].set_color("#C0392B"); axAc.get_xticklabels()[-1].set_fontweight("bold")
axAc.text(0.0, 1.62, f"minimap2 presence (qcov) + ALIGNED COORDINATE per genome — CONSERVED (present in {npres}/17) but at strain-specific coordinates "
          f"(carrier {STRAIN.replace('_','/')} chr{CH}:{ps/1e6:.1f}Mb ↔ GRCm39 {COORD.get('GRCm39',['?','?'])[0] if COORD.get('GRCm39') else '?'}:{COORD['GRCm39'][1]/1e6:.1f}Mb) → halLiftover fails = false 'non-reference'",
          transform=axAc.transAxes, fontsize=6.2, color="#555", fontweight="bold", ha="left", va="bottom")
# (subtitle removed — the suptitle already states conserved + coordinate-offset + liftover-failure)
# B: carrier coverage + TE + gene (reuse)
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()
_dc = DAT[CHOSEN]; reads = _dc["reads"] if _dc else []
_f5 = Counter(r[0] for r in reads if not r[2]); _r5 = Counter(r[1]-1 for r in reads if r[2]); _rn = lambda i: sum(_r5[j] for j in range(i-2, i+13))
_cd = [(min(_f5[i], _rn(i)), _f5[i]+_rn(i), i+5) for i in _f5 if _rn(i) > 0]
pk = max(_cd)[2] if _cd else ((Counter((r[1]-1 if r[2] else r[0]) for r in reads).most_common(1)[0][0]) if reads else ps + N//2); z0, z1 = pk-30, pk+50
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, DAT, tes, CH, ps, pe, FPM.loc[STRAIN], 0.0, z0, z1, name=STRAIN, is_top=True, wild=(STRAIN in WILD), TECOL=TECOL, only_tp=CHOSEN)
_fams = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for _,_,_,f in tes))[:5]
_lh = [Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#6a3d9a", label="solid = + strand"), Patch(facecolor=pc.pale("#6a3d9a", 0.55), label="pale = − strand"), Patch(facecolor="#efefef", label="non-TE piRNA"), Patch(facecolor="#C0392B", label="antisense-to-TE = silencing"), Patch(facecolor="#cfcfcf", label="sense-to-TE"), Patch(facecolor="#c9d6ea", edgecolor=pc.C_GENE, label="gene model")] + [Patch(facecolor=pc.famcol(f), label=f) for f in _fams]
pc.color_legend(axB.legend(handles=_lh, fontsize=6.0, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), columnspacing=1.0, handlelength=1.0), _lh)
pc.pbadge(axB, "B", f"{STRAIN.replace('_','/')} chr{CH}:{ps:,}–{pe:,} — per-timepoint sRNA coverage above TE + gene tracks · {_ntot:,} primary piRNAs · domTE {domTE}", fs=7.4, y=1.07)
# C: base resolution (reuse)
def anti_te(rs, re, isrev):
    st = testr((rs+re)//2); return None if st is None else ((st=="-") != isrev)
zr = [r for r in reads if r[0] < z1 and r[1] > z0]; pr, mr = Counter(), Counter()
for rs, re, isrev, seq in zr: (mr if isrev else pr)[(rs, re, seq, isrev)] += 1
def draw(items, y0, dirn):
    y = y0
    for (rs, re, seq, isrev), cnt in items:
        for k, ch in enumerate(seq):
            xx = rs + k
            if z0-2 <= xx <= z1+2: axC.text(xx, y, ch, fontsize=4.6, ha="center", va="center", family="monospace", color="white", bbox=dict(boxstyle="square,pad=0.02", fc=NT.get(ch,"#999"), ec="none"))
        a = anti_te(rs, re, isrev); acol = "#C0392B" if a else ("#888" if a is not None else "#ccc"); fp = re-1 if isrev else rs
        axC.annotate("", xy=(fp+(0.5 if not isrev else -0.5), y), xytext=(fp+(-2.6 if not isrev else 2.6), y), arrowprops=dict(arrowstyle="-|>", color=acol, lw=1.0))
        axC.text(z1+4, y, f"×{cnt}", fontsize=4.6, va="center", color="#666"); y += dirn
    return y
ytop = draw(pr.most_common(7), 1, 1); axC.axhline(0, color="#333", lw=0.7); ybot = draw(mr.most_common(7), -1, -1)
axC.set_xlim(z0-1, z1+8); axC.set_ylim(ybot-1.6, ytop+0.6)
for sp in ("top","left","right"): axC.spines[sp].set_visible(False)
axC.set_yticks([]); axC.spines["bottom"].set_position(("data", ybot-1.0))
tk = np.linspace(z0, z1, 5).astype(int); axC.set_xticks(tk); axC.set_xticklabels([f"{t:,}" for t in tk], fontsize=6.5); axC.tick_params(axis="x", length=3)
axC.set_xlabel(f"{STRAIN.replace('_','/')} chr{CH} position (bp)", fontsize=7)
axC.text(z0-0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0-0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
for x2 in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=((x2-ps)/N, -1.30), coordsA=axB.transData, xyB=(x2, ytop+0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0,(3,2))))
pc.pbadge(axC, "C", f"Base resolution, {STRAIN.replace('_','/')} at {TPLAB[CHOSEN]} · 5′ arrow RED = antisense-to-TE (silencing), grey = sense-to-TE", fs=7.6)
fig.suptitle(TITLE, fontsize=11.5, fontweight="bold", y=0.965)
# --- per-figure SourceData: panel-A (per-strain×tp FPM, minimap2 qcov, aligned coordinate) + summary ---
_sd = []
for s in ORDER:
    for tp in TPS:
        _sd.append(dict(panel="A_perstrain", strain=s, timepoint=tp, fpm_allread=round(float(FPM.loc[s, tp]), 3), minimap2_qcov=round(float(QCOV.get(s, 0)), 3), aligned_chrom=(COORD[s][0] if COORD.get(s) else ""), aligned_pos=(COORD[s][1] if COORD.get(s) else "")))
_sd.append(dict(panel="A_perstrain", strain="GRCm39", timepoint="(all)", fpm_allread=0.0, minimap2_qcov=round(float(QCOV.get("GRCm39", 0)), 3), aligned_chrom=(COORD["GRCm39"][0] if COORD.get("GRCm39") else ""), aligned_pos=(COORD["GRCm39"][1] if COORD.get("GRCm39") else "")))
for _k, _v in [("carrier", STRAIN), ("locus", f"chr{CH}:{ps}-{pe}"), ("dominant_TE", domTE), ("minimap2_present_of17", npres), ("carrier_primary_piRNAs", int(_ntot))]:
    _sd.append(dict(panel="summary", strain=str(_k), timepoint="", fpm_allread="", minimap2_qcov="", aligned_chrom="", aligned_pos=str(_v)))
import os as _os; _os.makedirs(f"{T22}/data/source_data", exist_ok=True); pd.DataFrame(_sd).to_csv(f"{T22}/data/source_data/SourceData_{OUT}.csv", index=False)
for e in ("pdf","svg","png"): fig.savefig(f"{T22}/figures/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png + SourceData_{OUT}.csv")
