#!/usr/bin/env python3
"""THEME 22 — NON-REFERENCE locus figure (make_pav_locus style), with the pangenome GRAPH (odgi inject + pav) as the
cross-strain evidence. Built for a piRNA cluster GENUINELY ABSENT from GRCm39.
  (A) Pangenome x timepoint — per-strain PICB-cluster FPM (log; via co-location to each strain) + GENOME-PRESENCE
      markers from the GRAPH (odgi coverage): ● full / ◐ partial-segmental / ○ absent (genetic loss). The graph gives
      QUANTITATIVE coverage, so partial/segmental presence is shown (liftover can only do binary). GRCm39 = ○.
      Tally: present in N/16 · partial in P · genetically absent in M · present-but-silent (regulatory) in K.
  (B) the carrier strain's per-timepoint sRNA coverage + TE + gene tracks; (C) base resolution.
Panels B/C reuse the strain-anchored pav_clusters blocks. Usage: make_nonref_locus.py <strain> <chrom> <start> <end> <graph_pav_name> <title> <outbase>"""
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
_fpmcache = {}
def cl_fpm(X):
    if X not in _fpmcache: _fpmcache[X] = pd.read_csv(f"{CP}/{X}.clusters_fpm.bed", sep="\t", header=None, names=["c", "s", "e", "allF", "uniqF", "st", "tp"], dtype={"c": str})
    return _fpmcache[X]
def fpm_at(X, c, s, e):   # per-tp all-read FPM of X's cluster(s) overlapping c:s-e
    f = cl_fpm(X); m = f[(f.c == str(c)) & (f.s < e) & (f.e > s)]; return {tp: float(m[m.tp == tp].allF.sum()) for tp in TPS}
# --- per-strain locus coords (co-location), FPM, and +/- STRAND proportion (sRNA) ---
COORDS = {STRAIN: (str(CH), ps, pe)}
for Y in ORDER:
    if Y == STRAIN: continue
    f = f"{D}/colo/{STRAIN}__{Y}.lifted.bed"
    if not os.path.exists(f): continue
    regs = [(c[0], int(c[1]), int(c[2])) for c in (l.rstrip("\n").split("\t") for l in open(f)) if c[3] == GPNAME]
    if not regs: continue
    dom = Counter(r[0] for r in regs).most_common(1)[0][0]; rr = [r for r in regs if r[0] == dom]
    COORDS[Y] = (dom, min(r[1] for r in rr), max(r[2] for r in rr))
DAT = {tp: (d if (d := fetch_primary(STRAIN, CH, ps, pe, tp, nb)) and d["ntot"] > 0 else None) for tp in TPS}
FPM = pd.DataFrame(0.0, index=ORDER, columns=TPS); PF = {s: {} for s in ORDER}   # PF[strain][tp] = + strand fraction
for Y, (yc, ys, ye) in COORDS.items():
    for tp, v in fpm_at(Y, yc, ys, ye).items(): FPM.loc[Y, tp] = v
    for tp in TPS:
        d = DAT[tp] if Y == STRAIN else fetch_primary(Y, yc, ys, ye, tp, nb)
        PF[Y][tp] = (d["nplus"] / d["ntot"]) if (d and d["ntot"]) else 0.5
gp = pd.read_csv(f"{D}/graph_check_pav3.tsv", sep="\t"); grow = gp[gp.name == GPNAME]
cov = {s: (float(grow.iloc[0].get(s, 0.0)) if len(grow) else 0.0) for s in ORDER}; gr_cov = float(grow.iloc[0].get("GRCm39", 0.0)) if len(grow) else 0.0
def state(c): return "present" if c >= 0.8 else ("partial" if c >= 0.2 else "absent")
PAV = {s: state(cov[s]) for s in ORDER}
npres = sum(PAV[s] == "present" for s in ORDER); npart = sum(PAV[s] == "partial" for s in ORDER); nabs = sum(PAV[s] == "absent" for s in ORDER)
nsil = sum(1 for s in ORDER if PAV[s] != "absent" and FPM.loc[s].max() == 0)
tes = next((DAT[tp]["tes"] for tp in TPS if DAT[tp]), te_at(STRAIN, CH, ps, pe)); domTE = dom_te_family(tes, ps, pe)
CHOSEN = max((tp for tp in TPS if DAT[tp]), key=lambda tp: DAT[tp]["ntot"], default=TPS[0])
print(f"{TITLE}: {STRAIN} chr{CH}:{ps:,}-{pe:,} domTE={domTE} | graph: present {npres} partial {npart} absent {nabs} silent {nsil}; GRCm39={gr_cov:.2f}")
def testr(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(12.6, 12.4), dpi=300); gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 0.36, 2.2, 1.35], hspace=0.55); fig.subplots_adjust(top=0.9, bottom=0.06)
axA = fig.add_subplot(gs[0]); axAc = fig.add_subplot(gs[1], sharex=axA); axB = fig.add_subplot(gs[2]); axC = fig.add_subplot(gs[3])
# === A: Pangenome x timepoint FPM (log), each bar split by HEIGHT: solid = + strand / pale = - strand ===
import matplotlib.cm as cm
x = np.arange(len(ORDER)); bw = 0.27; gx = len(ORDER) + 0.4
for j, tp in enumerate(TPS):
    h = FPM[tp].values
    axA.bar(x + (j - 1) * bw, np.maximum(h, 1e-3), width=bw, color=TPCOL[tp], edgecolor="white", linewidth=0.2, label=TPLAB[tp])   # solid = + strand base
    for xi in range(len(ORDER)):
        if h[xi] > 0:
            sp = max(h[xi] * PF[ORDER[xi]].get(tp, 0.5), 0.1)
            if h[xi] - sp > 1e-3: axA.bar(xi + (j - 1) * bw, h[xi] - sp, bottom=sp, width=bw, color=pc.pale(TPCOL[tp], 0.5), edgecolor="none", zorder=3)   # pale top = - strand portion
            axA.text(xi + (j - 1) * bw, h[xi] * 1.25, (f"{h[xi]:.0f}" if h[xi] >= 1 else f"{h[xi]:.1f}"), ha="center", va="bottom", fontsize=4.2, rotation=90, color=TPCOL[tp], fontweight="bold")
axA.bar([gx], [1e-3], width=bw, color="#ddd")
axA.set_yscale("log"); axA.set_ylim(0.1, max(FPM.values.max(), 1) * 24); axA.set_ylabel("PICB cluster\nFPM (log)", fontsize=8)
axA.spines[["top", "right"]].set_visible(False); plt.setp(axA.get_xticklabels(), visible=False)
lh = [Patch(facecolor=TPCOL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#555", label="solid = + strand"), Patch(facecolor=pc.pale("#555", 0.5), label="pale = − strand")]
axA.legend(handles=lh, fontsize=5.9, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.004, 1.02), handlelength=1.1, handletextpad=0.4, title="timepoint · strand", title_fontsize=6.2)
pc.pbadge(axA, "A", "Pangenome × timepoint — PICB-cluster FPM across 16 strains (log); each bar split by HEIGHT: solid = + strand / pale = − strand", fs=8.0, y=1.18)
pc.rtext(axA, 0.5, 0.945, [(f"present in {npres}/16 strains", pc.C_DUAL, True), ("·", "#bbb", False), (f"partial/segmental in {npart}", "#3182bd", bool(npart)), ("·", "#bbb", False), (f"genetically absent (incl. GRCm39) in {nabs+1}", pc.C_SILENCE, True), ("·", "#bbb", False), (f"present-but-silent in {nsil}", "#5a6b7a", bool(nsil))], fs=7.4, transform=fig.transFigure, center=True)
# === A-cov: GRAPH coverage strip (odgi pav) — HOW present/partial each strain is (continuous; liftover = binary) ===
def covcol(c): return "#C0392B" if c < 0.2 else cm.Blues(0.35 + 0.6 * min(c, 1.0))
for xi, X in enumerate(ORDER):
    axAc.bar(xi, max(cov[X], 0.012), width=0.74, color=covcol(cov[X]), edgecolor="white", linewidth=0.3)
    if 0.2 <= cov[X] < 0.92: axAc.text(xi, cov[X] + 0.05, f"{cov[X]*100:.0f}%", ha="center", va="bottom", fontsize=4.6, color="#2166ac", fontweight="bold")
axAc.bar(gx, 0.012, width=0.74, color="#C0392B", edgecolor="white"); axAc.text(gx, 0.07, "0", ha="center", va="bottom", fontsize=5.2, color="#C0392B", fontweight="bold")
axAc.axhline(0.5, ls=(0, (2, 2)), color="#bbb", lw=0.6); axAc.set_ylim(0, 1.12); axAc.set_yticks([0, 1]); axAc.tick_params(labelsize=6); axAc.set_ylabel("graph\ncov.", fontsize=6.6); axAc.spines[["top", "right"]].set_visible(False)
axAc.set_xticks(list(x) + [gx]); axAc.set_xticklabels([s.replace("_", "/") for s in ORDER] + ["GRCm39"], rotation=90, fontsize=7.0)
for t, X in zip(axAc.get_xticklabels()[:len(ORDER)], ORDER):
    if X in WILD: t.set_color("#C0392B")
axAc.get_xticklabels()[-1].set_color("#C0392B"); axAc.get_xticklabels()[-1].set_fontweight("bold")
axAc.text(0.0, 1.34, "GRAPH coverage (odgi pav) — HOW present each strain is: full ≈ 1.0 · partial/segmental (% labelled) · absent (genetic loss); GRCm39 = 0", transform=axAc.transAxes, fontsize=6.4, color="#555", fontweight="bold", ha="left", va="bottom")
# === B: carrier sRNA coverage + TE + gene (reuse draw_strain_block) ===
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()
_dc = DAT[CHOSEN]; reads = _dc["reads"] if _dc else []
_f5 = Counter(r[0] for r in reads if not r[2]); _r5 = Counter(r[1] - 1 for r in reads if r[2]); _rn = lambda i: sum(_r5[j] for j in range(i - 2, i + 13))
_cd = [(min(_f5[i], _rn(i)), _f5[i] + _rn(i), i + 5) for i in _f5 if _rn(i) > 0]
pk = max(_cd)[2] if _cd else ((Counter((r[1] - 1 if r[2] else r[0]) for r in reads).most_common(1)[0][0]) if reads else ps + N // 2); z0, z1 = pk - 30, pk + 50
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, DAT, tes, CH, ps, pe, FPM.loc[STRAIN], 0.0, z0, z1, name=STRAIN, is_top=True, wild=(STRAIN in WILD), TECOL=TECOL, only_tp=CHOSEN)
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
for x2 in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=((x2 - ps) / N, -1.30), coordsA=axB.transData, xyB=(x2, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
pc.pbadge(axC, "C", f"Base resolution, {STRAIN.replace('_','/')} at {TPLAB[CHOSEN]} (top tp)   ·   5′ arrow RED = antisense-to-TE (silencing), grey = sense-to-TE", fs=7.6)
fig.suptitle(TITLE, fontsize=11.5, fontweight="bold", y=0.978)
_OD = "24_strain_private_te_loci" if OUT.startswith("Fig_private_locus") else "23_nonref_locus_atlas"
# --- per-figure SourceData: panel-A quantitative content (per-strain×tp FPM, +strand fraction, graph coverage, PAV) + summary ---
_sd = []
for s in ORDER:
    for tp in TPS:
        _sd.append(dict(panel="A_perstrain", strain=s, timepoint=tp, fpm_allread=round(float(FPM.loc[s, tp]), 3), plus_strand_frac=round(float(PF[s].get(tp, 0.5)), 3), graph_coverage=round(float(cov[s]), 4), pav_state=PAV[s]))
_sd.append(dict(panel="A_perstrain", strain="GRCm39", timepoint="(all)", fpm_allread=0.0, plus_strand_frac=np.nan, graph_coverage=round(gr_cov, 4), pav_state=state(gr_cov)))
for _k, _v in [("carrier", STRAIN), ("locus", f"chr{CH}:{ps}-{pe}"), ("dominant_TE", domTE), ("present_in_of16", npres), ("partial_in", npart), ("genetically_absent_incl_GRCm39", nabs + 1), ("present_but_silent_in", nsil), ("carrier_primary_piRNAs", int(_ntot))]:
    _sd.append(dict(panel="summary", strain=str(_k), timepoint="", fpm_allread="", plus_strand_frac="", graph_coverage="", pav_state=str(_v)))
_dd = f"{ROOT}/figures/analysis_figures/{_OD}/data/source_data"; os.makedirs(_dd, exist_ok=True)
pd.DataFrame(_sd).to_csv(f"{_dd}/SourceData_{OUT}.csv", index=False)
for e in ("pdf", "svg", "png"): fig.savefig(f"{ROOT}/figures/analysis_figures/{_OD}/figures/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png + SourceData_{OUT}.csv")
