#!/usr/bin/env python3
"""Strain-private piRNA SOURCE LOCUS (individual piRNAs, NOT a PICB cluster) in the MULTI pangenome layout —
matches make_pav_locus_multi / Fig_pav_locus_chr1wildtrio_multi:
  (A) Pangenome × timepoint — source-locus expression (RPM, log) across the 16 strains + 3-state genome PAV
      (● locus present in genome / ○ absent = genetic loss; present-but-no-bar = present-but-silent) →
      genetic-loss vs regulatory-silencing 'why' (insertion-driven vs divergence-driven).
  (B) every EXPRESSING strain: per-timepoint PRIMARY-read coverage (deep ↑+ / pale ↓−) + TE + gene tracks.
  (C) base resolution in the top strain.
The locus is projected from the carrier into every strain via the minigraph-cactus HAL (presence; whole-genome
alignment, NOT chain-file liftover); each strain's strain-wise BAM gives expression. Own coordinates throughout.
Usage: make_source_pav_multi.py <carrier> <CHROM_own e.g. CAST_EiJ#1#chr10> <start> <end> <te_label> <te_strand> <outbase>"""
import sys, os, subprocess, tempfile, numpy as np, pysam
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, fetch_primary, te_at, dom_te_family
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"; RB = f"{ROOT}/results/STAR_srna_strain_wise"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}; comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb"}
CARRIER, CHROM, S, E, TELAB, TEST, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7]
PROJ = sys.argv[8] if len(sys.argv) > 8 else None   # optional precomputed projection (build_source_projection.py) — skips per-locus HAL liftover
CHc = CHROM.split("chr")[-1]; nb = 180; EXPR_MIN = 50   # >= EXPR_MIN pooled primary reads => expressed; present(maps)+below => present-but-silent
# ---- (1) genome presence + own coords per strain: precomputed projection if given, else live HAL ----
if PROJ and os.path.exists(PROJ):
    _lid = f"{CARRIER}:{CHc}:{S}"; _pt = pd.read_csv(PROJ, sep="\t")
    _pt = _pt[(_pt.locus_id == _lid) & (_pt.present.astype(bool))]
    owncoord = {CARRIER: (CHc, S, E)}
    for _, _r in _pt.iterrows():
        if _r.target != CARRIER: owncoord[_r.target] = (str(_r.tgt_chrom), int(_r.tgt_start), int(_r.tgt_end))
else:
    bed = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False); bed.write(f"{CHc}\t{S}\t{E}\tloc\n"); bed.close()
    owncoord = {CARRIER: (CHc, S, E)}
    for X in ORDER:
        if X == CARRIER: continue
        out = subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} {CARRIER} {bed.name} {X} /dev/stdout", shell=True, capture_output=True, text=True).stdout
        frs = [(ff[0], int(ff[1]), int(ff[2])) for ff in (ln.split("\t") for ln in out.splitlines()) if len(ff) >= 3]
        if frs:
            ch = Counter(f[0] for f in frs).most_common(1)[0][0]; fr = [f for f in frs if f[0] == ch]
            s0 = min(f[1] for f in fr); e0 = max(f[2] for f in fr); mapped = sum(f[2] - f[1] for f in fr)
            if mapped >= 0.30 * (E - S): owncoord[X] = (ch, s0, e0)   # >=30% of the window maps => present
    os.unlink(bed.name)
# ---- (1b) SNAP the carrier window to its actual PICB cluster near the breakpoint ----
# the master/projection coords are the one-directional TE-BREAKPOINT window (5 kb); the TE-created cluster frequently
# sits on the other side, so that window misses it (9/27 creation loci captured <60% of reads, 6 of them <25%).
# Snap the carrier window to its real PICB cluster(s) within ±15 kb of the breakpoint (own coords), padded ±3 kb.
_P = pc._pang(); _cc = _P[(_P.strain == CARRIER) & (_P.own_chrom.astype(str) == CHc)].copy()
_cc["_os"] = pd.to_numeric(_cc.own_start, errors="coerce"); _cc["_oe"] = pd.to_numeric(_cc.own_end, errors="coerce")
_cc = _cc.dropna(subset=["_os", "_oe"]); _cc = _cc[(_cc._oe > S - 15000) & (_cc._os < E + 15000)].sort_values("_os")
if len(_cc):                                                                       # gap-MERGE nearby clusters into contiguous groups (avoid over-broadening from scattered paralogs), pick the group overlapping the breakpoint (else densest by FPM)
    _grp = []; _cs = _ce = _cf = None
    for _, _r in _cc.iterrows():
        if _cs is None: _cs, _ce, _cf = _r._os, _r._oe, _r.all_primary_FPM
        elif _r._os <= _ce + 8000: _ce = max(_ce, _r._oe); _cf += _r.all_primary_FPM
        else: _grp.append((_cs, _ce, _cf)); _cs, _ce, _cf = _r._os, _r._oe, _r.all_primary_FPM
    _grp.append((_cs, _ce, _cf))
    _best = max(_grp, key=lambda g: (max(0, min(E, g[1]) - max(S, g[0])) > 0, g[2]))
    _gs, _ge = int(_best[0]) - 3000, int(_best[1]) + 3000
    if _ge - _gs > 45000: _gs, _ge = S - 12000, E + 12000   # giant paralog chain (e.g. LP_J chr5 spans 223 kb of clusters <8kb apart) → fall back to a bounded window at the breakpoint
    owncoord[CARRIER] = (CHc, max(1, _gs), _ge)
# ---- (2) per (present) strain: PRIMARY-read coverage per timepoint + RPM; classify expressed / silent ----
def libsize(X, tp):
    t = 0
    for r in (1, 2, 3):
        b = f"{RB}/{X}/{X}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
        if os.path.exists(b):
            try: t += sum(s.mapped for s in pysam.AlignmentFile(b, "rb").get_index_statistics())
            except Exception: pass
    return t
def collect(X, oc, ps, pe):
    per = {}; plus = np.zeros(nb); minus = np.zeros(nb); reads = []; nat = nte = n1u = npl = nmi = 0; tes = None; fps = {}
    for tp in TPS:
        d = fetch_primary(X, oc, ps, pe, tp, nb)
        if d and d["ntot"] > 0:
            per[tp] = d; plus += d["plus"]; minus += d["minus"]; reads += d["reads"]
            nat += d["nat"]; nte += d["nte"]; n1u += d["n1u"]; npl += d["nplus"]; nmi += d["nminus"]; tes = d["tes"]
            tm = libsize(X, tp); fps[tp] = (d["ntot"] / tm * 1e6) if tm else 0.0
        else: per[tp] = None
    if tes is None: tes = te_at(X, oc, ps, pe)
    ntot = npl + nmi
    fam_all = {}; fam_as = {}                                       # piRNA reads (+ antisense-to-TE subset) per TE family, pooled across tp (for 'TE-piRNA' + 'AS→TE' header tops)
    for _tp in TPS:
        if per[_tp]:
            for _f, _v in per[_tp]["fam_bd"].items(): fam_all[_f] = fam_all.get(_f, 0) + _v[0]; fam_as[_f] = fam_as.get(_f, 0) + _v[1]
    return dict(oc=oc, ps=ps, pe=pe, fps=fps, per=per, plus=plus, minus=minus, reads=reads, tes=tes, ntot=ntot, fam_as=fam_as,
                nplus=npl, nminus=nmi, n1u=n1u, nat=nat, nte=nte, domTE=dom_te_family(tes, ps, pe), fam_all=fam_all,
                pct_at=100 * nat / max(1, nte), arch="dual-strand" if min(nmi, ntot - nmi) / max(1, ntot) > 0.2 else "uni-strand")
COVall = {X: collect(X, *owncoord[X]) for X in owncoord}
present = [X for X in ORDER if X in COVall and COVall[X]["ntot"] >= EXPR_MIN]   # EXPRESSING strains (canonical order)
COV = {X: COVall[X] for X in present}
TOP = max(present, key=lambda X: COV[X]["ntot"]) if present else CARRIER
FPM = pd.DataFrame(0.0, index=ORDER, columns=TPS)
for X in present:
    for tp in TPS: FPM.loc[X, tp] = COV[X]["fps"].get(tp, 0.0)
PAV = {X: ("present" if X in owncoord else "absent") for X in ORDER}   # genome presence (mapped via HAL)
nabs = sum(1 for X in ORDER if PAV[X] == "absent"); nsil = sum(1 for X in ORDER if PAV[X] == "present" and X not in present)
print(f"{OUT}: carrier={CARRIER}; expressing={present}; top={TOP}; genome-present={sum(1 for X in ORDER if PAV[X]=='present')}; silent={nsil}; absent={nabs}")
# ---- (3) figure (layout reused from make_pav_locus_multi) ----
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
nP = max(1, len(present)); fig = plt.figure(figsize=(14, 9.8 + 1.7 * nP), dpi=300)
gs = fig.add_gridspec(3, 1, height_ratios=[1.9, 1.3 * nP + 1.0, 2.6 + 0.18 * nP], hspace=0.72); fig.subplots_adjust(top=0.88, bottom=0.07)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
def pbadge(ax, L, title, fs=8.2, y=1.05, tx=0.042, bx=0.012):
    ax.text(bx, y, f"{L}", transform=ax.transAxes, fontsize=9.5, fontweight="bold", color="white", ha="center", va="center", bbox=dict(boxstyle="circle,pad=0.30", fc="#2C3E50", ec="none"), zorder=10)
    ax.text(tx, y, title, transform=ax.transAxes, fontsize=fs, fontweight="bold", color="#1a1a1a", ha="left", va="center")
x = np.arange(len(ORDER)); bw = 0.26
for j, tp in enumerate(TPS):
    h = FPM[tp].values; axA.bar(x + (j - 1) * bw, np.maximum(h, 1e-3), width=bw, color=TPCOL[tp], edgecolor="white", linewidth=0.2, label=TPLAB[tp])
    for xi, X in enumerate(ORDER):       # STRAND: split each bar by HEIGHT → solid + (bottom) / pale − (top)
        if h[xi] > 0 and X in COV and COV[X]["per"][tp] is not None and COV[X]["per"][tp]["ntot"] > 0:
            pc.panelA_strand(axA, xi + (j - 1) * bw, h[xi], bw, COV[X]["per"][tp]["nplus"] / COV[X]["per"][tp]["ntot"], TPCOL[tp])
        if h[xi] > 0: axA.text(xi + (j - 1) * bw, h[xi] * 1.18, f"{h[xi]:.0f}" if h[xi] >= 1 else f"{h[xi]:.1f}", ha="center", va="bottom", fontsize=4.8, rotation=90, color=TPCOL[tp], fontweight="bold")
axA.set_yscale("log"); axA.set_ylim(0.1, max(FPM.values.max(), 1) * 20); axA.set_xticks(x); axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=7.5)
for t, X in zip(axA.get_xticklabels(), ORDER):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
for xi, X in enumerate(ORDER):
    pv = PAV[X]; fc = ("#333" if pv == "present" else "white"); ec = ("#333" if pv == "present" else "#C0392B")
    axA.plot(xi, 0.135, marker="o", markersize=4.2, markerfacecolor=fc, markeredgecolor=ec, markeredgewidth=0.9, clip_on=False)
axA.set_ylabel("source-locus\nRPM (log)", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
axA.legend([TPLAB[t] for t in TPS], fontsize=7.5, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.005, 1.0), title="timepoint", title_fontsize=7, handlelength=1.2)
axA.plot([], [], "o", markersize=4.2, markerfacecolor="#333", markeredgecolor="#333", label="locus in genome")
axA.plot([], [], "o", markersize=4.2, markerfacecolor="white", markeredgecolor="#C0392B", label="locus absent (genetic loss)")
axA.legend(fontsize=6.0, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.005, 0.42), handlelength=1.0, handletextpad=0.3)
pbadge(axA, "A", "Pangenome × timepoint — strain-private piRNA SOURCE LOCUS expression across 16 strains (RPM, log)   ·   ● present in genome   ○ absent", fs=8.0, y=1.16)
axA.text(0.012, 1.05, "each RPM bar split by HEIGHT (proportion):  solid = + strand  ·  pale = − strand", transform=axA.transAxes, fontsize=6.2, color="#777", style="italic", ha="left", va="center")
dT = COV[TOP]; psT, NT_ = dT["ps"], max(1, dT["pe"] - dT["ps"])
_f5 = Counter(r[0] for r in dT["reads"] if not r[2]); _r5 = Counter(r[1] - 1 for r in dT["reads"] if r[2]); _rn = lambda i: sum(_r5[j] for j in range(i - 2, i + 13))
_cd = [(min(_f5[i], _rn(i)), _f5[i] + _rn(i), i + 5) for i in _f5 if _rn(i) > 0]
pk = max(_cd)[2] if _cd else (Counter((r[1] - 1 if r[2] else r[0]) for r in dT["reads"]).most_common(1)[0][0] if dT["reads"] else psT + NT_ // 2); z0, z1 = pk - 30, pk + 50
ztps = [tp for tp in TPS if dT["per"][tp] is not None and any(r[0] < z1 and r[1] > z0 for r in dT["per"][tp]["reads"])]; ztpL = " + ".join(TPLAB[t] for t in ztps) or "—"
SP = 4.8; off_top = -(nP - 1) * SP
axB.set_xlim(0, 1); axB.set_ylim(off_top - 3.0, 2.1); axB.axis("off"); fig.canvas.draw()
for i, X in enumerate(present):
    d = COV[X]; off = -i * SP; N = max(1, d["pe"] - d["ps"]); xr = np.linspace(0, 1, nb)
    pm = max([max(d["per"][t]["plus"].max(), d["per"][t]["minus"].max()) for t in TPS if d["per"][t] is not None] + [1.0])
    BL = {"16.5dpc": off + 0.85, "12.5dpp": off + 0.0, "20.5dpp": off - 0.85}; HS = 0.38
    if X == TOP:
        zr0 = (z0 - d["ps"]) / N; zr1 = (z1 - d["ps"]) / N
        axB.add_patch(Rectangle((zr0, off - 1.30), zr1 - zr0, 2.60, facecolor="#FDB863", alpha=0.30, edgecolor="#E8A33D", lw=0.5, zorder=0.5))
    for _tp in TPS:
        _b = BL[_tp]; _pt = d["per"][_tp]
        axB.axhline(_b, color="#dddddd", lw=0.3, zorder=1)
        axB.text(-0.015, _b, TPLAB[_tp], fontsize=4.6, ha="right", va="center", color=(TPCOL[_tp] if _pt is not None else "#ccc"), fontweight="bold")
        if _pt is None:
            axB.text(0.5, _b, "OFF — not expressed (this tp)", fontsize=4.2, ha="center", va="center", color="#cccccc", style="italic"); continue
        axB.fill_between(xr, _b, _b + _pt["plus"] / pm * HS, color=pc.PLUS_COL[_tp], alpha=0.92, step="mid", lw=0)
        axB.fill_between(xr, _b, _b - _pt["minus"] / pm * HS, color=pc.MINUS_COL[_tp], alpha=0.92, step="mid", lw=0)
    axB.add_patch(Rectangle((0, off - 1.74), 1.0, 0.42, facecolor="#f6f6f6", edgecolor="#e3e3e3", lw=0.4, zorder=0))
    for (ts, te, st, fam) in d["tes"]:
        x0 = max(0.0, (ts - d["ps"]) / N); x1 = min(1.0, (te - d["ps"]) / N)
        if x1 > x0:
            fam2 = fam.split("|")[-1] if "|" in fam else fam
            axB.add_patch(Rectangle((x0, off - 1.70), x1 - x0, 0.34, facecolor=TECOL.get(fam2, "#dddddd"), edgecolor="white", lw=0.25, zorder=2))
            if x1 - x0 > 0.014:
                af, at = (x0, x1) if st == "+" else (x1, x0)
                axB.annotate("", xy=(at, off - 1.53), xytext=(af, off - 1.53), arrowprops=dict(arrowstyle="-|>", color="#111", lw=0.45, mutation_scale=5.0), zorder=3)
    axB.text(-0.075, off - 1.53, "TEs", fontsize=6.6, ha="right", va="center", color="#333", fontweight="bold", bbox=dict(boxstyle="round,pad=0.2", fc="#ececec", ec="#c9c9c9", lw=0.5))
    axB.add_patch(Rectangle((0, off - 2.24), 1.0, 0.40, facecolor="#edf2f9", edgecolor="#d8e0ec", lw=0.4, zorder=0))
    _genes = list(pc.genes_at(X, d["oc"], d["ps"], d["pe"]))
    _gw = [min(1.0, (ge - d["ps"]) / N) - max(0.0, (gs - d["ps"]) / N) for (gs, ge, gst, gn, gb) in _genes]
    _sym = lambda nm: not str(nm).startswith(("ENS", "MGP_"))
    _cand = [i for i in range(len(_genes)) if _sym(_genes[i][3]) or _gw[i] >= 0.04]
    _lbl = set(sorted(_cand, key=lambda i: (-int(_sym(_genes[i][3])), -_gw[i]))[:5]); _gi = 0
    for i, (gs, ge, gst, gname, gbt) in enumerate(_genes):
        gx0 = max(0.0, (gs - d["ps"]) / N); gx1 = min(1.0, (ge - d["ps"]) / N)
        if gx1 <= gx0: continue
        axB.add_patch(Rectangle((gx0, off - 2.20), gx1 - gx0, 0.32, facecolor="#c9d6ea", edgecolor=pc.C_GENE, lw=0.8, zorder=2))
        if gx1 - gx0 > 0.014:
            gaf, gat = (gx0, gx1) if gst == "+" else (gx1, gx0)
            axB.annotate("", xy=(gat, off - 2.04), xytext=(gaf, off - 2.04), arrowprops=dict(arrowstyle="-|>", color=pc.C_GENE, lw=0.6, mutation_scale=5.0), zorder=3)
        if i in _lbl:
            axB.text((gx0 + gx1) / 2, off - 2.40 - 0.17 * (_gi % 3), gname, fontsize=4.6, ha="center", va="top", color=pc.C_GENE, style="italic", fontweight="bold"); _gi += 1
    axB.text(-0.075, off - 2.04, "genes", fontsize=6.6, ha="right", va="center", color="white", fontweight="bold", bbox=dict(boxstyle="round,pad=0.2", fc=pc.C_GENE, ec=pc.C_GENE, lw=0.5))
    axB.text(-0.015, off + 1.80, X.replace("_", "/") + ("  ▼ zoom" if X == TOP else ""), fontsize=7.0, ha="right", va="center", fontweight="bold", color=pc.C_WILD if X in WILD else "#1a1a1a", bbox=(dict(boxstyle="round,pad=0.16", fc="#FDF3E3", ec=pc.C_ZOOM, lw=0.7) if X == TOP else None))
    axB.text(-0.015, off + 1.40, f"chr{d['oc']}:{d['ps']:,}", fontsize=4.6, ha="right", va="center", color=pc.C_META, style="italic")
    _fam = d["domTE"] or "—"
    pc.rtext(axB, 1.015, off + 1.80, [(f"{d['ntot']:,} primary reads", "#444", True), ("·  TE-DNA (genome)", "#999", False), (f"{_fam}", pc.famcol_text(_fam), True)], fs=6.0)   # strain-level genomic only; TE-piRNA + AS→TE are PER-TIMEPOINT (in each row)
    tx0 = 1.035; TW = 0.26                                    # right margin (snug after coverage): ALL-locus-piRNA bar = TE families ranked by coverage + non-TE, each split by STRAND (solid + / pale −)
    for _tp in TPS:
        _y = BL[_tp]; pt = d["per"][_tp]
        if pt is None or pt["ntot"] == 0: continue
        _p = pt["nplus"]; _m = pt["nminus"]; _t = _p + _m
        pc.te_strand_bar(axB, tx0, TW, _y, 0.16, pt["fam_bd"], pt["ntot"], pt["nte"])
        _pct_te = 100 * pt["nte"] / pt["ntot"]; _as_te = 100 * sum(v[1] for v in pt["fam_bd"].values()) / max(1, pt["nte"]); _dual = min(_p, _m) / max(1, _t) > 0.2; _1u_tp = 100 * pt["n1u"] / max(1, _t); _rpm = d["fps"].get(_tp, 0); _u1 = pt["u1"]; _u1pk = ((np.argmax(_u1) + 0.5) / len(_u1) * (d["pe"] - d["ps"]) / 1000) if _u1.sum() > 0 else None   # 1U-peak offset (kb) from cluster start = dominant primary-processing site
        _topf = max(pt["fam_bd"], key=lambda f: pt["fam_bd"][f][0]) if pt["fam_bd"] else "—"                                                  # per-tp top TE-piRNA family
        _asf_tp = max(pt["fam_bd"], key=lambda f: pt["fam_bd"][f][1]) if (pt["fam_bd"] and any(v[1] for v in pt["fam_bd"].values())) else "—"   # per-tp top antisense-to-TE (silenced) family
        _inc = (X == TOP) and any(r[0] < z1 and r[1] > z0 for r in pt["reads"])
        _tc = pc.TPCOL[_tp]   # timepoint-coloured stats row (family + AS→TE keep semantic colours)
        _segs = [(f"{_pct_te:.0f}% on TE", _tc, _pct_te >= 40), ("piRNA", _tc, False), (f"{_topf}", pc.famcol_text(_topf), True), (f"{_as_te:.0f}% AS→", "#C0392B" if _as_te >= 40 else pc.pale("#C0392B", 0.35), True), (f"{_asf_tp}", pc.famcol_text(_asf_tp), True), ("↑↓DUAL" if _dual else "→UNI", _tc, _dual), (f"RPM {_rpm:.1f}", _tc, True), (f"{_1u_tp:.0f}%1U", _tc, False)] + ([(f"1U▲+{_u1pk:.1f}kb", _tc, True)] if _u1pk is not None else []) + ([("→C", pc.C_ZOOM, True)] if _inc else [])
        pc.rtext(axB, tx0, _y - 0.17, _segs, fs=4.1, sep=" · ", va="top")   # per-tp stats UNDER the bar (left-aligned at bar start), with a clear gap below the bar → no overlap, no right-edge cut-off
axB.set_xlim(0, 1); axB.set_ylim(off_top - 3.0, 2.1); axB.axis("off")
pbadge(axB, "B", "Per-timepoint sRNA coverage at the source locus in each EXPRESSING strain — height = expression, colour = timepoint (deep ↑+ / pale ↓−)   ·   TE + gene tracks below each strain   ·   example → zoom C", fs=7.3, y=1.07)
axB.text(0.012, 1.02, "‘primary reads’ = each sRNA read (24–32 nt) counted once at its STAR primary locus (multimappers kept, not double-counted) · architecture (genomic strand) ≠ sense/antisense (relative to TE)", transform=axB.transAxes, fontsize=5.2, color="#8a8a8a", style="italic", ha="left", va="center")
famset = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for X in present for (_, _, _, f) in COV[X]["tes"]))[:6]
_tpkey = [Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS]
_stkey = [Patch(facecolor="#6a3d9a", label="solid = + strand"), Patch(facecolor=pc.pale("#6a3d9a", 0.55), label="pale = − strand"), Patch(facecolor="#efefef", label="non-TE piRNA"), Patch(facecolor="#C0392B", label="antisense-to-TE = silencing (red outline); bars = TE family"), Patch(facecolor="#cfcfcf", label="sense-to-TE"), Patch(facecolor="#e9edf3", edgecolor=pc.C_GENE, label="gene model (GFF)")]
_legH = _tpkey + _stkey + [Patch(facecolor=pc.famcol(f), label=f) for f in famset]
pc.color_legend(axB.legend(handles=_legH, fontsize=6.4, frameon=False, loc="upper center", ncol=9, bbox_to_anchor=(0.5, -0.015), title="per-tp right-margin element (all locus piRNA):  ▲ TOP = TE families RANKED by coverage (colours below), each split by STRAND (solid + / pale −) + non-TE; tick = TE|non-TE   ·   ▼ BOTTOM = antisense-to-TE broken by TE family (ranked; RED underline = silencing span) + sense-to-TE (grey)   ·   Panel A bars split ±strand by height   ·   timepoint colour + gene models", title_fontsize=6.3), _legH)   # legend labels coloured by swatch
d = COV[TOP]; ps = d["ps"]; reads = d["reads"]; N = max(1, d["pe"] - d["ps"]); CH = d["oc"]
def tst(p):
    for ts, te, s, fam in d["tes"]:
        if ts <= p < te: return s
    return None
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
axC.set_xlabel(f"{TOP.replace('_','/')} chr{CH} position (bp) — every base at its true genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
pbadge(axC, "C", f"Base resolution, {TOP.replace('_','/')} — pooled {ztpL}   ·   5′-U = 1U   ·   5′ arrow RED = antisense-to-TE (silencing), grey = sense-to-TE", fs=7.6)
for (xbr, xcc) in [((z0 - psT) / NT_, z0), ((z1 - psT) / NT_, z1)]:
    fig.add_artist(ConnectionPatch(xyA=(xbr, off_top - 1.30), coordsA=axB.transData, xyB=(xcc, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.9, ls=(0, (3, 2))))
fig.suptitle(f"{TELAB} → strain-private piRNA SOURCE LOCUS (individual piRNAs, NOT a PICB cluster) across the 16-strain pangenome   ·   {CARRIER.replace('_','/')} chr{CHc}:{S:,}-{E:,}", fontsize=11.0, fontweight="bold", y=0.986, color="#1a1a1a")
pc.rtext(axA, 0.5, 0.953, [(f"expressed in {len(present)}/16 strains", pc.C_DUAL, True), ("·", "#bbb", False), (f"genetically absent in {nabs}", (pc.C_SILENCE if nabs else pc.C_META), bool(nabs)), ("·", "#bbb", False), (f"present-but-silent (regulatory) in {nsil}", "#5a6b7a", False)], fs=7.6, transform=fig.transFigure, center=True)
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/12_creation_source_loci/data/source_data"; _os.makedirs(_SD,exist_ok=True); FPM.to_csv(f"{_SD}/SourceData_{OUT}.csv")
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png ({len(present)} expressing strains)")
