#!/usr/bin/env python3
"""Shared helpers for the PICB-cluster locus figures. DATA SOURCES (only these, NO per-figure liftover):
  - pangenome cluster projection  : cluster_pav/picb_pangenome_clusters.tsv  (PICB-COMBINED clusters carrying
                                     BOTH own-genome and GRCm39 coords; built once by build_picb_pangenome_fpm.py)
  - sRNA BAMs                      : results/STAR_srna_strain_wise/{X}/{X}-{tp}.{rep}  (PRIMARY alignments only)
  - TE annotation                  : resources/repeatMasker/{X}_*.out  (own genome)
  - genome presence/absence (PAV)  : cluster_pav/locus_genome_pav.tsv  (pangenome; built once)
Reads are fetched ONLY inside a strain's PICB cluster coordinates, per timepoint; a strain/timepoint with no
PICB cluster is reported as OFF (the developmental/strain 'why' signal)."""
import os, glob, subprocess, pysam, numpy as np, pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; CP = f"{U}/cluster_pav"
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]; TPLAB = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}
TPCOL = {"16.5dpc": "#4393C3", "12.5dpp": "#E8852B", "20.5dpp": "#B2182B"}   # timepoint hue (E16.5 blue / P12.5 amber / P20.5 red) — Nature-style sequential
import matplotlib.colors as _mc
def _adj(c, f):   # f>0 lighten toward white, f<0 deepen toward black — keep the SAME hue
    r, g, b = _mc.to_rgb(c)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f) if f >= 0 else (r * (1 + f), g * (1 + f), b * (1 + f))
PLUS_COL = {t: _adj(TPCOL[t], -0.22) for t in TPS}    # + strand: DEEP shade of the timepoint hue (↑)
MINUS_COL = {t: _adj(TPCOL[t], 0.50) for t in TPS}    # − strand: PALE shade of the same hue (↓)
# semantic accent palette (used consistently across all locus figures)
C_SILENCE = "#C0392B"   # antisense-to-TE = silencing-competent (the key biology)
C_1U = "#B8860B"        # 1U signature = piRNA hallmark (gold)
C_DUAL = "#1B7837"      # dual-strand / ping-pong-competent architecture
C_ZOOM = "#E8A33D"      # zoom-callout accent (ties text to the highlighted bar + connector)
C_META = "#8a8a8a"      # secondary metadata (counts, coords)
C_WILD = "#C0392B"      # wild-derived strain emphasis
C_GENE = "#34495E"      # gene-model track (slate; distinct from TE family fills)
# ---- shared TE-family colours + the per-figure STRAND / TE-composition bar helpers (used by ALL locus figures) ----
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb"}
TEPAL = ["#8dd3c7", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f", "#d9d9d9"]
def famcol(f): return TECOL.get(f) or TEPAL[sum(map(ord, str(f))) % len(TEPAL)]   # stable TE-family colour (fallback palette for families not in TECOL)
def pale(c, f=0.62): return _adj(c, f)                                            # lighten toward white (− strand / pale shade)
def famcol_text(f):                                                               # VIBRANT, readable text colour for a TE family — deepen pale family fills so the label pops on white
    r, g, b = _mc.to_rgb(famcol(f)); lum = 0.299 * r + 0.587 * g + 0.114 * b
    return (r * 0.5 / max(lum, 0.05), g * 0.5 / max(lum, 0.05), b * 0.5 / max(lum, 0.05)) if lum > 0.5 else (r, g, b)
def top_pirna_fam(fam):                                                           # top TE family by piRNA reads (header 'TE-piRNA'); accepts {fam:count} or {fam:[count,...]}
    if not fam: return "—"
    g = lambda v: v[0] if isinstance(v, (list, tuple)) else v
    return max(fam, key=lambda k: g(fam[k]))
def color_legend(leg, handles):                                                  # colour each legend LABEL in its swatch colour; very-light swatches are darkened (hue kept) so text stays readable on white
    for t, h in zip(leg.get_texts(), handles):
        try: r, g, b = _mc.to_rgb(h.get_facecolor())
        except Exception: continue
        lum = 0.299 * r + 0.587 * g + 0.114 * b
        if lum > 0.6: s = 0.42 / max(lum, 0.05); r, g, b = r * s, g * s, b * s    # too pale on white → scale down to a readable luminance, keep the hue
        t.set_color((r, g, b)); t.set_fontweight("bold")
def panelA_strand(ax, xc, h, bw, pf, col, ylo=0.1):                              # Panel A: split a (log) FPM bar by HEIGHT → solid bottom (+ strand, frac pf) keeps `col`, PALE top (− strand)
    if h <= 0: return
    sp = 10 ** (np.log10(ylo) + pf * (np.log10(max(h, ylo * 1.001)) - np.log10(ylo)))
    ax.bar(xc, h - sp, bottom=sp, width=bw, color=_adj(col, 0.62), edgecolor="none", zorder=3)
def te_strand_bar(ax, tx0, TW, y, h, fam_bd, ntot, nte):                          # Panel B: stacked ALL-locus-piRNA element — TOP sub-bar (full width) = TE families RANKED by coverage (each split by STRAND solid +/pale −) + non-TE (grey); BOTTOM sub-bar = TE-piRNA SILENCING composition, FULL width (same as TOP, colours always readable): ANTISENSE-to-TE by TE FAMILY (ranked, family colours = silenced TEs) + sense-to-TE (grey). WIDTH ∝ % on TE (= teW, 0.10·TW floor); HEIGHT = full FIXED (same for every second bar) so colours stay readable. tick = TE|non-TE
    from matplotlib.patches import Rectangle
    nat = sum(v[1] for v in fam_bd.values())                                     # antisense-to-TE reads (silencing) among the TE-mapped piRNA
    h1 = 0.085; y1 = y + 0.105; h2 = 0.075; y2 = y - 0.075; teW = TW * nte / max(1, ntot)   # TOP (families×strand) raised, BOTTOM (antisense-by-family) lowered → clear vertical GAP; h2 tall enough that even a floored low-on-TE band shows its colours
    xx = tx0                                                                      # TOP sub-bar: TE-family composition × strand, then non-TE
    for f2, v in sorted(fam_bd.items(), key=lambda kv: -kv[1][0]):
        c = v[0]; fp = v[2]; fm = c - fp; fc = famcol(f2)
        if fp > 0: ax.add_patch(Rectangle((xx, y1 - h1 / 2), TW * fp / ntot, h1, facecolor=fc, ec="white", lw=0.2, clip_on=False, zorder=5)); xx += TW * fp / ntot
        if fm > 0: ax.add_patch(Rectangle((xx, y1 - h1 / 2), TW * fm / ntot, h1, facecolor=_adj(fc, 0.55), ec="white", lw=0.2, clip_on=False, zorder=5)); xx += TW * fm / ntot
    if ntot - nte > 0: ax.add_patch(Rectangle((xx, y1 - h1 / 2), TW * (ntot - nte) / ntot, h1, facecolor="#efefef", ec="white", lw=0.2, clip_on=False, zorder=5))
    ax.add_patch(Rectangle((tx0, y1 - h1 / 2), TW, h1, facecolor="none", ec="#888", lw=0.35, clip_on=False, zorder=6))
    ax.plot([tx0 + teW, tx0 + teW], [y1 - h1 / 2 - 0.015, y1 + h1 / 2 + 0.015], color="#222", lw=0.5, clip_on=False, zorder=7)   # TE | non-TE tick (top sub-bar)
    bw2 = teW; h2s = h2; _yb = y2 + h2 / 2              # BOTTOM-bar LENGTH (width) = the FIRST/top bar's TE-coverage = %on-TE × TW (= teW), left-aligned UNDER the first bar's TE part; HEIGHT = full FIXED. So second-bar length ∝ %on-TE w.r.t. the first bar; shows the antisense/sense (silencing) split of that TE-piRNA
    if nte > 0:                                                                   # BOTTOM sub-bar = TE-piRNA SILENCING composition (spans bw2): ANTISENSE-to-TE by TE FAMILY (ranked, family colours = silenced TEs) | sense-to-TE (grey)
        xx2 = tx0
        for f2, v in sorted(fam_bd.items(), key=lambda kv: -kv[1][1]):
            if v[1] > 0: ax.add_patch(Rectangle((xx2, _yb - h2s), bw2 * v[1] / nte, h2s, facecolor=famcol(f2), ec="white", lw=0.2, clip_on=False, zorder=5)); xx2 += bw2 * v[1] / nte
        if nte - nat > 0: ax.add_patch(Rectangle((xx2, _yb - h2s), bw2 * (nte - nat) / nte, h2s, facecolor="#cfcfcf", ec="white", lw=0.2, clip_on=False, zorder=5))   # sense-to-TE (grey)
        ax.add_patch(Rectangle((tx0, _yb - h2s), bw2, h2s, facecolor="none", ec="#999", lw=0.3, clip_on=False, zorder=6))   # outline = bw2 (TE-coverage width, floored)
def rtext(ax, x, y, segs, fs, va="center", sep="   ", transform=None, zorder=6, center=False):
    """Draw a logical line as left-to-right coloured/weighted SEGMENTS (rich text). segs=[(text,colour,bold),...].
    Anchored ha='left' at x (or centred on x if center=True); returns the x just past the last segment. Requires the
    axes' limits fixed first (stable transData). Each biological token gets its own colour/weight."""
    tr = transform or ax.transData; fig = ax.figure
    try: r = fig.canvas.get_renderer()
    except Exception: fig.canvas.draw(); r = fig.canvas.get_renderer()
    cx = x; arts = []
    for i, (t, c, b) in enumerate(segs):
        to = ax.text(cx, y, (sep if i else "") + t, fontsize=fs, ha="left", va=va, color=c, fontweight=("bold" if b else "normal"), transform=tr, zorder=zorder)
        arts.append(to); cx = tr.inverted().transform((to.get_window_extent(renderer=r).x1, 0))[0]
    if center:
        shift = (x - cx) / 2.0
        for to in arts: to.set_x(to.get_position()[0] + shift)
    return cx
def pbadge(ax, L, title, fs=8.2, y=1.05, tx=0.042, bx=0.012):
    """Nature-style panel header: filled circular letter chip + bold title (axes-fraction coords)."""
    ax.text(bx, y, f"{L}", transform=ax.transAxes, fontsize=9.5, fontweight="bold", color="white", ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.30", fc="#2C3E50", ec="none"), zorder=10)
    ax.text(tx, y, title, transform=ax.transAxes, fontsize=fs, fontweight="bold", color="#1a1a1a", ha="left", va="center")
def draw_strain_block(axB, per, tes, oc, ps, pe, fps, off, z0, z1, name=None, is_top=True, wild=False, TECOL=None, only_tp=None, nb=200):
    """Render ONE strain's Panel-B block in 0–1 fractional x — the CANONICAL multi-figure style, shared by the
    multi / main / single locus figures so they look identical. Draws: zoom highlight over the profiles (is_top),
    per-timepoint coverage (height = expression, colour = timepoint, deep ↑+ / pale ↓−), framed TE + gene lanes
    with chip labels, strain name + coord + summary, and the per-tp right margin (strand-split bar + DUAL/reads/
    FPM/AS→TE/1U). Pooled values derived from `per`. only_tp: if set, ▸ in C marks just that tp (single variant).
    Returns (ntot, domTE, arch, pct_AS, pct_1U)."""
    from matplotlib.patches import Rectangle
    TECOL = TECOL or {}; N = max(1, pe - ps); xr = np.linspace(0, 1, nb)
    on = [per[t] for t in TPS if per[t] is not None and per[t]["ntot"] > 0]
    npl = sum(p["nplus"] for p in on); nmi = sum(p["nminus"] for p in on); ntot = npl + nmi
    nat = sum(p["nat"] for p in on); nte = sum(p["nte"] for p in on); n1u = sum(p["n1u"] for p in on)
    domTE = dom_te_family(tes, ps, pe); arch = "dual-strand" if min(nmi, ntot - nmi) / max(1, ntot) > 0.2 else "uni-strand"
    pm = max([max(per[t]["plus"].max(), per[t]["minus"].max()) for t in TPS if per[t] is not None] + [1.0])
    BL = {"16.5dpc": off + 0.85, "12.5dpp": off + 0.0, "20.5dpp": off - 0.85}; HS = 0.38
    if is_top:
        zr0 = (z0 - ps) / N; zr1 = (z1 - ps) / N
        axB.add_patch(Rectangle((zr0, off - 1.30), zr1 - zr0, 2.60, facecolor="#FDB863", alpha=0.30, edgecolor="#E8A33D", lw=0.5, zorder=0.5))   # ▸ in C marker moved to the right-margin per-tp summary (was overlapping the coverage peak)
    for _tp in TPS:
        _b = BL[_tp]; _pt = per[_tp]
        axB.axhline(_b, color="#dddddd", lw=0.3, zorder=1)
        axB.text(-0.015, _b, TPLAB[_tp], fontsize=4.6, ha="right", va="center", color=(TPCOL[_tp] if (_pt is not None and _pt["ntot"]) else "#ccc"), fontweight="bold")
        if _pt is None or _pt["ntot"] == 0:
            axB.text(0.5, _b, "OFF — no PICB cluster", fontsize=4.2, ha="center", va="center", color="#cccccc", style="italic"); continue
        axB.fill_between(xr, _b, _b + _pt["plus"] / pm * HS, color=PLUS_COL[_tp], alpha=0.92, step="mid", lw=0)
        axB.fill_between(xr, _b, _b - _pt["minus"] / pm * HS, color=MINUS_COL[_tp], alpha=0.92, step="mid", lw=0)
    axB.add_patch(Rectangle((0, off - 1.74), 1.0, 0.42, facecolor="#f6f6f6", edgecolor="#e3e3e3", lw=0.4, zorder=0))
    for (ts, te, st, fam) in tes:
        x0 = max(0.0, (ts - ps) / N); x1 = min(1.0, (te - ps) / N)
        if x1 > x0:
            fam2 = fam.split("|")[-1] if "|" in fam else fam
            axB.add_patch(Rectangle((x0, off - 1.70), x1 - x0, 0.34, facecolor=TECOL.get(fam2, "#dddddd"), edgecolor="white", lw=0.25, zorder=2))
            if x1 - x0 > 0.014:
                af, at = (x0, x1) if st == "+" else (x1, x0)
                axB.annotate("", xy=(at, off - 1.53), xytext=(af, off - 1.53), arrowprops=dict(arrowstyle="-|>", color="#111", lw=0.45, mutation_scale=5.0), zorder=3)
    axB.text(-0.075, off - 1.53, "TEs", fontsize=6.6, ha="right", va="center", color="#333", fontweight="bold", bbox=dict(boxstyle="round,pad=0.2", fc="#ececec", ec="#c9c9c9", lw=0.5))
    axB.add_patch(Rectangle((0, off - 2.24), 1.0, 0.40, facecolor="#edf2f9", edgecolor="#d8e0ec", lw=0.4, zorder=0))
    _genes = genes_at(name, oc, ps, pe) if name else []
    _gw = [min(1.0, (ge - ps) / N) - max(0.0, (gs - ps) / N) for (gs, ge, gst, gn, gb) in _genes]
    _sym = lambda nm: not str(nm).startswith(("ENS", "MGP_"))   # real gene symbol vs bare Ensembl/MGP stable ID (covers ENSMUS/ENSMSP + strain prefixes ENSTCUG[CAST]/ENSLUMG[PWK]/...)
    _cand = [i for i in range(len(_genes)) if _sym(_genes[i][3]) or _gw[i] >= 0.04]
    _lbl = set(sorted(_cand, key=lambda i: (-int(_sym(_genes[i][3])), -_gw[i]))[:5])   # NAME real symbols always; predicted IDs only if wide; ≤5 — declutter dense loci (not every gene)
    _gi = 0
    for i, (gs, ge, gst, gname, gbt) in enumerate(_genes):
        gx0 = max(0.0, (gs - ps) / N); gx1 = min(1.0, (ge - ps) / N)
        if gx1 <= gx0: continue
        axB.add_patch(Rectangle((gx0, off - 2.20), gx1 - gx0, 0.32, facecolor="#c9d6ea", edgecolor=C_GENE, lw=0.8, zorder=2))
        if gx1 - gx0 > 0.014:
            gaf, gat = (gx0, gx1) if gst == "+" else (gx1, gx0)
            axB.annotate("", xy=(gat, off - 2.04), xytext=(gaf, off - 2.04), arrowprops=dict(arrowstyle="-|>", color=C_GENE, lw=0.6, mutation_scale=5.0), zorder=3)
        if i in _lbl:                                           # gene NAME centred under the gene, staggered over 3 rows
            axB.text((gx0 + gx1) / 2, off - 2.40 - 0.17 * (_gi % 3), gname, fontsize=4.3, ha="center", va="top", color=C_GENE, style="italic", fontweight="bold"); _gi += 1
    axB.text(-0.075, off - 2.04, "genes", fontsize=6.6, ha="right", va="center", color="white", fontweight="bold", bbox=dict(boxstyle="round,pad=0.2", fc=C_GENE, ec=C_GENE, lw=0.5))
    if name is not None:
        axB.text(-0.015, off + 1.80, name.replace("_", "/") + ("  ▼ zoom" if is_top else ""), fontsize=7.0, ha="right", va="center", fontweight="bold", color=(C_WILD if wild else "#1a1a1a"),
                 bbox=(dict(boxstyle="round,pad=0.16", fc="#FDF3E3", ec=C_ZOOM, lw=0.7) if is_top else None))
        axB.text(-0.015, off + 1.40, f"chr{oc}:{ps:,}", fontsize=4.6, ha="right", va="center", color=C_META, style="italic")
    rtext(axB, 1.015, off + 1.80, [(f"{ntot:,} primary reads", "#444", True), ("·  TE-DNA (genome)", "#999", False), (f"{domTE or '—'}", famcol_text(domTE or '—'), True)], fs=6.0)   # strain-level genomic only; TE-piRNA source + AS→TE silenced family are PER-TIMEPOINT (in each row)
    tx0 = 1.035; TW = 0.26                                      # right margin (snug after coverage): ALL-locus-piRNA bar = TE families ranked by coverage + non-TE, each split by STRAND (solid + / pale −)
    for _tp in TPS:
        _y = BL[_tp]; pt = per[_tp]
        if pt is None or pt["ntot"] == 0: continue
        _p = pt["nplus"]; _m = pt["nminus"]; _t = _p + _m
        te_strand_bar(axB, tx0, TW, _y, 0.16, pt["fam_bd"], pt["ntot"], pt["nte"])
        _pct_te = 100 * pt["nte"] / pt["ntot"]; _as_te = 100 * sum(v[1] for v in pt["fam_bd"].values()) / max(1, pt["nte"]); _dual = min(_p, _m) / max(1, _t) > 0.2; _1u = 100 * pt["n1u"] / max(1, _t); _fpm = None if fps is None else fps.get(_tp, 0); _u1 = pt["u1"]; _u1pk = ((np.argmax(_u1) + 0.5) / len(_u1) * (pe - ps) / 1000) if _u1.sum() > 0 else None   # 1U-peak offset (kb) from cluster start = dominant primary-processing site
        _topf = max(pt["fam_bd"], key=lambda f: pt["fam_bd"][f][0]) if pt["fam_bd"] else "—"                                                  # per-tp top TE-piRNA family
        _asf_tp = max(pt["fam_bd"], key=lambda f: pt["fam_bd"][f][1]) if (pt["fam_bd"] and any(v[1] for v in pt["fam_bd"].values())) else "—"   # per-tp top antisense-to-TE (silenced) family
        _inc = is_top and ((only_tp == _tp) if only_tp else any(r[0] < z1 and r[1] > z0 for r in pt["reads"]))
        _tc = TPCOL[_tp]   # colour-code this timepoint's stats row in the timepoint colour (family names + AS→TE keep their semantic colours)
        _segs = [(f"{_pct_te:.0f}% on TE", _tc, _pct_te >= 40), ("piRNA", _tc, False), (f"{_topf}", famcol_text(_topf), True), (f"{_as_te:.0f}% AS→", "#C0392B" if _as_te >= 40 else pale("#C0392B", 0.35), True), (f"{_asf_tp}", famcol_text(_asf_tp), True), ("↑↓DUAL" if _dual else "→UNI", _tc, _dual)] + ([(f"FPM {_fpm:.0f}" if _fpm >= 1 else f"FPM {_fpm:.1f}", _tc, True)] if _fpm is not None else []) + [(f"{_1u:.0f}%1U", _tc, False)] + ([(f"1U▲+{_u1pk:.1f}kb", _tc, True)] if _u1pk is not None else []) + ([("→C", C_ZOOM, True)] if _inc else [])
        rtext(axB, tx0, _y - 0.17, _segs, fs=4.1, sep=" · ", va="top")   # per-tp stats UNDER the bar (left-aligned at bar start), with a clear gap below the bar → no overlap, no right-edge cut-off
    return ntot, domTE, arch, 100 * nat / max(1, nte), 100 * n1u / max(1, ntot)
COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
_PANG = None; _PAV = None

def _pang():
    global _PANG
    if _PANG is None:
        tsv = f"{CP}/picb_pangenome_clusters.tsv"; pkl = tsv + ".pkl"   # pickle cache: ~2 s vs ~15 s CSV load
        if os.path.exists(pkl) and os.path.getmtime(pkl) >= os.path.getmtime(tsv):
            _PANG = pd.read_pickle(pkl)
        else:
            _PANG = pd.read_csv(tsv, sep="\t", dtype={"g39_chrom": str, "own_chrom": str}); _PANG.to_pickle(pkl)
    return _PANG

def clusters_at(g39c, g39s, g39e):
    """All PICB-combined clusters projecting into the GRCm39 window [g39s,g39e), with own-genome coords."""
    P = _pang()
    return P[(P.g39_chrom == str(g39c)) & (P.end > g39s) & (P.start < g39e)].copy()

def present_strains(sub, order):
    """Strains (in canonical `order`) that have a PICB cluster anywhere in the window."""
    have = set(sub.strain.unique())
    return [X for X in order if X in have]

def fpm_by_tp(sub, strain):
    """{tp: max all_primary_FPM} for a strain at this locus (Panel-A bars / per-tp on/off)."""
    d = sub[sub.strain == strain]
    return {tp: float(d[d.tp == tp].all_primary_FPM.max()) for tp in TPS if (d.tp == tp).any()}

def cluster_extent(sub, strain, tp=None, gap=30000):
    """Own-genome extent of a strain's cluster at this locus (optionally one tp). The pangenome projection is
    MANY-to-ONE for repetitive sequence (several of the strain's clusters / paralogous repeats project into one
    g39 window), so a raw union(min,max) over-broadens (gives 0-read regions or multi-Mb spans). Instead: merge
    projected fragments into contiguous groups (gap-capped, per chrom) and return the group with the MOST projected
    coverage (the true colinear ortholog, not a scattered paralog). Returns (own_chrom,start,end,allFPM,uniqFPM,strand) or None."""
    d = sub[sub.strain == strain]
    if tp is not None: d = d[d.tp == tp]
    if len(d) == 0: return None
    d = d.sort_values(["own_chrom", "own_start"]); groups = []; last = None
    for _, r in d.iterrows():
        oc, s, e, f, u, st = str(r.own_chrom), int(r.own_start), int(r.own_end), float(r.all_primary_FPM), float(r.uniq_FPM), r.strand
        if last is not None and last[0] == oc and s <= last[2] + gap:
            last[2] = max(last[2], e); last[3] = max(last[3], f); last[4] = max(last[4], u); last[5].append(st); last[6] += (e - s)
        else:
            last = [oc, s, e, f, u, [st], e - s]; groups.append(last)
    if len(groups) == 1:
        g = groups[0]
    else:                                          # disambiguate paralogs by ACTUAL READS (spurious projections have ~0)
        g = max(groups, key=lambda x: (_quick_reads(strain, x[0], x[1], x[2]), x[6], x[3]))
    return (g[0], g[1], g[2], g[3], g[4], max(set(g[5]), key=g[5].count))

def _quick_reads(strain, oc, s, e):
    """Primary 25-32 nt sRNA read count in a candidate cluster region (pooled over tp/rep) — to pick the real
    cluster over a spurious repeat-paralog projection that carries no reads."""
    n = 0; BAMC = f"{strain}#1#chr{oc}"
    for tp in TPS:
        for r in (1, 2, 3):
            b = f"{ROOT}/results/STAR_srna_strain_wise/{strain}/{strain}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
            if not os.path.exists(b): continue
            try:
                bam = pysam.AlignmentFile(b, "rb")
                for a in bam.fetch(BAMC, s, e):
                    if a.is_unmapped or a.is_secondary: continue
                    if 25 <= a.reference_end - a.reference_start <= 32: n += 1
                bam.close()
            except (OSError, ValueError):
                try: bam.close()
                except Exception: pass
    return n

def te_at(strain, own_chrom, ps, pe):
    """RepeatMasker TEs overlapping the own-genome cluster (own coords). [(start,end,strand,'class/family')]
    (the 4th field is RM .out col 11 = class/family, e.g. 'LTR/ERVK'; the legacy split('|')[-1] downstream is a no-op here)."""
    outf = glob.glob(f"{ROOT}/resources/repeatMasker/{strain}_*.out"); tes = []
    if outf:
        BAMC = f"{strain}#1#chr{own_chrom}"
        aw = subprocess.run("awk -v c=%s -v s=%d -v e=%d 'NR>3 && $5==c && $7>s && $6<e{st=($9==\"C\")?\"-\":\"+\"; print $6\"\\t\"$7\"\\t\"st\"\\t\"$11}' %s" % (BAMC, ps, pe, outf[0]), shell=True, capture_output=True, text=True).stdout
        for ln in aw.splitlines():
            f = ln.split("\t"); tes.append((int(f[0]), int(f[1]), f[2], f[3]))
    return sorted(tes)

def dom_te_family(tes, ps, pe):
    """Dominant TE family (by bp covered) across the cluster — the 'TE annotation' label."""
    from collections import defaultdict
    cov = defaultdict(int)
    for s, e, st, fam in tes:
        fam2 = fam.split("|")[-1] if "|" in fam else fam
        cov[fam2] += max(0, min(e, pe) - max(s, ps))
    return max(cov, key=cov.get) if cov else None

def genes_at(strain, own_chrom, ps, pe):
    """GENE-model features from the strain's gene GFF (resources/annotation/{strain}_v3.3.gff3) overlapping [ps,pe].
    SAME chr-level assembly as the BAM/PICB coords — VERIFIED: the GFF 'region' lengths for {strain}#1#chr{N} match
    the BAM @SQ byte-for-byte, so chrom = f'{strain}#1#chr{own_chrom}' and coordinates are identical (NO liftover).
    Returns [(start, end, strand, name, biotype), ...] for gene / pseudogene / ncRNA_gene features."""
    import subprocess
    gff = next((f"{ROOT}/resources/annotation/{strain}_v{v}.gff3" for v in ("3.3", "3.5", "3.2")
                if os.path.exists(f"{ROOT}/resources/annotation/{strain}_v{v}.gff3")), None)
    if gff is None: return []
    chrom = f"{strain}#1#chr{own_chrom}"
    try:
        out = subprocess.run(["grep", "-P", rf"^{chrom}\t[^\t]*\t(gene|pseudogene|ncRNA_gene)\t", gff], capture_output=True, text=True, timeout=180).stdout
    except Exception:
        return []
    genes = []
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) < 9: continue
        gs, ge = int(p[3]), int(p[4])
        if ge < ps or gs > pe: continue
        attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        genes.append((gs, ge, p[6], attrs.get("Name") or attrs.get("gene_id") or attrs.get("ID", "?"), attrs.get("biotype", p[2])))
    return genes

def fetch_primary(strain, own_chrom, ps, pe, tp, nb=200):
    """PRIMARY sRNA reads (25-32 nt, exclude secondary/multimap spillover) inside the cluster, one timepoint.
    Strand split = GENOMIC strand (architecture). AS->TE = antisense relative to the overlapping TE (sense/antisense)."""
    BAMC = f"{strain}#1#chr{own_chrom}"; N = max(1, pe - ps); tes = te_at(strain, own_chrom, ps, pe)
    def tst(p):
        for ts, te, s, fam in tes:
            if ts <= p < te: return s, (fam.split("|")[-1] if "|" in fam else fam)   # (TE strand, TE class/family) of the TE under p
        return None, None
    plus = np.zeros(nb); minus = np.zeros(nb); u1 = np.zeros(nb); reads = []; nat = nte = n1u = npl = nmi = 0; fam_bd = {}   # u1 = per-position density of 1U (5′-U) read 5′-ENDS → its argmax = dominant primary-piRNA processing site within the cluster
    for r in (1, 2, 3):
        b = f"{ROOT}/results/STAR_srna_strain_wise/{strain}/{strain}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
        if not os.path.exists(b): continue
        try:
            bam = pysam.AlignmentFile(b, "rb")
            for a in bam.fetch(BAMC, ps, pe):
                if a.is_unmapped or a.is_secondary or not a.query_sequence: continue   # PRIMARY only
                L = a.reference_end - a.reference_start
                if not (25 <= L <= 32): continue
                b0 = max(0, min(nb - 1, int((a.reference_start - ps) / N * nb))); b1 = max(0, min(nb, int((a.reference_end - ps) / N * nb)))
                seq = a.query_sequence.upper(); rev = a.is_reverse
                (minus if rev else plus)[b0:b1] += 1
                reads.append((a.reference_start, a.reference_end, rev, seq))
                is1u = ((COMP.get(seq[-1], "N") if rev else seq[0]) == "T"); n1u += is1u
                if is1u:                                                          # bin the 1U read by its 5′-END (the piRNA processing cut site), not its whole span
                    p5 = (a.reference_end - 1) if rev else a.reference_start; u1[max(0, min(nb - 1, int((p5 - ps) / N * nb)))] += 1
                if rev: nmi += 1
                else: npl += 1
                s, fam2 = tst((a.reference_start + a.reference_end) // 2)
                if s is not None:
                    is_as = (s == "-") != rev; nte += 1; nat += is_as
                    e = fam_bd.setdefault(fam2, [0, 0, 0]); e[0] += 1; e[1] += is_as; e[2] += (0 if rev else 1)   # per-TE-family: [reads over it, antisense-to-TE subset, +strand subset] — for the ranked TE-composition bar (split by strand)
            bam.close()
        except (OSError, ValueError):
            try: bam.close()
            except Exception: pass
    return dict(plus=plus, minus=minus, u1=u1, reads=reads, nplus=npl, nminus=nmi, ntot=npl + nmi,
               n1u=n1u, nat=nat, nte=nte, fam_bd=fam_bd, tes=tes, ps=ps, pe=pe, N=N, BAMC=BAMC, own_chrom=own_chrom)

def genome_pav(g39c, g39s, g39e, strain):
    """Is the locus SEQUENCE present in the strain's genome (pangenome PAV)? 'present'/'absent'/'?' .
    Distinguishes genetic LOSS (absent) from regulatory SILENCING (present but no cluster). Reads the
    precomputed locus_genome_pav.tsv (built once via pangenome). Returns '?' if not yet computed."""
    global _PAV
    if _PAV is None:
        _parts = [pd.read_csv(p, sep="\t", dtype={"g39_chrom": str}) for p in (f"{CP}/locus_genome_pav.tsv", f"{CP}/locus_genome_pav_divergence.tsv") if os.path.exists(p)]
        if not _parts: return "?"
        _PAV = pd.concat(_parts, ignore_index=True)
    r = _PAV[(_PAV.g39_chrom == str(g39c)) & (_PAV.g39_start <= g39s) & (_PAV.g39_end >= g39e) & (_PAV.strain == strain)]
    if len(r) == 0: return "?"
    return "present" if bool(r.present.iloc[0]) else "absent"
