#!/usr/bin/env python3
"""piRNA PANGENOME ATLAS of 16 inbred mouse strains — a Figure-1-style overview INSPIRED BY Helmy et al.
(Cell Genomics 2026; the 17-genome mouse reference pangenome). It is the piRNA counterpart of that
genomic pangenome: where the paper maps NON-REFERENCE SEQUENCE across the genome, we map STRAIN-PRIVATE
piRNA loci — the accessory piRNA repertoire.
(A) Genome-wide strain-private piRNA landscape: per-2Mb-bin density along all chromosomes for the four
    wild-derived strains (CAST/PWK/SPRET/WSB), the piRNA analog of the paper's per-chromosome
    non-reference tracks (Fig 1A/C).
(B) Genuinely-unique piRNA yield per strain (conserved-but-silent + strain-private), wild >> classical.
(C) piRNA locus FREQUENCY SPECTRUM (how many strains carry the homologous locus) — a population-genetics
    site-frequency-spectrum analog: a conserved CORE + a large strain-PRIVATE tail (cf. Fig 1D size dist).
(D) TE-family drivers of the strain-private loci — the transposon substrate seeding new piRNA source loci
    (the piRNA analog of the paper's defense/immunity functional enrichment, Fig 1E).
(F) TOTAL genuinely-unique piRNA expression across the frequency spectrum, by category — antisense-to-TE (silencing) / sense-to-TE / non-TE (other loci).
(G) genomic regions covered by TE — the TE CLASS (LTR/LINE/SINE) covering the locus, across the spectrum.
(H) per-TE-family piRNA output split by strand — antisense (silencing) vs sense.
Panels F-H use the per-candidate sense/antisense-to-TE calls (sense_antisense/…_percand) for genuinely-unique
(conserved-but-silent + strain-private) TE-annotated candidates, weighted by per-SEQUENCE expression
(Σ RPM = read count / 24-32 nt library size × 1e6, per strain × timepoint; compute_percand_rpm.py) — the
biological standard for piRNA abundance, with antisense-to-TE = the TE-silencing species."""
import sys, os, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
TH = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/25_pirna_pangenome_atlas"
SD = f"{TH}/data/source_data"; os.makedirs(SD, exist_ok=True)
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
WILD_ORD = [s for s in CANON if s in WILD]          # canonical order among wild
BIN = 2_000_000
CHROMS = [str(i) for i in range(1, 20)] + ["X"]
CLEN = {"1":195154279,"2":181755017,"3":159745316,"4":156860686,"5":151758149,"6":149588044,"7":144995196,
        "8":130127694,"9":124359700,"10":130530862,"11":121973369,"12":120092757,"13":120883175,"14":125139656,
        "15":104073951,"16":98008968,"17":95294699,"18":90720763,"19":61420004,"X":169476592}
WCOL = {"CAST_EiJ":"#009E73","PWK_PhJ":"#D55E00","SPRET_EiJ":"#7a3b9a","WSB_EiJ":"#0072B2"}

# ---- klass5: genuinely-unique sets + per-strain yield + locus frequency spectrum ----
k = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",
                usecols=["sequence","strain","klass5","homolog_strains"])
PRIV = "unique: strain-private locus"; CBS = "unique: conserved-but-silent"
priv_by_strain = {X: set(k.loc[(k.strain==X)&(k.klass5.isin([PRIV,CBS])),"sequence"]) for X in CANON}   # genuinely-unique loci that project to GRCm39 (strain-private new loci are mostly OFF-reference, like the paper's non-reference sequence)
uniq = k[k.klass5.isin([PRIV,CBS])]
yield_tab = (uniq.groupby(["strain","klass5"]).size().unstack(fill_value=0).reindex(index=CANON, columns=[CBS,PRIV]).fillna(0))
# locus frequency = # strains carrying the homologous locus (from homolog_strains), de-duplicated per (klass, sequence)
_fs = uniq.drop_duplicates(["klass5","sequence"]).copy()
_fs["nstrains"] = _fs.homolog_strains.fillna("").apply(lambda s: len([x for x in str(s).split(",") if x]) if s else 1).clip(1,16)
spec_priv = np.bincount(_fs.loc[_fs.klass5==PRIV,"nstrains"], minlength=17)[1:17]
spec_cbs  = np.bincount(_fs.loc[_fs.klass5==CBS ,"nstrains"], minlength=17)[1:17]

# ---- per-2Mb-bin strain-private density from GRCm39-projected loci ----
def strain_density(X):
    bed = f"{U}/unique16/loci/{X}.cand_GRCm39.bed"; dens = {c: np.zeros(int(np.ceil(CLEN[c]/BIN))) for c in CHROMS}
    seen = set()
    if os.path.exists(bed):
        want = priv_by_strain[X]
        for ln in open(bed):
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4: continue
            c = f[0]
            if c not in CHROMS: continue
            seq = f[3].split("|")[-1]
            if seq not in want: continue
            key = (c, int(f[1]))
            if key in seen: continue                 # count each locus once (a seq maps once)
            seen.add(key)
            dens[c][int(f[1])//BIN] += 1
    return dens
DENS = {X: strain_density(X) for X in CANON}
perchrom = pd.DataFrame({X: {c: DENS[X][c].sum() for c in CHROMS} for X in CANON}).T.reindex(CANON)[CHROMS]

# ---- TE-family drivers of strain-private loci ----
te = pd.read_csv(f"{U}/pangenome_te/SourceData_TE_private_families16_byclass.csv")
tep = te[te.klass=="strain-private"].copy()
tep["grp"] = tep.strain.map(lambda s: "wild" if s in WILD else "classical")
topfam = (tep[~tep.family.astype(str).str.startswith("__")]           # drop __nseen__ (unannotated) / __nte__ (no-TE overlap) placeholders
          .groupby("family")["count"].sum().sort_values(ascending=False).head(8).index.tolist())
famtab = tep[tep.family.isin(topfam)].groupby(["family","grp"])["count"].sum().unstack(fill_value=0).reindex(topfam)

# ---- Panel E: a strain-private TE-driven piRNA locus (SPRET/EiJ, antisense to a young ERVK — the TE-insertion-gain origin) ----
import pysam as _ps
_ELOC = ("SPRET_EiJ", "13", 46306947, 46306978)                # strain-private klass5 piRNA, pachytene 20.5dpp, antisense to the TE (silencing)
_ETE  = (46306428, 46307217, "RLTR22_Mur", "LTR/ERVK")         # young ERVK it sits in (SPRET RepeatMasker); the piRNA locus is absent in all classical strains (locus gain)
_ea, _eb = _ETE[0]-260, _ETE[1]+260
_ebam = _ps.AlignmentFile(f"{U.split('/analysis/')[0]}/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam", "rb")
_ecov = np.zeros(_eb-_ea)
for _r in _ebam.fetch("SPRET_EiJ#1#chr13", _ea, _eb):
    if _r.is_unmapped or not 25 <= _r.reference_end-_r.reference_start <= 32: continue
    for _p in range(max(_ea, _r.reference_start), min(_eb, _r.reference_end)): _ecov[_p-_ea] += 1
_ebam.close(); _emax = max(_ecov.max(), 1.0)

# ---- per-piRNA-SEQUENCE expression (RPM) x strand (sense/antisense to TE) x TE, genuinely-unique repertoire ----
# RPM = read count / library size (24-32 nt window) x 1e6, per (sequence, strain, timepoint), mean over reps —
# the biological standard for piRNA abundance; antisense-to-TE piRNAs are the silencing species. Panels weight
# by this per-sequence expression (Sum RPM), NOT by locus count. Precomputed by compute_percand_rpm.py.
# _rpm = ALL genuinely-unique candidates (per-sequence RPM); _sa = the TE-overlapping subset with sense/antisense-to-TE.
_rpm = pd.read_csv(f"{U}/unique16/percand_rpm_expr.csv.gz").rename(columns={"timepoint": "tp"})
_sa = pd.read_csv(f"{U}/sense_antisense/SourceData_sense_antisense16_percand.csv.gz")
_sa["tp"] = _sa.id.str.split("|").str[1]; _sa["sequence"] = _sa.id.str.split("|").str[-1]
_sa = _sa.drop_duplicates(["sequence","strain","tp"])                                          # one TE-orientation per candidate
_all = _rpm.merge(_sa[["sequence","strain","tp","orientation","family"]], on=["sequence","strain","tp"], how="left")   # non-TE candidates -> NaN orientation
_all = _all.merge(k.drop_duplicates("sequence")[["sequence","homolog_strains","klass5"]], on="sequence", how="left")
_all["nstr"] = _all.homolog_strains.fillna("").apply(lambda s: len([x for x in str(s).split(",") if x]) if s else 1).clip(1,16)
_all["cat"] = _all.orientation.fillna("non-TE")                                                # antisense (to TE, silencing) / sense (to TE) / non-TE (other loci)
# ---- Panel F: TOTAL genuinely-unique piRNA expression by category across the frequency spectrum ----
strand_spec = _all.groupby(["nstr","cat"])["rpm"].sum().unstack(fill_value=0).reindex(index=range(1,17), fill_value=0)
for _c in ("antisense","sense","non-TE"):
    if _c not in strand_spec.columns: strand_spec[_c] = 0
anti_pct = (100*strand_spec["antisense"]/(strand_spec["antisense"]+strand_spec["sense"]).replace(0, np.nan)).reindex(range(1,17))   # silencing share of the TE-overlapping part — Panel F line
tp_cat = _all.groupby(["tp","cat"])["rpm"].sum().unstack(fill_value=0)                          # Panel I: TOTAL by timepoint x category (antisense-to-TE / sense-to-TE / non-TE)
# ---- TE-overlapping genuinely-unique subset for Panels G / H / I ----
_gu = _all[_all.orientation.notna()].copy(); _gu["tesuper"] = _gu.family.astype(str).str.split("/").str[0]              # LTR / LINE / SINE / ...
SUPERS = ["LTR","LINE","SINE"]
te_spec = _gu.groupby(["nstr","tesuper"])["rpm"].sum().unstack(fill_value=0).reindex(index=range(1,17), fill_value=0)
te_spec["other"] = te_spec.drop(columns=[c for c in SUPERS if c in te_spec], errors="ignore").sum(axis=1)
famH = _gu.groupby("family")["rpm"].sum().sort_values(ascending=False).head(7).index.tolist()
fam_strand = _gu[_gu.family.isin(famH)].groupby(["family","orientation"])["rpm"].sum().unstack(fill_value=0).reindex(famH)
for _o in ("antisense","sense"):
    if _o not in fam_strand.columns: fam_strand[_o] = 0
TP_ORD = ["16.5dpc","12.5dpp","20.5dpp"]; TP_LAB = {"16.5dpc":"E16.5\n(fetal)","12.5dpp":"P12.5\n(early postnatal)","20.5dpp":"P20.5\n(pachytene)"}
# (Panel I uses tp_cat — the total-by-timepoint x category aggregation above)

print(f"genuinely-unique loci (any strain, deduped seq): {sum(len(v) for v in priv_by_strain.values()):,}")
print(f"wild per-chrom private density total: {perchrom.loc[WILD_ORD].values.sum():.0f}")

# =====================  FIGURE  =====================
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})
fig = plt.figure(figsize=(24, 13.8), dpi=300)
gs = fig.add_gridspec(4, 3, width_ratios=[1.35, 1.0, 1.0], height_ratios=[1.05, 1.0, 1.0, 0.82], hspace=0.5, wspace=0.24,
                      left=0.045, right=0.99, top=0.92, bottom=0.05)

# ---- Panel A: genome-wide strain-private piRNA landscape (4 wild strains) ----
axA = fig.add_subplot(gs[0:2, 0]); axA.set_xlim(0, 1); axA.set_ylim(-1.9, len(CHROMS)); axA.axis("off")
axA.set_title("A   Genome-wide strain-specific (genuinely-unique) piRNA landscape — the four wild-derived strains",
              fontsize=10.5, fontweight="bold", loc="left")
_alld = np.array([v for X in WILD_ORD for c in CHROMS for v in DENS[X][c] if v>0]); GMAX = float(np.percentile(_alld, 88)) if _alld.size else 1.0
row_h = 1.0; sub = row_h/ (len(WILD_ORD)+0.4)
L, R = 0.10, 0.985
for ci, c in enumerate(CHROMS):
    y0 = len(CHROMS) - 1 - ci
    axA.text(L-0.012, y0+row_h*0.5, c, ha="right", va="center", fontsize=7.2, fontweight="bold")
    xr = (R-L) * CLEN[c]/CLEN["1"]                                  # scale each chrom by true length
    axA.plot([L, L+xr], [y0+0.04, y0+0.04], color="#ccc", lw=0.6, zorder=0)
    nb = len(DENS[WILD_ORD[0]][c]); xs = L + (np.arange(nb)+0.5)/nb * xr
    for wi, X in enumerate(WILD_ORD):
        base = y0 + 0.06 + wi*sub; d = DENS[X][c]
        h = np.clip(d/GMAX, 0, 1) * sub*0.60
        axA.fill_between(xs, base, base+h, step="mid", color=WCOL[X], edgecolor=WCOL[X], linewidth=0.15, zorder=2)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(facecolor=WCOL[X], label=X.replace("_","/")) for X in WILD_ORD],
           fontsize=7.5, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.004),
           title="genuinely-unique piRNA loci per 2-Mb bin (filled profile = density)", title_fontsize=8)

# ---- Panel E: a strain-private TE-driven piRNA locus (SPRET/EiJ) — real pachytene coverage over a young ERVK ----
axE = fig.add_subplot(gs[2:4, 0]); _exk = np.arange(_ea, _eb)/1e3
axE.fill_between(_exk, 0, _ecov, step="mid", color=WCOL["SPRET_EiJ"], alpha=0.62, lw=0, zorder=3)   # SPRET 20.5dpp piRNA coverage (25-32 nt)
axE.set_xlim(_ea/1e3, _eb/1e3); axE.set_ylim(-_emax*0.32, _emax*1.32)
axE.add_patch(plt.Rectangle((_ETE[0]/1e3, -_emax*0.18), (_ETE[1]-_ETE[0])/1e3, _emax*0.10, facecolor="#9e9e9e", edgecolor="#555", lw=0.5, clip_on=False, zorder=4))   # the ERVK TE span
axE.text((_ETE[0]+_ETE[1])/2e3, -_emax*0.215, f"{_ETE[2]}  ({_ETE[3]})", ha="center", va="top", fontsize=6.3, color="#555", fontweight="bold")
axE.add_patch(plt.Rectangle((_ELOC[2]/1e3, _emax*1.02), max((_ELOC[3]-_ELOC[2])/1e3, 0.008), _emax*0.08, facecolor="#C0392B", edgecolor="none", zorder=5))   # the strain-private piRNA
axE.text(_ea/1e3+0.015, _emax*1.18, "strain-private piRNA — antisense (silencing)", ha="left", va="center", fontsize=5.9, color="#C0392B", fontweight="bold")
axE.set_title("E   A strain-private TE-driven piRNA locus (SPRET/EiJ)", fontsize=9.6, fontweight="bold", loc="left")
axE.set_xlabel("chr13 position (kb)", fontsize=8.5); axE.set_ylabel("piRNA coverage (25-32 nt)", fontsize=8); axE.tick_params(labelsize=6.6)
axE.ticklabel_format(axis="x", useOffset=False, style="plain"); axE.locator_params(axis="x", nbins=5)
axE.spines[["top", "right"]].set_visible(False)
axE.text(0.5, -0.33, "A young ERVK (RLTR22_Mur) insertion private to SPRET/EiJ spawns an antisense (silencing) pachytene piRNA — the TE-insertion-gain origin of the strain-private accessory repertoire (cf. panel D). The locus is absent in all 12 classical strains.",
         transform=axE.transAxes, ha="center", fontsize=6.3, color="#666", style="italic", wrap=True)

# ---- Panel B: genuinely-unique yield per strain ----
axB = fig.add_subplot(gs[0, 1]); x = np.arange(len(CANON))
axB.bar(x, yield_tab[CBS].values, color="#9ecae8", edgecolor="white", linewidth=0.3, label="conserved-but-silent")
axB.bar(x, yield_tab[PRIV].values, bottom=yield_tab[CBS].values, color="#7a3b9a", edgecolor="white", linewidth=0.3, label="strain-private locus")
tot = (yield_tab[CBS]+yield_tab[PRIV]).values
for xi, t in zip(x, tot): axB.text(xi, t+tot.max()*0.01, f"{t/1000:.1f}k", ha="center", va="bottom", fontsize=4.6, rotation=90, color="#333")
axB.set_xticks(x); axB.set_xticklabels([s.replace("_","/") for s in CANON], rotation=90, fontsize=6.2)
for t, s in zip(axB.get_xticklabels(), CANON):
    if s in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axB.set_ylim(0, tot.max()*1.16); axB.set_ylabel("genuinely-unique piRNAs (n)", fontsize=8.5)
axB.legend(fontsize=6.6, frameon=False, loc="upper left"); axB.spines[["top","right"]].set_visible(False); axB.tick_params(labelsize=7)
axB.set_title("B   Genuinely-unique piRNA yield per strain — wild-derived (red) dominate", fontsize=9.4, fontweight="bold", loc="left")

# ---- Panel C: locus frequency spectrum ----
axC = fig.add_subplot(gs[1, 1]); xs2 = np.arange(1, 17)
axC.bar(xs2, spec_priv, color="#7a3b9a", label="strain-private locus", edgecolor="white", linewidth=0.3)
axC.bar(xs2, spec_cbs, bottom=spec_priv, color="#9ecae8", label="conserved-but-silent", edgecolor="white", linewidth=0.3)
axC.set_yscale("log"); axC.set_xticks(xs2); axC.tick_params(labelsize=7)
axC.set_ylim(top=axC.get_ylim()[1]*2.4)                              # headroom for bar-top count labels
for _xi,_tv in zip(xs2, spec_priv+spec_cbs):                         # creative: total loci labelled on every bar
    if _tv<=0: continue
    _lab = f"{_tv/1000:.1f}k".replace(".0k","k") if _tv>=1000 else f"{int(_tv)}"
    axC.text(_xi, _tv*1.22, _lab, ha="center", va="bottom", fontsize=5.2, fontweight="bold", rotation=90,
             color="#7a3b9a" if _xi==1 else "#0072B2" if _xi==16 else "#555")
axC.set_xlabel("number of strains carrying the homologous locus  (1 = private … 16 = core)", fontsize=8)
axC.set_ylabel("unique piRNA loci (log)", fontsize=8.5)
axC.axvspan(0.5,1.5,color="#7a3b9a",alpha=0.06); axC.axvspan(15.5,16.5,color="#0072B2",alpha=0.06)
axC.text(1,axC.get_ylim()[1]*0.4,"PRIVATE",ha="center",fontsize=6.3,color="#7a3b9a",fontweight="bold")
axC.text(16,axC.get_ylim()[1]*0.4,"CORE",ha="center",fontsize=6.3,color="#0072B2",fontweight="bold")
axC.legend(fontsize=6.6, frameon=False, loc="upper center"); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C   piRNA locus frequency spectrum — a conserved core + a large private tail", fontsize=9.4, fontweight="bold", loc="left")

# ---- Panel D: TE-family drivers ----
axD = fig.add_subplot(gs[2, 1]); xf = np.arange(len(topfam)); w = 0.4
axD.bar(xf-w/2, famtab.get("wild",pd.Series(0,index=topfam)).values, w, color="#C0392B", label="wild-derived", edgecolor="white", linewidth=0.3)
axD.bar(xf+w/2, famtab.get("classical",pd.Series(0,index=topfam)).values, w, color="#4393C3", label="classical", edgecolor="white", linewidth=0.3)
axD.set_xticks(xf); axD.set_xticklabels(topfam, rotation=40, ha="right", fontsize=6.8); axD.tick_params(labelsize=7)
axD.set_ylabel("strain-private piRNA loci (n)", fontsize=8.5)
axD.legend(fontsize=6.8, frameon=False, loc="upper right"); axD.spines[["top","right"]].set_visible(False)
axD.set_title("D   TE-family drivers of strain-private piRNA loci — young active retrotransposons (ERVK/L1)", fontsize=9.4, fontweight="bold", loc="left")

# ---- Panel F: piRNA strand (sense/antisense to TE) across the frequency spectrum ----
axF = fig.add_subplot(gs[0, 2])
_fcum = np.zeros(16)   # TOTAL genuinely-unique piRNA, stacked: antisense-to-TE (silencing) + sense-to-TE + non-TE
for _lab, _col, _key in [("antisense-to-TE (silencing)","#C0392B","antisense"), ("sense-to-TE","#9e9e9e","sense"), ("non-TE (other loci)","#a9c7dd","non-TE")]:
    _v = strand_spec[_key].values; axF.bar(xs2, _v, bottom=_fcum, color=_col, label=_lab, edgecolor="white", linewidth=0.3); _fcum = _fcum + _v
axF2 = axF.twinx()                                                        # % antisense of the TE-overlapping part (silencing share) on a secondary axis
axF2.plot(xs2, anti_pct.values, color="#111111", lw=1.4, marker="o", ms=2.6, zorder=6)
axF2.axhline(50, color="#111111", lw=0.5, ls=(0,(3,2)), alpha=0.5)
axF2.set_ylim(0, 100); axF2.set_ylabel("% antisense of TE-overlapping (silencing)", fontsize=7.0, color="#111111"); axF2.tick_params(labelsize=6.4)
axF2.spines[["top"]].set_visible(False)
axF.set_yscale("log"); axF.set_xticks(xs2); axF.tick_params(labelsize=7)
axF.set_xlabel("number of strains carrying the locus  (1 = private … 16 = core)", fontsize=8)
axF.set_ylabel("TOTAL genuinely-unique piRNA\nexpression (Σ RPM, log)", fontsize=8.2)
axF.axvspan(0.5,1.5,color="#7a3b9a",alpha=0.06); axF.axvspan(15.5,16.5,color="#0072B2",alpha=0.06)
axF.text(1, axF.get_ylim()[1]*0.4, "PRIVATE", ha="center", fontsize=6.3, color="#7a3b9a", fontweight="bold")
axF.text(16, axF.get_ylim()[1]*0.4, "CORE", ha="center", fontsize=6.3, color="#0072B2", fontweight="bold")
axF.legend(fontsize=6.0, frameon=False, loc="upper center"); axF.spines[["top","right"]].set_visible(False)
axF.set_title("F   Total piRNA expression across the spectrum — antisense-to-TE (silencing) / sense-to-TE / non-TE", fontsize=8.8, fontweight="bold", loc="left")

# ---- Panel G: TE class covering the piRNA locus, across the frequency spectrum ----
axG = fig.add_subplot(gs[1, 2]); _gcol = {"LTR":"#6a3d9a","LINE":"#E69F00","SINE":"#1f78b4","other":"#bbbbbb"}
_bot = np.zeros(16)
for s in SUPERS + ["other"]:
    _v = te_spec[s].values if s in te_spec.columns else np.zeros(16)
    axG.bar(xs2, _v, bottom=_bot, color=_gcol[s], label=s, edgecolor="white", linewidth=0.2); _bot = _bot + _v
axG.set_yscale("log"); axG.set_xticks(xs2); axG.tick_params(labelsize=7)
axG.set_xlabel("number of strains carrying the locus  (1 = private … 16 = core)", fontsize=8)
axG.set_ylabel("TE-covered piRNA expression\n(Σ RPM, log)", fontsize=8.2)
axG.axvspan(0.5,1.5,color="#7a3b9a",alpha=0.06); axG.axvspan(15.5,16.5,color="#0072B2",alpha=0.06)
axG.text(1, axG.get_ylim()[1]*0.4, "PRIVATE", ha="center", fontsize=6.3, color="#7a3b9a", fontweight="bold")
axG.text(16, axG.get_ylim()[1]*0.4, "CORE", ha="center", fontsize=6.3, color="#0072B2", fontweight="bold")
axG.legend(fontsize=6.3, frameon=False, loc="upper center", ncol=2); axG.spines[["top","right"]].set_visible(False)
axG.set_title("G   Genomic regions covered by TE — expression by TE class across the spectrum", fontsize=9.2, fontweight="bold", loc="left")

# ---- Panel H: per-TE-family sense vs antisense (TE-derived piRNA strand) ----
axH = fig.add_subplot(gs[2, 2]); xh = np.arange(len(famH)); wh = 0.4
axH.bar(xh-wh/2, fam_strand["antisense"].values, wh, color="#C0392B", label="antisense (silencing)", edgecolor="white", linewidth=0.3)
axH.bar(xh+wh/2, fam_strand["sense"].values, wh, color="#9e9e9e", label="sense", edgecolor="white", linewidth=0.3)
axH.set_xticks(xh); axH.set_xticklabels(famH, rotation=40, ha="right", fontsize=6.6); axH.tick_params(labelsize=7)
axH.set_ylabel("TE-family piRNA expression\n(Σ RPM)", fontsize=8.2)
axH.legend(fontsize=6.5, frameon=False, loc="upper right"); axH.spines[["top","right"]].set_visible(False)
axH.set_title("H   TE-family piRNA expression — antisense (silencing) vs sense", fontsize=9.2, fontweight="bold", loc="left")

# ---- Panel I: developmental timepoint — the SAME total (3-category) across spermatogenesis (fetal -> pachytene) ----
axI = fig.add_subplot(gs[3, 1:3]); axI2 = axI.twinx()
xI = np.arange(len(TP_ORD)); wI = 0.5; _cum = np.zeros(len(TP_ORD))
for _key, _col, _lab in [("antisense","#C0392B","antisense-to-TE (silencing)"), ("sense","#9e9e9e","sense-to-TE"), ("non-TE","#a9c7dd","non-TE (other loci)")]:
    _v = np.array([float(tp_cat.loc[tp, _key]) if (tp in tp_cat.index and _key in tp_cat.columns) else 0.0 for tp in TP_ORD])
    axI.bar(xI, _v, wI, bottom=_cum, color=_col, edgecolor="white", lw=0.3, label=_lab); _cum = _cum + _v
_ant = np.array([float(tp_cat.loc[tp,"antisense"]) if (tp in tp_cat.index and "antisense" in tp_cat.columns) else 0.0 for tp in TP_ORD])
_sen = np.array([float(tp_cat.loc[tp,"sense"]) if (tp in tp_cat.index and "sense" in tp_cat.columns) else 0.0 for tp in TP_ORD])
axI2.plot(xI, np.where((_ant+_sen) > 0, 100*_ant/(_ant+_sen), np.nan), color="#111111", lw=1.6, marker="o", ms=5.5, zorder=6)
axI2.set_ylim(0, 100); axI2.set_ylabel("% antisense of TE-overlapping (silencing)", fontsize=7.3); axI2.tick_params(labelsize=6.4); axI2.spines[["top"]].set_visible(False)
axI.set_xticks(xI); axI.set_xticklabels([TP_LAB[t] for t in TP_ORD], fontsize=9)
axI.set_ylabel("TOTAL piRNA expression (Σ RPM)", fontsize=8.5); axI.tick_params(labelsize=7.5); axI.margins(x=0.22)
axI.spines[["top", "right"]].set_visible(False)
axI.legend(fontsize=7.2, frameon=False, ncol=1, loc="upper right")
axI.set_title("I   Developmental timepoint — TOTAL genuinely-unique piRNA expression (fetal → pachytene), by category; line = % antisense-to-TE (silencing)",
              fontsize=9.0, fontweight="bold", loc="left")

fig.suptitle("The piRNA PANGENOME of 16 inbred mouse strains — a conserved core piRNA-ome and a large, wild-derived-dominated, TE-driven strain-private accessory repertoire\n"
             "(the piRNA counterpart of the 17-genome mouse reference pangenome, Helmy et al., Cell Genomics 2026)",
             fontsize=11.5, fontweight="bold", y=0.975, linespacing=1.5)
for e in ("pdf","svg","png"): fig.savefig(f"{TH}/figures/Fig_pirna_pangenome_atlas.{e}", bbox_inches="tight")
# ---- source data ----
perchrom.to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_perchrom_density.csv")
yield_tab.assign(total=tot).to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_yield.csv")
pd.DataFrame({"strains_carrying":xs2,"strain_private":spec_priv,"conserved_but_silent":spec_cbs}).to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_freq_spectrum.csv",index=False)
famtab.to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_TE_drivers.csv")
pd.DataFrame({"chr13_position": np.arange(_ea,_eb), "SPRET_piRNA_coverage_25_32nt": _ecov.astype(int)}).to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_panelE_ERVK_locus.csv",index=False)   # SPRET strain-private piRNA over RLTR22_Mur (LTR/ERVK), chr13:46,306,428-46,307,217
strand_spec.rename_axis("strains_carrying").reset_index().to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_strand_spectrum.csv", index=False)     # Panel F
te_spec.rename_axis("strains_carrying").reset_index().to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_TEclass_spectrum.csv", index=False)       # Panel G
fam_strand.rename_axis("family").reset_index().to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_TEfamily_strand.csv", index=False)               # Panel H
anti_pct.rename_axis("strains_carrying").rename("pct_antisense").reset_index().to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_antisense_pct.csv", index=False)   # Panel F line
tp_cat.rename_axis("timepoint").reset_index().to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_timepoint_category.csv", index=False)               # Panel I (total by timepoint x category)
print("wrote Fig_pirna_pangenome_atlas.{png,pdf,svg} + 10 source_data files")
