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
    (the piRNA analog of the paper's defense/immunity functional enrichment, Fig 1E)."""
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
topfam = tep.groupby("family")["count"].sum().sort_values(ascending=False).head(8).index.tolist()
famtab = tep[tep.family.isin(topfam)].groupby(["family","grp"])["count"].sum().unstack(fill_value=0).reindex(topfam)

# ---- structural-diversity zoom: the chr17 ~27.5 Mbp strain-private piRNA hotspot (each wild strain a DISTINCT private locus) ----
ZC, ZS, ZE = "17", 27_500_000, 27_620_000
_pc = pd.read_csv(f"{U}/cluster_pav/picb_pangenome_clusters.tsv", sep="\t", dtype={"g39_chrom":str}, low_memory=False)
_z = _pc[(_pc.g39_chrom==ZC) & (_pc.start < ZE) & (_pc.end > ZS)].copy()
zoom = _z.sort_values("all_primary_FPM", ascending=False).groupby("strain", as_index=False).first()   # best locus per strain in the window

print(f"genuinely-unique loci (any strain, deduped seq): {sum(len(v) for v in priv_by_strain.values()):,}")
print(f"wild per-chrom private density total: {perchrom.loc[WILD_ORD].values.sum():.0f}")

# =====================  FIGURE  =====================
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})
fig = plt.figure(figsize=(17, 11), dpi=300)
gs = fig.add_gridspec(3, 2, width_ratios=[1.35, 1.0], height_ratios=[1.05, 1.0, 1.0], hspace=0.42, wspace=0.20,
                      left=0.055, right=0.985, top=0.90, bottom=0.06)

# ---- Panel A: genome-wide strain-private piRNA landscape (4 wild strains) ----
axA = fig.add_subplot(gs[0:2, 0]); axA.set_xlim(0, 1); axA.set_ylim(0, len(CHROMS)); axA.axis("off")
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
        h = np.clip(d/GMAX, 0, 1) * sub*0.92
        axA.fill_between(xs, base, base+h, step="mid", color=WCOL[X], edgecolor=WCOL[X], linewidth=0.15, zorder=2)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(facecolor=WCOL[X], label=X.replace("_","/")) for X in WILD_ORD],
           fontsize=7.5, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.02),
           title="genuinely-unique piRNA loci per 2-Mb bin (filled profile = density)", title_fontsize=8)

# ---- Panel E: structural-diversity locus zoom (chr17 strain-private piRNA hotspot) ----
axE = fig.add_subplot(gs[2, 0]); axE.set_xlim(ZS/1e6, ZE/1e6); axE.set_ylim(-0.7, len(WILD_ORD)-0.3)
axE.set_title("E   A high-diversity locus (chr17 ≈ 27.5 Mbp) — each wild strain carries a DISTINCT strain-private piRNA locus",
              fontsize=9.6, fontweight="bold", loc="left")
for yi, X in enumerate(WILD_ORD[::-1]):
    axE.text(ZS/1e6-0.003, yi, X.replace("_","/"), ha="right", va="center", fontsize=7.5, fontweight="bold", color="#C0392B")
    axE.axhline(yi, color="#eeeeee", lw=0.5, zorder=0)
    r = zoom[zoom.strain==X]
    if len(r):
        s, e, fpm = r.start.iloc[0]/1e6, r.end.iloc[0]/1e6, r.all_primary_FPM.iloc[0]
        axE.add_patch(plt.Rectangle((s, yi-0.30), max(e-s, 0.0015), 0.60, facecolor=WCOL[X], edgecolor="#222", lw=0.4, zorder=3))
        axE.annotate(f"{fpm/1000:.0f}k FPM", ((s+e)/2, yi+0.33), ha="center", va="bottom", fontsize=6, color=WCOL[X], fontweight="bold")
axE.set_yticks([]); axE.set_xlabel("chr17 position (Mbp)", fontsize=8.5); axE.tick_params(labelsize=7.5)
axE.spines[["top", "right", "left"]].set_visible(False)
axE.text(0.5, -0.34, "The piRNA analog of a high-diversity genomic locus (cf. the GBP cluster, Fig 1A): four private, non-overlapping piRNA source loci within ~80 kb — one per wild strain, silent in all 12 classical strains.",
         transform=axE.transAxes, ha="center", fontsize=6.3, color="#666", style="italic")

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
axD.set_xticks(xf); axD.set_xticklabels(topfam, rotation=35, ha="right", fontsize=7); axD.tick_params(labelsize=7)
axD.set_ylabel("strain-private piRNA loci (n)", fontsize=8.5)
axD.legend(fontsize=6.8, frameon=False, loc="upper right"); axD.spines[["top","right"]].set_visible(False)
axD.set_title("D   TE-family drivers of strain-private piRNA loci — young active retrotransposons (ERVK/L1)", fontsize=9.4, fontweight="bold", loc="left")

fig.suptitle("The piRNA PANGENOME of 16 inbred mouse strains — a conserved core piRNA-ome and a large, wild-derived-dominated, TE-driven strain-private accessory repertoire\n"
             "(the piRNA counterpart of the 17-genome mouse reference pangenome, Helmy et al., Cell Genomics 2026)",
             fontsize=11.5, fontweight="bold", y=0.975, linespacing=1.5)
for e in ("pdf","svg","png"): fig.savefig(f"{TH}/figures/Fig_pirna_pangenome_atlas.{e}", bbox_inches="tight")
# ---- source data ----
perchrom.to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_perchrom_density.csv")
yield_tab.assign(total=tot).to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_yield.csv")
pd.DataFrame({"strains_carrying":xs2,"strain_private":spec_priv,"conserved_but_silent":spec_cbs}).to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_freq_spectrum.csv",index=False)
famtab.to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_TE_drivers.csv")
zoom[["strain","g39_chrom","start","end","all_primary_FPM","strand"]].to_csv(f"{SD}/SourceData_Fig_pirna_pangenome_atlas_chr17_zoom.csv",index=False)
print("wrote Fig_pirna_pangenome_atlas.{png,pdf,svg} + 5 source_data files")
