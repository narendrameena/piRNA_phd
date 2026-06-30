#!/usr/bin/env python3
"""How much of TOTAL piRNA EXPRESSION (read-mass) is strain-private — coverage decomposition, 16 strains.
Read-mass-weighted % of each strain's total piRNA reads that fall in the strain-private set, under the
absence-rule ladder loose(>=2/3 reps) -> ADOPTED >=2 reads -> strict(0 reads). DA-only coverage panel is
added automatically once edger16/*.coverage_full.csv (edger16_coverage.R) exists.
A: per-strain strain-private coverage (3 absence rules);  B: by-group + overall summary (+DA-only when ready);
C: DA-only vs strain-private per strain (enriched-but-not-specific vs specific) — populated when coverage_full ready.
Source: edger16/*.coverage_probe.csv (+ *.coverage_full.csv)."""
import sys, glob
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
E16 = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
PG  = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
SD  = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
df = pd.concat([pd.read_csv(f) for f in glob.glob(f"{E16}/*.coverage_probe.csv")], ignore_index=True)
CANON = [s for s in STRAIN_ORDER if s in set(df.strain)]; WPOS = [i for i,s in enumerate(CANON) if s in WILD]
full = glob.glob(f"{E16}/*.coverage_full.csv"); HAVE_DA = len(full) == 3
if HAVE_DA:
    fd = pd.concat([pd.read_csv(f) for f in full], ignore_index=True)
    df = df.merge(fd[["strain","timepoint","cov_DA"]], on=["strain","timepoint"], how="left")
def w(d, col):  # read-mass-weighted mean coverage
    return (d[col]*d.total_reads).sum()/d.total_reads.sum()
def per_strain(col):
    if col not in df.columns: return pd.Series(0.0, index=CANON)
    return pd.Series({s: w(df[df.strain==s], col) for s in CANON})
loose, twor, strict = per_strain("cov_presence_loose2"), per_strain("cov_presence_2read2"), per_strain("cov_presence_strict2")
inter = per_strain("cov_intersection")                                  # DA ∩ loose (current pipeline)
inter2, interS = per_strain("cov_intersection_2read"), per_strain("cov_intersection_strict")  # DA ∩ ≥2-read (adopted), DA ∩ strict
# per-REPLICATE coverage (coverage_perlib.csv) -> per-strain SD across the 9 libraries (3 tp x 3 reps) for error bars
plf = glob.glob(f"{E16}/*.coverage_perlib.csv")
pl = pd.concat([pd.read_csv(f) for f in plf], ignore_index=True) if plf else None
def sd_strain(col):
    if pl is None or col not in pl.columns: return pd.Series(0.0, index=CANON)
    return pl.groupby("strain")[col].std(ddof=1).reindex(CANON).fillna(0.0)
sdL, sd2, sdS, sdI = sd_strain("cov_loose"), sd_strain("cov_2read"), sd_strain("cov_strict"), sd_strain("cov_intersection")
sdI2, sdIS = sd_strain("cov_intersection_2read"), sd_strain("cov_intersection_strict")
def yerr(ser, sd):  # asymmetric, clipped so the lower cap stays positive on a log axis
    v=ser.values; s=np.nan_to_num(sd.values); return [np.minimum(s, v*0.9), s]
COL = {"loose":"#bdbdbd","2read":"#1565a8","strict":"#7a3b9a","DA":"#E69F00","inter":"#009E73"}
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig = plt.figure(figsize=(13, 10.5), dpi=300)
gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 0.85, 1.0], hspace=0.55)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
x = np.arange(len(CANON)); bw = 0.26
# ---- Panel A: per-strain strain-private coverage, absence ladder (log y) ----
if WPOS: axA.axvspan(min(WPOS)-0.5, max(WPOS)+0.5, color="#C0392B", alpha=0.06, zorder=0)
for j,(lab,ser,key,sd) in enumerate([("absence: loose (≥2/3 reps)",loose,"loose",sdL),("≥2 reads (ADOPTED)",twor,"2read",sd2),("strict (0 reads)",strict,"strict",sdS)]):
    xs=x+(j-1)*bw
    axA.bar(xs, ser.values, bw, color=COL[key], edgecolor="white", linewidth=0.3, label=lab, zorder=3)
    axA.errorbar(xs, ser.values, yerr=yerr(ser,sd), fmt="none", ecolor="#333", elinewidth=0.4, capsize=1.0, capthick=0.4, zorder=5)
    for xi,v,e in zip(xs, ser.values, np.nan_to_num(sd.values)):
        if v>0: axA.text(xi, (v+e)*1.22, f"{v:.2g}", ha="center", va="bottom", fontsize=4.0, rotation=90, color=COL[key], zorder=4)
axA.set_yscale("log"); axA.set_ylim(0.005, 20)
axA.set_xticks(x); labs=axA.set_xticklabels([s.replace("_","/") for s in CANON], rotation=45, ha="right", fontsize=7.5)
for l,s in zip(labs,CANON): l.set_color("#C0392B" if s in WILD else "#333")
axA.set_ylabel("strain-private piRNAs\n(% of total piRNA read-mass, log)", fontsize=8.5)
axA.legend(fontsize=7.5, frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(0,1.005), columnspacing=1.4)
axA.set_title("A  Strain-private piRNA expression as % of each strain's total piRNA read-mass (absence-rule ladder)", fontsize=9.6, fontweight="bold", loc="left", pad=22)
axA.text(np.mean(WPOS), 13, "wild-derived", ha="center", va="top", fontsize=8, fontweight="bold", color="#C0392B")
axA.spines[["top","right"]].set_visible(False); axA.grid(axis="y", alpha=0.25, which="both")
# ---- Panel B: by-group + overall summary ----
def grp(mask, col): d=df[df.strain.isin([s for s in CANON if (s in WILD)==mask])]; return w(d,col)
groups = [("classical", False), ("wild", True), ("all", None)]
cats = [("loose","cov_presence_loose2"),("2read","cov_presence_2read2"),("strict","cov_presence_strict2")]
if HAVE_DA: cats = [("DA","cov_DA")] + cats
gx = np.arange(len(groups)); gbw = 0.8/len(cats)
for j,(key,col) in enumerate(cats):
    vals=[]
    for gname,gm in groups:
        d = df if gm is None else df[df.strain.isin([s for s in CANON if (s in WILD)==gm])]
        vals.append(w(d,col))
    axB.bar(gx+(j-(len(cats)-1)/2)*gbw, vals, gbw, color=COL[key], edgecolor="white", linewidth=0.3,
            label={"loose":"loose","2read":"≥2 reads (adopted)","strict":"strict","DA":"DA-only (enriched, not specific)"}[key], zorder=3)
    for gi,v in zip(gx, vals): axB.text(gi+(j-(len(cats)-1)/2)*gbw, v*1.12, f"{v:.2f}", ha="center", va="bottom", fontsize=5.2, rotation=90, color=COL[key])
axB.set_yscale("log"); axB.set_ylim(0.01, 100 if HAVE_DA else 20)
axB.set_xticks(gx); axB.set_xticklabels(["classical (12)","wild (4)","all 16"], fontsize=8.5)
axB.set_ylabel("% of total piRNA\nread-mass (log)", fontsize=8.5)
axB.legend(fontsize=7, frameon=False, ncol=4, loc="lower left", bbox_to_anchor=(0,1.005))
axB.set_title(f"B  Strain-private coverage by group{' (+ DA-only)' if HAVE_DA else ''} — wild ≈ {grp(True,'cov_presence_2read2'):.1f}% vs classical ≈ {grp(False,'cov_presence_2read2'):.2f}% (≥2-read)", fontsize=9.4, fontweight="bold", loc="left", pad=20)
axB.spines[["top","right"]].set_visible(False); axB.grid(axis="y", alpha=0.25, which="both")
# ---- Panel C: DA-only vs strain-private per strain (the enriched-vs-specific contrast) ----
if HAVE_DA:
    da = per_strain("cov_DA")
    if WPOS: axC.axvspan(min(WPOS)-0.5, max(WPOS)+0.5, color="#C0392B", alpha=0.06, zorder=0)
    seriesC=[("DA-only (strain-ENRICHED, not specific)",da,"DA",None),
             ("DA ∩ loose absence (current pipeline)",inter,"inter",sdI),
             ("DA ∩ ≥2-read absence (ADOPTED)",inter2,"2read",sdI2),
             ("DA ∩ strict absence (0 reads)",interS,"strict",sdIS)]
    nbC=len(seriesC); bwC=0.82/nbC
    for j,(lab,ser,key,sd) in enumerate(seriesC):
        xs=x+(j-(nbC-1)/2)*bwC
        axC.bar(xs, ser.values, bwC, color=COL[key], edgecolor="white", linewidth=0.3, label=lab, zorder=3)
        if sd is not None: axC.errorbar(xs, ser.values, yerr=yerr(ser,sd), fmt="none", ecolor="#333", elinewidth=0.35, capsize=0.7, capthick=0.35, zorder=5)
        ev = np.nan_to_num(sd.values) if sd is not None else np.zeros(len(ser))
        for xi,v,e in zip(xs, ser.values, ev):
            if v>0: axC.text(xi, (v+e)*1.25, f"{v:.2g}", ha="center", va="bottom", fontsize=3.2, rotation=90, color=COL[key], zorder=4)
    axC.set_yscale("log"); axC.set_ylim(0.005, 350)
    axC.set_xticks(x); labs=axC.set_xticklabels([s.replace("_","/") for s in CANON], rotation=45, ha="right", fontsize=7.5)
    for l,s in zip(labs,CANON): l.set_color("#C0392B" if s in WILD else "#333")
    axC.set_ylabel("% of total piRNA\nread-mass (log)", fontsize=8.5)
    axC.legend(fontsize=6.6, frameon=False, ncol=4, loc="lower left", bbox_to_anchor=(0,1.005), columnspacing=1.0)
    axC.set_title("C  DA-only (≈50%, enriched not specific) vs the DA-intersection under the absence ladder (loose → ≥2-read ADOPTED → strict) — the specific set shrinks as absence tightens", fontsize=8.6, fontweight="bold", loc="left", pad=22)
    axC.spines[["top","right"]].set_visible(False); axC.grid(axis="y", alpha=0.25, which="both")
else:
    axC.text(0.5,0.5,"Panel C (DA-only vs strain-private coverage) pending —\nedger16_coverage.R (full tagwise edgeR) still computing; re-run this script when *.coverage_full.csv is ready.",
             ha="center", va="center", transform=axC.transAxes, fontsize=9, color="#999", style="italic")
    axC.axis("off")
fig.suptitle("Expression coverage of strain-private piRNAs across 16 strains (% of total piRNA read-mass; ≥2-read absence adopted)", fontsize=10.6, fontweight="bold", y=1.0)
fig.text(0.5,-0.01,"read-mass-weighted over 3 timepoints · total piRNA expression = all reads in a strain's libraries over every sequence · strain-private ≈ presence/absence-only (DA filter redundant) · "
    f"{'DA-only included' if HAVE_DA else 'DA-only panel pending the ~2h edger16_coverage run'} · red = wild-derived", ha="center", fontsize=6, color="#666")
fig.tight_layout(rect=[0,0,1,0.985])
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_strain_specific_coverage16.{e}", bbox_inches="tight")
# source data
out = pd.DataFrame({"strain":CANON,"cov_loose_%":loose.values,"cov_2read_%":twor.values,"cov_strict_%":strict.values})
if HAVE_DA: out["cov_DA_%"]=per_strain("cov_DA").values
out.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data/Fig_strain_specific_coverage16.csv", index=False)
print("wrote Fig_strain_specific_coverage16 (DA-only included:", HAVE_DA, ")")
print(out.to_string(index=False))
