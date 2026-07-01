#!/usr/bin/env python3
"""THEME 19 — WHY do the 4,394 SNP-variant piRNAs concentrate in WILD-derived strains? Data-driven decomposition.
It is NOT sequencing depth (~1.5x wild vs classical) but PHYLOGENETIC DIVERGENCE from the C57BL/6 reference: the
more a strain has diverged, the more SNPs sit in its piRNAs -> exact-matching calls them 'unique' though they are
1-3 SNP variants of shared alleles. SPRET/EiJ (a separate species, M. spretus) is the extreme.
Panels: A divergence dose-response, depth-normalised (SNP-variant piRNAs / million reads; wild >> classical);
        B depth ruled out (library size vs SNP-variant count: depth ~1.5x but count ~132x);
        C detection (more unique piRNAs detected per strain, ~70x) DWARFS propensity (fraction that are SNP-variant, ~1.4x).
Reproducible from data/exact_stagepeak_classified.csv.gz + resources/log_stat/strain_srna_mapping_log_stat.tab.
div_rank = subspecies divergence rank from C57BL/6 (M. spretus 5 > M.m.castaneus/musculus 4 > M.m.domesticus-wild 3 > classical 1)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
GU=["unique: conserved-but-silent","unique: strain-private locus","unique: stage-shifted (heterochronic)"]
WILD=["SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"]
SUB={"SPRET_EiJ":"M. spretus (species)","CAST_EiJ":"M. m. castaneus","PWK_PhJ":"M. m. musculus","WSB_EiJ":"M. m. domesticus"}
DIVRANK={"SPRET_EiJ":5,"CAST_EiJ":4,"PWK_PhJ":4,"WSB_EiJ":3}
WCOL={"SPRET_EiJ":"#7a0177","CAST_EiJ":"#c51b8a","PWK_PhJ":"#f768a1","WSB_EiJ":"#fbb4b9"}; CG="#9aa0a6"
# ---- per-strain table from the classification ----
d=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz")
g=d.groupby("strain").agg(N_exact_unique=("klass_exact",lambda s:s.isin(GU).sum()),n_snpvar=("was_snp_variant","sum")).reset_index()
g["frac_snpvar"]=g.n_snpvar/g.N_exact_unique; g["wild"]=g.strain.isin(WILD); g["div_rank"]=g.strain.map(DIVRANK).fillna(1).astype(int)
# ---- library size from the strain STAR mapping log ----
log=pd.read_csv(f"{ROOT}/resources/log_stat/strain_srna_mapping_log_stat.tab",sep=";")
log["strain"]=log.iloc[:,0].str.split("/").str[-2]
g=g.merge(log.groupby("strain")["Number of input reads"].sum().rename("lib_reads"),on="strain",how="left")
g["snpvar_per_Mread"]=g.n_snpvar/g.lib_reads*1e6
g.to_csv(f"{T}/data/source_data/SourceData_Fig_why_wild_divergence_not_depth.csv",index=False)
wd=g[g.wild]; cl=g[~g.wild]
r82=wd.snpvar_per_Mread.mean()/cl.snpvar_per_Mread.mean(); depthx=wd.lib_reads.mean()/cl.lib_reads.mean()
countx=wd.n_snpvar.mean()/cl.n_snpvar.mean(); detx=wd.N_exact_unique.mean()/cl.N_exact_unique.mean(); propx=wd.frac_snpvar.mean()/cl.frac_snpvar.mean()
rho_dr,p_dr=spearmanr(g.div_rank,g.n_snpvar); gg=g.dropna(subset=["lib_reads"]); rho_b,_=spearmanr(gg.lib_reads,gg.n_snpvar)
fig,(axA,axB,axC)=plt.subplots(1,3,figsize=(16.5,5.2),dpi=300,gridspec_kw=dict(width_ratios=[1.25,1,1.1],wspace=0.32))
# A: depth-normalised dose-response
gA=g.sort_values(["div_rank","snpvar_per_Mread"],ascending=False)
cols=[WCOL.get(s,CG) for s in gA.strain]; xb=np.arange(len(gA))
axA.bar(xb,gA.snpvar_per_Mread,color=cols,edgecolor="white",lw=0.4)
axA.set_yscale("log"); axA.set_xticks(xb); axA.set_xticklabels([s.replace("_","/") for s in gA.strain],rotation=55,ha="right",fontsize=6.6)
for t,s in zip(axA.get_xticklabels(),gA.strain):
    if s in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.set_ylabel("SNP-variant piRNAs / million reads (log)",fontsize=9); axA.spines[["top","right"]].set_visible(False)
key="Wild-derived (by divergence):\n"+"\n".join(f"{s.replace('_','/')} — {SUB[s]}" for s in WILD)
axA.text(0.975,0.97,key,transform=axA.transAxes,ha="right",va="top",fontsize=6.6,color="#7a0177",
         bbox=dict(boxstyle="round",fc="white",ec="#c51b8a",alpha=0.9,lw=0.6))
axA.set_title(f"A  Divergence dose-response (depth-normalised)\nstill {r82:.0f}x wild vs classical; rho(rank,count)={rho_dr:.2f} p={p_dr:.0e}",fontsize=9.5,fontweight="bold",loc="left")
# B: depth ruled out
axB.scatter(cl.lib_reads/1e6,cl.n_snpvar,s=26,color=CG,edgecolor="white",lw=0.4,zorder=2)
for _,r in wd.iterrows():
    axB.scatter(r.lib_reads/1e6,r.n_snpvar,s=60,color=WCOL[r.strain],edgecolor="#222",lw=0.6,zorder=4)
    axB.annotate(r.strain.replace("_","/"),(r.lib_reads/1e6,r.n_snpvar),xytext=(5,3),textcoords="offset points",fontsize=7,fontweight="bold",color=WCOL[r.strain])
axB.set_xscale("log"); axB.set_yscale("log"); axB.set_xlabel("library size (million reads, log)",fontsize=9); axB.set_ylabel("SNP-variant piRNAs (log)",fontsize=9)
axB.text(0.04,0.96,f"depth only {depthx:.1f}x wild vs classical\nyet count {countx:.0f}x  (rho={rho_b:.2f})\n-> NOT a depth effect",transform=axB.transAxes,va="top",fontsize=8,fontweight="bold",color="#333")
axB.spines[["top","right"]].set_visible(False); axB.set_title("B  Depth ruled out\n(deeper libraries do not explain the wild excess)",fontsize=9.5,fontweight="bold",loc="left")
# C: detection vs propensity
gC=g.sort_values("N_exact_unique"); yb=np.arange(len(gC)); colsC=[WCOL.get(s,CG) for s in gC.strain]
axC.barh(yb,gC.N_exact_unique,color=colsC,edgecolor="white",lw=0.4)
axC.set_xscale("log"); axC.set_yticks(yb); axC.set_yticklabels([s.replace("_","/") for s in gC.strain],fontsize=6.6)
for t,s in zip(axC.get_yticklabels(),gC.strain):
    if s in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
for i,(s,v) in enumerate(zip(gC.strain,gC.N_exact_unique)): axC.text(v*1.15,i,f"{int(v):,}",va="center",fontsize=6,color="#555")
axC.set_xlabel("exact-unique piRNAs detected per strain (log)",fontsize=9); axC.set_xlim(right=gC.N_exact_unique.max()*2.4); axC.spines[["top","right"]].set_visible(False)
axC.set_title(f"C  Detection {detx:.0f}x  >>  propensity {propx:.1f}x\n(wild strains simply HAVE more unique piRNAs to be SNP-variants)",fontsize=9.5,fontweight="bold",loc="left")
fig.suptitle("WHY wild-derived strains dominate the SNP-variant piRNAs — phylogenetic divergence, not sequencing depth",fontsize=12.5,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.95])
out=f"{T}/figures/Fig_why_wild_divergence_not_depth"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print(f"wrote {out}")
print(f"82x={r82:.1f} depth={depthx:.2f}x count={countx:.0f}x detection={detx:.0f}x propensity={propx:.2f}x rho(rank,count)={rho_dr:.3f}")
print(g[["strain","N_exact_unique","n_snpvar","frac_snpvar","lib_reads","snpvar_per_Mread","div_rank"]].sort_values("n_snpvar",ascending=False).to_string(index=False))
