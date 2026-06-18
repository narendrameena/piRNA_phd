#!/usr/bin/env python3
"""Pangenome piRNA-cluster PAV — message-driven, integrative (16 strains, GRCm39 frame). Three
biological points: (A) the piRNA-cluster pangenome is OPEN — small conserved core, large private fraction
(frequency spectrum); (B) cluster GAIN tracks strain DIVERGENCE — private-cluster bp highest in wild
strains; (C) conservation is STRATIFIED BY DEVELOPMENTAL STAGE — pachytene clusters are core, pre-
pachytene/hybrid are variable (Zamore stages). Ties PAV + divergence + Ozata-2020 stage biology together."""
import os,subprocess,sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
STR16=["WSB_EiJ","CAST_EiJ","BALB_cJ","C57BL_6NJ","A_J","DBA_2J","NOD_ShiLtJ","SPRET_EiJ","AKR_J","C3H_HeJ","CBA_J","PWK_PhJ","NZO_HlLtJ","FVB_NJ","129S1_SvImJ","LP_J"]

cat=pd.read_csv(f"{U}/cluster_PAV_catalogue.csv.gz",dtype={"chrom":str})
spec=cat.groupby("n_strains").bp.sum()/1e6   # Mb present in exactly k strains
core=cat[cat.n_strains==16].bp.sum()/1e6; acc=cat[(cat.n_strains>=2)&(cat.n_strains<=15)].bp.sum()/1e6; priv=cat[cat.n_strains==1].bp.sum()/1e6
tot=core+acc+priv
# private-by-strain via multiinter
merged=[f"{U}/{X}.GRCm39.merged.bed" for X in STR16]
mi=subprocess.run([BT,"multiinter","-names",*STR16,"-i",*merged],capture_output=True,text=True).stdout
from collections import Counter
pbp=Counter()
for ln in mi.splitlines():
    f=ln.split("\t")
    if len(f)<5 or int(f[3])!=1: continue
    pbp[f[4]]+=(int(f[2])-int(f[1]))
zs=pd.read_csv(f"{U}/zamore_loci_stage_conservation.csv").dropna(subset=["mean_nstrains"])

plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(1,3,figsize=(13.5,4.3),dpi=300)
# A: frequency spectrum (openness)
a=ax[0]; ks=list(range(1,17)); ys=[spec.get(k,0) for k in ks]
cols=["#C0392B"]+["#bdbdbd"]*14+["#0072B2"]
a.bar(ks,ys,color=cols,edgecolor="white",linewidth=0.4,zorder=3)
a.set_xlabel("# strains carrying the cluster (1 = private … 16 = core)",fontsize=8.5)
a.set_ylabel("piRNA-cluster sequence (Mb)",fontsize=8.5); a.set_xticks(ks); a.tick_params(labelsize=6.5)
a.set_title("A  An OPEN piRNA-cluster pangenome",fontsize=9.2,fontweight="bold",loc="left")
a.text(0.97,0.95,f"private {priv:.0f} Mb ({100*priv/tot:.0f}%)\naccessory {acc:.0f} Mb ({100*acc/tot:.0f}%)\ncore {core:.0f} Mb ({100*core/tot:.0f}%)",
       transform=a.transAxes,ha="right",va="top",fontsize=7,
       bbox=dict(boxstyle="round",fc="white",ec="#ccc"))
a.spines[['top','right']].set_visible(False)
# B: private by strain, canonical order, wild highlighted (divergence)
b=ax[1]; order=[s for s in STRAIN_ORDER if s in pbp]; vals=[pbp[s]/1e6 for s in order]
colsB=["#D55E00" if s in WILD else "#7fb3d5" for s in order]
b.barh(range(len(order))[::-1],vals,color=colsB,edgecolor="white",height=0.72,zorder=3)
for i,v in enumerate(vals): b.text(v+0.2,len(order)-1-i,f"{v:.1f}",va="center",fontsize=6)
b.set_yticks(range(len(order))[::-1]); b.set_yticklabels([s.replace("_","/") for s in order],fontsize=6.6)
b.set_xlabel("strain-private cluster sequence (Mb)",fontsize=8.5)
b.set_title("B  Cluster gain tracks divergence (wild strains, orange)",fontsize=9.2,fontweight="bold",loc="left")
b.spines[['top','right']].set_visible(False)
# C: conservation by Zamore stage (stage stratification)
c=ax[2]; ORD=["prepachytene","hybrid","pachytene"]; LAB=["pre-\npachytene","hybrid","pachytene"]; CCOL={"prepachytene":"#0072B2","hybrid":"#009E73","pachytene":"#D55E00"}
data=[zs.loc[zs.stage==s,"mean_nstrains"].values for s in ORD]
bp_=c.boxplot(data,patch_artist=True,widths=0.6,showfliers=False,medianprops=dict(color="black",lw=1.1))
for patch,s in zip(bp_["boxes"],ORD): patch.set_facecolor(CCOL[s]); patch.set_alpha(0.75)
for i,s in enumerate(ORD):
    y=zs.loc[zs.stage==s,"mean_nstrains"].values; x=np.random.default_rng(0).normal(i+1,0.05,len(y))
    c.scatter(x,y,s=6,color=CCOL[s],edgecolor="white",linewidth=0.2,alpha=0.55,zorder=3)
    c.text(i+1,16.5,f"n={len(y)}",ha="center",fontsize=6.5,color=CCOL[s])
c.axhline(16,ls=":",lw=0.6,color="#888"); c.set_ylim(4,17.5); c.set_xticks([1,2,3]); c.set_xticklabels(LAB,fontsize=8)
c.set_ylabel("# strains carrying the cluster (/16)",fontsize=8.5)
c.set_title("C  Conservation is stage-stratified (pachytene = core)",fontsize=9.2,fontweight="bold",loc="left")
c.spines[['top','right']].set_visible(False)
fig.suptitle("piRNA-cluster pangenome across 16 mouse strains: open, divergence-driven, and developmentally stratified",fontsize=10.5,fontweight="bold",y=1.02)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/Fig_pangenome_pav.{e}",bbox_inches="tight")
pd.DataFrame({"strain":order,"private_Mb":[round(v,2) for v in vals]}).to_csv(f"{U}/SourceData_private_by_strain.csv",index=False)
print(f"core {core:.1f} / accessory {acc:.1f} / private {priv:.1f} Mb (total {tot:.1f})")
print("private Mb by strain:",{s:round(pbp[s]/1e6,1) for s in order})
print("wrote Fig_pangenome_pav.{png,pdf,svg}")
