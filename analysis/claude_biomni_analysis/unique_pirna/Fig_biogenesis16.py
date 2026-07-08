#!/usr/bin/env python3
"""16-strain quantitative biogenesis signatures (the per-strain version of the single-locus Fig_biogenesis).
(A) ping-pong: the 10-nt 5'-overlap z-score per strain x timepoint (MILI/MIWI2 amplification, pre-pachytene-
biased); (B) phasing: the +1-nt head-to-tail fraction per strain x timepoint (Zucchini/PLD6, pachytene-biased).
Across all 16 strains the two signatures show the expected developmental split — ping-pong relatively stronger
early (E16.5), phasing rising to a pachytene (P20.5) maximum. Data = results/pingpong/all_zscore.tab (offset 10)
+ phasing_allstrains_1random frac_plus1; means over replicates."""
import warnings; warnings.filterwarnings("ignore")
import sys,glob,os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis"; PG=f"{U}/unique_pirna/pangenome_te"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
TPO=["E16.5","P12.5","P20.5"]; TCOL={"E16.5":"#0072B2","P12.5":"#009E73","P20.5":"#D55E00"}
# ping-pong: 10-nt z-score per strain x tp
z=pd.read_csv(f"{ROOT}/results/pingpong/all_zscore.tab",sep="\t",header=None,names=["chrom","offset","count","z","src"])
z=z[z.offset==10].copy(); z["strain"]=z.src.str.split("/").str[0]
z["samp"]=z.src.str.split("/").str[1]; z["tp"]=z.samp.str.rsplit("-",n=1).str[1].str.rsplit(".",n=1).str[0].map(TPMAP)
ppz=z.groupby(["strain","tp"]).z.mean().reset_index()
# phasing: +1 fraction per strain x tp
rows=[]
for f in glob.glob(f"{U}/phasing_allstrains_1random/*.phasing_summary.csv"):
    d=pd.read_csv(f)
    if "direction" in d and (d.direction=="follow").any(): d=d[d.direction=="follow"]
    base=os.path.basename(f).replace(".phasing_summary.csv",""); strain=base.split("-")[0]; rest=base.split("-",1)[1]; tp=rest.rsplit(".",1)[0]
    rows.append(dict(strain=strain,tp=TPMAP.get(tp),frac=float(d.frac_plus1.mean())*100))
ph=pd.DataFrame(rows).groupby(["strain","tp"]).frac.mean().reset_index()
ppP=ppz.pivot(index="strain",columns="tp",values="z").reindex(CANON); phP=ph.pivot(index="strain",columns="tp",values="frac").reindex(CANON)
ppP.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/data/source_data/SourceData_biogenesis16_pingpong.csv"); phP.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/15_pirna_biology_misc/data/source_data/SourceData_biogenesis16_phasing.csv")
print("ping-pong z10 mean by tp:",{t:round(ppP[t].mean(),2) for t in TPO if t in ppP})
print("phasing %  mean by tp:",{t:round(phP[t].mean(),1) for t in TPO if t in phP})
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,(axA,axB)=plt.subplots(2,1,figsize=(13,8.2),dpi=300); x=np.arange(len(CANON)); bw=0.27
for i,t in enumerate(TPO):
    if t in ppP:
        axA.bar(x+(i-1)*bw,ppP[t].values,bw,color=TCOL[t],label=t)
        for xi,v in zip(x+(i-1)*bw,ppP[t].values):
            if pd.notna(v): axA.text(xi,v-0.06,f"{v:.1f}",ha="center",va="top",fontsize=3.8,rotation=90,color="white",fontweight="bold")   # number INSIDE its own bar → separated by bar width, no inter-bar overlap/clipping
    if t in phP:
        axB.bar(x+(i-1)*bw,phP[t].values,bw,color=TCOL[t],label=t)
        for xi,v in zip(x+(i-1)*bw,phP[t].values):
            if pd.notna(v): axB.text(xi,v-1.0,f"{v:.0f}",ha="center",va="top",fontsize=3.8,rotation=90,color="white",fontweight="bold")   # number INSIDE its own bar → separated by bar width, no inter-bar overlap/clipping
axA.axhline(1.96,ls="--",lw=0.8,color="#888"); axA.set_ylim(0,4.95); axB.set_ylim(0,69)
axA.set_xticks(x); axA.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.5)
axA.set_ylabel("ping-pong 10-nt z-score",fontsize=9); axA.legend(fontsize=7.5,frameon=True,facecolor="white",framealpha=0.9,edgecolor="none",title="timepoint",ncol=3,loc="upper left")
axA.set_title("A  Ping-pong (10-nt 5′ overlap) z-score per strain × timepoint — MILI/MIWI2 (dashed = z 1.96, p≈0.05; all strains far above)",fontsize=9.6,fontweight="bold",loc="left")
axA.spines[['top','right']].set_visible(False)
axB.set_xticks(x); axB.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.5)
axB.set_ylabel("+1-nt phasing fraction (%)",fontsize=9); axB.legend(fontsize=7.5,frameon=True,facecolor="white",framealpha=0.9,edgecolor="none",title="timepoint",ncol=3,loc="upper left")
axB.set_title("B  Phasing (+1-nt head-to-tail) fraction per strain × timepoint — Zucchini/PLD6, rises to a pachytene (P20.5) maximum",fontsize=9.8,fontweight="bold",loc="left")
axB.spines[['top','right']].set_visible(False)
fig.suptitle("piRNA biogenesis signatures across all 16 strains: ping-pong vs phasing, by developmental timepoint",fontsize=10.6,fontweight="bold",y=1.0)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_biogenesis16.{e}",bbox_inches="tight")
print("wrote Fig_biogenesis16.{png,pdf,svg} + source data")
