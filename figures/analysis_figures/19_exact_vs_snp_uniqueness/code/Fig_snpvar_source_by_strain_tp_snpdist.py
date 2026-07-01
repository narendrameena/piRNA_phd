#!/usr/bin/env python3
"""THEME 19 — SOURCE of the 4,394 EXACT-'unique' SNP-variant piRNAs: which STRAINS, which TIMEPOINTS, and how many
SNPs from the nearest expressed allele. They are 1-3 SNP variants (88% a single SNP) of alleles expressed in other
strains, concentrated in the wild-derived strains (SPRET/CAST/PWK/WSB) and peaking at P12.5.
Panels: A counts per strain x timepoint (stacked); B SNP-distance to the nearest allele (donut, 1/2/3 SNP);
C SNP-distance x timepoint (log). Reproducible from data/exact_stagepeak_classified.csv.gz +
unique16/snp_variant_refinement_withintp.csv (SNP-distance cached to SourceData on first run)."""
import warnings, os; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#2ca25f"}
DCOL={1:"#7a0030",2:"#d6336c",3:"#f7a8c4"}; WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
d=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz"); sv=d[d.was_snp_variant].copy()
# --- SNP-distance (mm) per candidate: cache on first run (else read the 1.1GB refinement) ---
cache=f"{T}/data/source_data/SourceData_Fig_snpvar_source_by_strain_tp_snpdist.csv.gz"
if os.path.exists(cache):
    cd=pd.read_csv(cache)
else:
    snpv=set(sv.cand_id); keep=[]
    for ch in pd.read_csv(f"{U}/unique16/snp_variant_refinement_withintp.csv",usecols=["cand_id","mm"],chunksize=1_000_000):
        keep.append(ch[ch.cand_id.isin(snpv)])
    best=pd.concat(keep).sort_values("mm").drop_duplicates("cand_id")
    cd=sv[["cand_id","strain","timepoint"]].merge(best,on="cand_id").rename(columns={"mm":"snp_dist"})
    cd.to_csv(cache,index=False)
N=len(sv)
fig=plt.figure(figsize=(13,6.6),dpi=300); gs=fig.add_gridspec(2,2,width_ratios=[1.35,1],hspace=0.55,wspace=0.32)
axA=fig.add_subplot(gs[:,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,1])
# A: strain x timepoint (stacked horizontal, sorted by total)
piv=sv.groupby(["strain","timepoint"]).size().unstack(fill_value=0).reindex(columns=TPS,fill_value=0)
piv=piv.loc[piv.sum(1).sort_values().index]; y=np.arange(len(piv)); left=np.zeros(len(piv))
for tp in TPS:
    axA.barh(y,piv[tp],left=left,color=TPCOL[tp],label=TPN[tp],edgecolor="white",lw=0.3); left+=piv[tp].values
for i,(st,tot) in enumerate(zip(piv.index,piv.sum(1))): axA.text(tot+N*0.012,i,f"{int(tot):,}",va="center",fontsize=7,fontweight="bold",color="#333")
axA.set_yticks(y); axA.set_yticklabels([s.replace("_","/") for s in piv.index],fontsize=7.5)
for t,st in zip(axA.get_yticklabels(),piv.index):
    if st in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.text(0.97,0.05,"$\\bigstar$ wild-derived strains",transform=axA.transAxes,ha="right",fontsize=8,color="#C0392B",fontstyle="italic")
axA.set_xlabel("EXACT-unique piRNAs that are 1-3 SNP variants (n)",fontsize=9.5); axA.set_xlim(right=piv.sum(1).max()*1.14)
axA.legend(title="timepoint",fontsize=8,title_fontsize=8,frameon=False,loc="lower right",bbox_to_anchor=(1.0,0.12))
axA.spines[["top","right"]].set_visible(False); axA.set_title("A  Source by strain & timepoint",fontsize=10.5,fontweight="bold",loc="left")
# B: SNP-distance donut
vc=cd.snp_dist.value_counts().reindex([1,2,3],fill_value=0)
w,_=axB.pie(vc.values,colors=[DCOL[k] for k in [1,2,3]],startangle=90,counterclock=False,wedgeprops=dict(width=0.42,edgecolor="white"))
axB.text(0,0,f"{N:,}\npiRNAs",ha="center",va="center",fontsize=11,fontweight="bold")
axB.legend(w,[f"{k} SNP — {int(vc[k]):,} ({100*vc[k]/N:.0f}%)" for k in [1,2,3]],fontsize=7.6,frameon=False,loc="center left",bbox_to_anchor=(1.0,0.5))
axB.set_title("B  SNP-distance to nearest allele",fontsize=10.5,fontweight="bold",loc="left")
# C: SNP-distance x timepoint (log)
xt=np.arange(len(TPS)); bw=0.26
for j,k in enumerate([1,2,3]):
    vals=[((cd.timepoint==tp)&(cd.snp_dist==k)).sum() for tp in TPS]
    axC.bar(xt+(j-1)*bw,vals,bw,color=DCOL[k],label=f"{k} SNP")
    for xi,v in zip(xt,vals): axC.text(xi+(j-1)*bw,v,f"{v:,}",ha="center",va="bottom",fontsize=6,color="#444")
axC.set_yscale("log"); axC.set_ylim(bottom=8); axC.set_xticks(xt); axC.set_xticklabels([TPN[t] for t in TPS],fontsize=9)
axC.set_ylabel("piRNAs (log)",fontsize=9); axC.legend(fontsize=7,frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(0.5,1.18),columnspacing=1.0)
axC.spines[["top","right"]].set_visible(False); axC.set_title("C  SNP-distance x timepoint",fontsize=10.5,fontweight="bold",loc="left")
fig.suptitle(f"Theme 19 · the {N:,} \"EXACT-unique\" piRNAs are 1-3 SNP variants ({100*(cd.snp_dist==1).mean():.0f}% single-SNP) of shared alleles\nconcentrated in wild-derived strains (SPRET/CAST/PWK/WSB) and peaking at P12.5",fontsize=11.5,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.92])
out=f"{T}/figures/Fig_snpvar_source_by_strain_tp_snpdist"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print(f"wrote {out} | n={N}; per-strain top:",dict(sv.strain.value_counts().head(4)),"; snp_dist",dict(cd.snp_dist.value_counts().sort_index()))
