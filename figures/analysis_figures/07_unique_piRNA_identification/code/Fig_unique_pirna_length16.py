#!/usr/bin/env python3
"""16-strain version of the unique-piRNA length-distribution figure. The genuinely-unique (strain-specific)
piRNAs = final_classified (16-strain). Panels: (A) pooled unique vs bulk reads with the data-driven FWHM window;
(B) per-timepoint (developmental length shift); (C) per-strain (all 16, canonical order); (D) timepoint
composition per length (100% stacked). Distinct sequences per length. Window justified data-drivenly (FWHM)."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
cand=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint","length","klass5"])
cand=cand[cand.klass5.str.startswith("unique")].copy()   # GENUINELY unique only = conserved-but-silent + mm0-clean strain-private (excludes exact, SNP-variant, AND low-quality)
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
cand["tp"]=cand.timepoint.map(TPMAP)
CANON=[s for s in STRAIN_ORDER if s in set(cand.strain)]
PAL=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#9edae5']
SCOL={s:PAL[i] for i,s in enumerate(CANON)}; TCOL={"E16.5":"#0072B2","P12.5":"#009E73","P20.5":"#D55E00"}
LR=range(18,41)
def windows(counts):
    pct=counts/counts.sum()*100; mode=int(pct.idxmax()); half=pct.max()/2
    lo=mode
    while lo-1 in pct.index and pct[lo-1]>=half: lo-=1
    hi=mode
    while hi+1 in pct.index and pct[hi+1]>=half: hi+=1
    return pct,mode,half,lo,hi
pc=cand.drop_duplicates("sequence").groupby("length").size().reindex(LR,fill_value=0)
pct,mode,half,lo,hi=windows(pc[pc>0])
bulk=pd.read_csv(f"{U}/QC_length_distribution.csv").groupby("length")["reads"].sum(); bulk=(bulk/bulk.sum()*100).reindex(LR,fill_value=0)
tp_pct={}; tp_stats={}
for t in TPO:
    c=cand[cand.tp==t].drop_duplicates("sequence").groupby("length").size().reindex(LR,fill_value=0)
    tp_pct[t]=c/c.sum()*100; _,m,_,l,h=windows(c[c>0]); tp_stats[t]=(m,l,h,int(c.sum()))
st_pct={s:(cand[cand.strain==s].drop_duplicates("sequence").groupby("length").size().reindex(LR,fill_value=0).pipe(lambda x:x/x.sum()*100)) for s in CANON}
comp=pd.DataFrame({t:cand[cand.tp==t].groupby("length").size().reindex(LR,fill_value=0) for t in TPO})
comp_frac=comp.div(comp.sum(axis=1).replace(0,np.nan),axis=0)*100
# ---- by genuinely-unique SUBCATEGORY (klass5): conserved-but-silent vs strain-private ----
KORD=["unique: conserved-but-silent","unique: strain-private locus"]
KCOL={"unique: conserved-but-silent":"#0072B2","unique: strain-private locus":"#7a3b9a"}
KLAB={"unique: conserved-but-silent":"conserved-but-silent","unique: strain-private locus":"strain-private locus"}
cu=cand.drop_duplicates("sequence"); kl_pct={}; kl_stats={}
for k in KORD:
    cc=cu[cu.klass5==k].groupby("length").size().reindex(LR,fill_value=0)
    kl_pct[k]=cc/cc.sum()*100; _,m,_,l,h=windows(cc[cc>0]); kl_stats[k]=(m,l,h,int(cc.sum()))
kcomp=pd.DataFrame({k:cu[cu.klass5==k].groupby("length").size().reindex(LR,fill_value=0) for k in KORD})
kcomp_frac=kcomp.div(kcomp.sum(axis=1).replace(0,np.nan),axis=0)*100
print(f"pooled unique: mode {mode} nt, FWHM {lo}-{hi} nt | 24-32 captures {pct[(pct.index>=24)&(pct.index<=32)].sum():.1f}%")
for k in KORD: print(f"  {KLAB[k]}: mode {kl_stats[k][0]} FWHM {kl_stats[k][1]}-{kl_stats[k][2]} n={kl_stats[k][3]:,}")
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(2,2,figsize=(10.4,7.4),dpi=300); xs=np.array(list(LR))
a=ax[0,0]
a.bar(xs,pct.reindex(xs,fill_value=0).values,width=0.82,color=["#E69F00" if lo<=x<=hi else "#d0d0d0" for x in xs],edgecolor="white",linewidth=0.4,zorder=3,label="unique piRNAs")
a.plot(xs,bulk.values,color="#333",lw=1.4,marker="o",ms=2.3,zorder=4,label="bulk reads")
a.axhline(half,ls="--",lw=0.7,color="#999"); a.axvspan(23.5,32.5,color="#0072B2",alpha=0.05,zorder=0)
a.set_title(f"A  Pooled unique piRNA length (mode {mode}, FWHM {lo}-{hi} nt)",fontsize=8.6,fontweight="bold",loc="left")
a.set_xlabel("piRNA length (nt)",fontsize=8.5); a.set_ylabel("% of distinct unique seqs / reads",fontsize=8); a.set_xticks(xs[::2]); a.tick_params(labelsize=6.3); a.legend(fontsize=6.4,frameon=False,loc="upper left")
b=ax[0,1]
for t in TPO:
    b.fill_between(xs,tp_pct[t].reindex(xs,fill_value=0).values,color=TCOL[t],alpha=0.18,zorder=2)
    b.plot(xs,tp_pct[t].reindex(xs,fill_value=0).values,color=TCOL[t],lw=1.8,marker="o",ms=2.6,zorder=3,label=f"{t} (mode {tp_stats[t][0]}, FWHM {tp_stats[t][1]}-{tp_stats[t][2]})")
b.axvspan(23.5,32.5,color="#888",alpha=0.05,zorder=0); b.set_title("B  Unique piRNA length by timepoint (developmental shift)",fontsize=8.6,fontweight="bold",loc="left")
b.set_xlabel("piRNA length (nt)",fontsize=8.5); b.set_ylabel("% of that timepoint's unique seqs",fontsize=8); b.set_xticks(xs[::2]); b.tick_params(labelsize=6.3); b.legend(fontsize=6.1,frameon=False,loc="upper right")
c=ax[1,0]
for k in KORD:
    c.fill_between(xs,kl_pct[k].reindex(xs,fill_value=0).values,color=KCOL[k],alpha=0.16,zorder=2)
    c.plot(xs,kl_pct[k].reindex(xs,fill_value=0).values,color=KCOL[k],lw=1.9,marker="o",ms=2.6,zorder=3,label=f"{KLAB[k]} (mode {kl_stats[k][0]}, FWHM {kl_stats[k][1]}-{kl_stats[k][2]}; n={kl_stats[k][3]:,})")
c.axvspan(23.5,32.5,color="#888",alpha=0.05,zorder=0); c.set_title("C  Unique piRNA length by genuinely-unique subcategory (CBS vs strain-private)",fontsize=8.6,fontweight="bold",loc="left")
c.set_xlabel("piRNA length (nt)",fontsize=8.5); c.set_ylabel("% of that subcategory's unique seqs",fontsize=8); c.set_xticks(xs[::2]); c.tick_params(labelsize=6.3); c.legend(fontsize=6.0,frameon=False,loc="upper right")
d=ax[1,1]; mask=(comp_frac.index>=23)&(comp_frac.index<=33); xd=comp_frac.index[mask].values
d.stackplot(xd,*[comp_frac.loc[mask,t].values for t in TPO],colors=[TCOL[t] for t in TPO],labels=TPO,alpha=0.9,edgecolor="white",linewidth=0.3)
d.set_ylim(0,100); d.set_xlim(xd.min(),xd.max()); d.set_title("D  Timepoint composition at each length (100% stacked)",fontsize=8.6,fontweight="bold",loc="left")
d.set_xlabel("piRNA length (nt)",fontsize=8.5); d.set_ylabel("% of unique seqs of that length",fontsize=8); d.set_xticks(xd[::1]); d.tick_params(labelsize=6.3); d.legend(fontsize=6.4,frameon=False,loc="lower center",ncol=3)
fig.suptitle("Unique (strain-specific) piRNA length distribution across 16 strains — pooled, by timepoint, and by genuinely-unique subcategory (conserved-but-silent vs strain-private)",fontsize=9.2,fontweight="bold",y=1.005)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_unique_pirna_length16.{e}",bbox_inches="tight")
out=pd.DataFrame({"length":list(LR)}); out["pooled_unique_pct"]=pct.reindex(LR,fill_value=0).round(3).values; out["bulk_reads_pct"]=bulk.reindex(LR,fill_value=0).round(3).values
for t in TPO: out[f"{t}_unique_pct"]=tp_pct[t].reindex(LR,fill_value=0).round(3).values
for s in CANON: out[f"{s}_unique_pct"]=st_pct[s].reindex(LR,fill_value=0).round(3).values
for k in KORD: out[f"{KLAB[k]}_pct"]=kl_pct[k].reindex(LR,fill_value=0).round(3).values
out.to_csv(f"{U}/pangenome_te/SourceData_unique_pirna_length16.csv",index=False)
print("wrote Fig_unique_pirna_length16.{png,pdf,svg} + source data")
