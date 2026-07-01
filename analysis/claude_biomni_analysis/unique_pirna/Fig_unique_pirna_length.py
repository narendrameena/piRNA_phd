#!/usr/bin/env python3
"""Length distribution of the UNIQUE (strain-specific) piRNAs -- with TIMEPOINT contribution.
Each unique piRNA = one distinct strain-specific sequence; we count DISTINCT sequences per length.
Panels: (A) pooled unique vs bulk read pool, data-driven FWHM window; (B) per-timepoint distributions
(developmental length shift: E16.5 pre-pachytene shorter -> P20.5 pachytene longer); (C) per-strain;
(D) timepoint COMPOSITION at each length (100% stacked) -- which stage dominates short vs long unique
piRNAs. Window justified data-drivenly (FWHM = lengths >= half the modal count; no magic number).
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"

cand=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz", usecols=["strain", "timepoint", "length", "sequence"])   # ADOPTED >=2-read set (klass5 5-class); synced from removed OLD_unrestricted_candidates_snapshot
cand=cand[cand.strain.isin(["C57BL_6NJ", "CAST_EiJ", "SPRET_EiJ"])]   # 3-strain pilot subset (matches the SCOL palette below)
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
cand["tp"]=cand.timepoint.map(TPMAP)
PILOT=[s for s in STRAIN_ORDER if s in set(cand.strain)]
SCOL={"C57BL_6NJ":"#0072B2","CAST_EiJ":"#009E73","SPRET_EiJ":"#D55E00"}
TCOL={"E16.5":"#0072B2","P12.5":"#009E73","P20.5":"#D55E00"}   # Okabe-Ito, CB-safe
LR=range(18,41)

def windows(counts):
    pct=counts/counts.sum()*100; mode=int(pct.idxmax()); half=pct.max()/2
    lo=mode
    while lo-1 in pct.index and pct[lo-1]>=half: lo-=1
    hi=mode
    while hi+1 in pct.index and pct[hi+1]>=half: hi+=1
    return pct,mode,half,lo,hi

# pooled (distinct sequences) + bulk
pc=cand.drop_duplicates("sequence").groupby("length").size().reindex(LR,fill_value=0)
pct,mode,half,lo,hi=windows(pc[pc>0])
bulk=pd.read_csv(f"{U}/QC_length_distribution.csv").groupby("length")["reads"].sum()
bulk=(bulk/bulk.sum()*100).reindex(LR,fill_value=0)
# per timepoint (distinct seqs within each timepoint; strain-specific => 1 strain/seq so no dup within tp)
tp_pct={}; tp_stats={}
for t in TPO:
    c=cand[cand.tp==t].groupby("length").size().reindex(LR,fill_value=0)
    tp_pct[t]=c/c.sum()*100
    _,m,_,l,h=windows(c[c>0]); tp_stats[t]=(m,l,h,int(c.sum()))
# per strain
st_pct={}
for s in PILOT:
    c=cand[cand.strain==s].drop_duplicates("sequence").groupby("length").size().reindex(LR,fill_value=0)
    st_pct[s]=c/c.sum()*100
# timepoint composition per length (counts, then 100% stacked)
comp=pd.DataFrame({t:cand[cand.tp==t].groupby("length").size().reindex(LR,fill_value=0) for t in TPO})
comp_frac=comp.div(comp.sum(axis=1).replace(0,np.nan),axis=0)*100

print(f"pooled unique: mode {mode} nt, FWHM {lo}-{hi} nt | 24-32 captures {pct[(pct.index>=24)&(pct.index<=32)].sum():.1f}%")
print("per timepoint (n unique, mode, FWHM):")
for t in TPO: m,l,h,n=tp_stats[t]; print(f"  {t}: n={n:,} | mode {m} nt | FWHM {l}-{h} nt")

# ---------------- figure 2x2 ----------------
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(2,2,figsize=(9.8,7.1),dpi=300)
xs=np.array(list(LR))
# A pooled vs bulk
a=ax[0,0]
a.bar(xs,pct.reindex(xs,fill_value=0).values,width=0.82,
      color=["#E69F00" if lo<=x<=hi else "#d0d0d0" for x in xs],edgecolor="white",linewidth=0.4,zorder=3,label="unique piRNAs")
a.plot(xs,bulk.values,color="#333",lw=1.4,marker="o",ms=2.3,zorder=4,label="bulk reads")
a.axhline(half,ls="--",lw=0.7,color="#999"); a.axvspan(23.5,32.5,color="#0072B2",alpha=0.05,zorder=0)
a.text(32.4,a.get_ylim()[1]*0.95,"locked\n24-32 nt",ha="right",va="top",fontsize=6,color="#0072B2")
a.set_title(f"A  Pooled unique piRNA length (mode {mode}, FWHM {lo}-{hi} nt)",fontsize=8.6,fontweight="bold",loc="left")
a.set_xlabel("piRNA length (nt)",fontsize=8.5); a.set_ylabel("% of distinct unique seqs / reads",fontsize=8)
a.set_xticks(xs[::2]); a.tick_params(labelsize=6.3); a.legend(fontsize=6.4,frameon=False,loc="upper left")
# B per timepoint (developmental shift) -- filled curves
b=ax[0,1]
for t in TPO:
    b.fill_between(xs,tp_pct[t].reindex(xs,fill_value=0).values,color=TCOL[t],alpha=0.18,zorder=2)
    b.plot(xs,tp_pct[t].reindex(xs,fill_value=0).values,color=TCOL[t],lw=1.8,marker="o",ms=2.6,zorder=3,
           label=f"{t} (mode {tp_stats[t][0]}, FWHM {tp_stats[t][1]}-{tp_stats[t][2]})")
b.axvspan(23.5,32.5,color="#888",alpha=0.05,zorder=0)
b.set_title("B  Unique piRNA length by timepoint (developmental shift)",fontsize=8.6,fontweight="bold",loc="left")
b.set_xlabel("piRNA length (nt)",fontsize=8.5); b.set_ylabel("% of that timepoint's unique seqs",fontsize=8)
b.set_xticks(xs[::2]); b.tick_params(labelsize=6.3); b.legend(fontsize=6.1,frameon=False,loc="upper right")
# C per strain
c=ax[1,0]
for s in PILOT:
    c.plot(xs,st_pct[s].reindex(xs,fill_value=0).values,color=SCOL.get(s,"#888"),lw=1.7,marker="o",ms=2.5,label=s.replace("_","/"))
c.axvspan(23.5,32.5,color="#888",alpha=0.05,zorder=0)
c.set_title("C  Unique piRNA length by strain",fontsize=8.6,fontweight="bold",loc="left")
c.set_xlabel("piRNA length (nt)",fontsize=8.5); c.set_ylabel("% of that strain's unique seqs",fontsize=8)
c.set_xticks(xs[::2]); c.tick_params(labelsize=6.3); c.legend(fontsize=6.5,frameon=False)
# D timepoint composition per length (100% stacked area) -- focus 23-33 where data is dense
d=ax[1,1]
mask=(comp_frac.index>=23)&(comp_frac.index<=33)
xd=comp_frac.index[mask].values
d.stackplot(xd,*[comp_frac.loc[mask,t].values for t in TPO],colors=[TCOL[t] for t in TPO],
            labels=TPO,alpha=0.9,edgecolor="white",linewidth=0.3)
d.set_ylim(0,100); d.set_xlim(xd.min(),xd.max())
d.set_title("D  Timepoint composition at each length (100% stacked)",fontsize=8.6,fontweight="bold",loc="left")
d.set_xlabel("piRNA length (nt)",fontsize=8.5); d.set_ylabel("% of unique seqs of that length",fontsize=8)
d.set_xticks(xd[::1]); d.tick_params(labelsize=6.3); d.legend(fontsize=6.4,frameon=False,loc="lower center",ncol=3)
fig.suptitle("Unique (strain-specific) piRNA length distribution and its timepoint structure",fontsize=10,fontweight="bold",y=1.005)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/Fig_unique_pirna_length.{e}",bbox_inches="tight")

# source data (incl. timepoint breakdown)
out=pd.DataFrame({"length":list(LR)})
out["pooled_unique_pct"]=pct.reindex(LR,fill_value=0).round(3).values
out["bulk_reads_pct"]=bulk.reindex(LR,fill_value=0).round(3).values
for t in TPO: out[f"{t}_unique_pct"]=tp_pct[t].reindex(LR,fill_value=0).round(3).values
for t in TPO: out[f"{t}_composition_pct"]=comp_frac.reindex(LR)[t].round(2).values
for s in PILOT: out[f"{s}_unique_pct"]=st_pct[s].reindex(LR,fill_value=0).round(3).values
out.to_csv(f"{U}/SourceData_unique_pirna_length.csv",index=False)
print("wrote Fig_unique_pirna_length.{png,pdf,svg} + SourceData_unique_pirna_length.csv")
