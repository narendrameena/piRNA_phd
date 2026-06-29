#!/usr/bin/env python3
"""TEST FIGURE — justifies the ADOPTED >=2-read absence threshold for strain-private piRNA calling.
Across the absence ladder (0-tolerance/strict -> >=2-read adopted -> loose), at the (moot) >=2/3-rep
presence rule, shows (A) pooled strain-private count and (B) per-(strain x timepoint) read-mass coverage
(% of that strain's total piRNA expression, from coverage_probe.R = 100*reads_of_set/total).
Message: 0-tolerance discards real signal (a single index-hop read elsewhere wrongly disqualifies a locus);
loose admits piRNAs expressed elsewhere; >=2-read removes only the single-read noise floor.
Source: edger16/{tp}.coverage_probe.csv  (coverage_probe.R)."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
df=pd.concat([pd.read_csv(f"{U}/{tp}.coverage_probe.csv") for tp in ["16.5dpc","12.5dpp","20.5dpp"]],ignore_index=True)
# absence ladder at fixed >=2/3-rep presence (presence rep-rule is moot: strict1==strict2, 2read1==2read2)
LAD=[("strict2","0-tolerance\n(0 reads elsewhere)","#7a3b9a"),
     ("2read2","≥2-read  (ADOPTED)\n(tolerate 1 stray read)","#009E73"),
     ("loose2","loose\n(<2/3 reps elsewhere)","#9e9e9e")]
keys=[k for k,_,_ in LAD]; labs=[l for _,l,_ in LAD]; cols=[c for _,_,c in LAD]
counts=[int(df[f"n_presence_{k}"].sum()) for k in keys]
covs=[df[f"cov_presence_{k}"].values for k in keys]          # per strain x tp, already % of read mass
means=[float(np.mean(c)) for c in covs]
x=np.arange(3); rng=np.random.default_rng(0)
fig,(axA,axB)=plt.subplots(1,2,figsize=(9.6,4.4),dpi=300,gridspec_kw={"width_ratios":[1,1.15]})

# Panel A — pooled strain-private count
axA.bar(x,counts,width=0.66,color=cols,edgecolor="white",linewidth=0.5,zorder=3)
for xi,v in zip(x,counts): axA.text(xi,v*1.01,f"{v:,}",ha="center",va="bottom",fontsize=8.5,fontweight="bold")
axA.set_xticks(x); axA.set_xticklabels(labs,fontsize=8); axA.set_ylim(0,max(counts)*1.15)
axA.set_ylabel("strain-private piRNAs (Σ over 16 strains × 3 tp)",fontsize=9)
axA.set_title("A  candidate count vs absence threshold",fontsize=9.5,fontweight="bold",loc="left")
axA.spines[["top","right"]].set_visible(False)
# delta annotations (what each step removes)
axA.annotate(f"+{counts[2]-counts[1]:,}\n(low-read noise)",xy=(1.5,(counts[1]+counts[2])/2),ha="center",va="center",fontsize=6.6,color="#555")
axA.annotate(f"−{counts[1]-counts[0]:,}\n(real signal)",xy=(0.5,(counts[0]+counts[1])/2),ha="center",va="center",fontsize=6.6,color="#b2182b")

# Panel B — per (strain x tp) read-mass coverage (% of total piRNA expression)
for xi,c,col in zip(x,covs,cols):
    axB.scatter(np.full(len(c),xi)+rng.uniform(-0.12,0.12,len(c)),c,s=14,color=col,alpha=0.6,edgecolors="white",linewidths=0.3,zorder=3)
axB.plot(x,means,"-o",color="#222",lw=1.4,ms=5,zorder=4)
for xi,m in zip(x,means):
    axB.annotate(f"mean {m:.2f}%",xy=(xi,m),xytext=(xi,m+2.4),ha="center",va="bottom",fontsize=7.5,fontweight="bold",color="#222",zorder=6,
                 bbox=dict(boxstyle="round,pad=0.18",fc="white",ec="#222",lw=0.4),
                 arrowprops=dict(arrowstyle="-",lw=0.5,color="#222"))
axB.set_xticks(x); axB.set_xticklabels(labs,fontsize=8); axB.set_xlim(-0.5,2.7)
axB.set_ylabel("strain-private read-mass coverage\n(% of strain's total piRNA expression)",fontsize=9)
axB.set_title("B  read-mass captured vs absence threshold (each point = strain × timepoint)",fontsize=9.5,fontweight="bold",loc="left")
axB.spines[["top","right"]].set_visible(False)

fig.suptitle("Choosing the absence threshold for strain-private piRNAs — the adopted ≥2-read rule",fontsize=11,fontweight="bold",y=1.0)
fig.text(0.5,-0.05,"0-tolerance loses real strain-private piRNAs (a single index-hopping read elsewhere disqualifies a well-expressed locus): "
  f"−{counts[1]-counts[0]:,} candidates / −{means[1]-means[0]:.2f} pp read-mass vs ≥2-read. loose admits piRNAs expressed elsewhere "
  f"(+{counts[2]-counts[1]:,} candidates for only +{means[2]-means[1]:.2f} pp read-mass = low-expression strays). "
  "≥2-read = principled middle (tolerate the 1-read noise floor). Presence rep-rule is moot (filterByExpr forces all 3 reps).",
  ha="center",fontsize=6.6,color="#555")
fig.tight_layout(rect=[0,0.02,1,0.99])
base=f"{OUT}/Fig_absence_threshold_choice"
for e in ("pdf","svg","png"): fig.savefig(f"{base}.{e}",bbox_inches="tight")
# source data
src=pd.DataFrame({"absence_threshold":["strict_0tolerance","2read_ADOPTED","loose"],
                  "sum_strain_private":counts,"mean_readmass_coverage_pct":[round(m,4) for m in means],
                  "median_readmass_coverage_pct":[round(float(np.median(c)),4) for c in covs]})
src.to_csv(f"{SD}/SourceData_Fig_absence_threshold_choice.csv",index=False)
print(src.to_string(index=False)); print("wrote",base)
