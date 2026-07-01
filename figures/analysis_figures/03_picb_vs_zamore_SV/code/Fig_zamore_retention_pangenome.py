#!/usr/bin/env python3
"""
Theme 03 (REFRAMED) — pangenome RETENTION of Zamore pachytene piRNA loci across 16
strains + per-strain UCSC-liftOver QC comparison. Nature-Genetics style.

fraction_lifted = fraction of a locus's 100-bp windows projecting into a strain via
the cactus pangenome graph (halLiftover). RETAINED = >=0.5.

A  per-strain mean span retained (bar) + LOCI COUNT lifted atop each bar.
B  fraction_lifted distribution with a BROKEN y-axis (tiny disrupted tail + the
   huge retained bar both visible).
C  QC: per-strain loci lifted, pangenome vs UCSC liftOver (counts atop bars).
"""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":7.5,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"xtick.major.width":0.6,
    "ytick.major.width":0.6,"xtick.major.size":2.5,"ytick.major.size":2.5,
    "pdf.fonttype":42,"svg.fonttype":"none"})

CR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
BLUE="#0072B2"; VERM="#D55E00"; SKY="#56B4E9"; GREY="#9aa0a6"
import matplotlib.colors as mc
def pale(c,f=0.45):
    r,g,b=mc.to_rgb(c); return (1-f*(1-r),1-f*(1-g),1-f*(1-b))   # mix toward white

frac=pd.read_csv(f"{CR}/zamore_fraction_lifted.csv")
qc=pd.read_csv(f"{CR}/SourceData_ucsc_vs_pangenome_lifted.csv")
N=frac.locus.nunique()
order=[s for s in STRAIN_ORDER if s in set(frac.strain)]
qc=qc.set_index("strain").reindex(order).reset_index()

fig=plt.figure(figsize=(9.6,6.8),dpi=300,constrained_layout=True)
gs=fig.add_gridspec(2,2,height_ratios=[1,1],width_ratios=[2.1,1.0])
axA=fig.add_subplot(gs[0,0])
gsB=gs[0,1].subgridspec(2,1,height_ratios=[1,2.4],hspace=0.10)
axBt=fig.add_subplot(gsB[0]); axBb=fig.add_subplot(gsB[1])
axC=fig.add_subplot(gs[1,:])

# ── A: mean span retained (bar) + loci count lifted (atop) ────────────────────
m=frac.groupby("strain")["fraction_lifted"].mean().reindex(order)*100
x=np.arange(len(order)); cols=[VERM if s in WILD else BLUE for s in order]
axA.bar(x, m.values, width=0.74, color=cols, edgecolor="white", linewidth=0.3, zorder=3)
for xi,v,nl in zip(x,m.values,qc.pan_lifted):                       # number atop = loci lifted (count)
    axA.text(xi, v+0.2, f"{int(nl)}", rotation=90, ha="center", va="bottom", fontsize=5.0, color="#222", zorder=4)
axA.set_xticks(x)
labs=axA.set_xticklabels(order, rotation=55, ha="right", fontsize=6.4)
for lab,s in zip(labs,order): lab.set_color("#C0392B" if s in WILD else "#222222")
axA.set_ylim(85,104); axA.set_yticks([85,90,95,100])
axA.set_ylabel("locus span retained (%, mean)\nnumber atop = loci lifted (of 214)", fontsize=7.6)
axA.yaxis.grid(True, lw=0.3, color="#eee", zorder=0); axA.set_axisbelow(True)
axA.legend(handles=[Patch(fc=BLUE,label="classical"),Patch(fc=VERM,label="wild-derived")],
           fontsize=6.4, frameon=False, loc="lower right", bbox_to_anchor=(1.0,1.005), ncol=2)
axA.set_title("A   Zamore piRNA loci retained across 16 strains",
              fontsize=8.5, fontweight="bold", loc="left", pad=16)

# ── B: fraction_lifted distribution with BROKEN y-axis ────────────────────────
bins=np.linspace(0,100,26); centers=(bins[:-1]+bins[1:])/2; w=(bins[1]-bins[0])*0.9
fc=frac[~frac.strain.isin(WILD)].fraction_lifted*100
fw=frac[frac.strain.isin(WILD)].fraction_lifted*100
cc,_=np.histogram(fc,bins=bins); cw,_=np.histogram(fw,bins=bins)
hi=(cc+cw).max()
for ax in (axBt,axBb):
    ax.bar(centers,cc,width=w,color=BLUE,edgecolor="white",linewidth=0.2,zorder=3,label="classical")
    ax.bar(centers,cw,bottom=cc,width=w,color=VERM,edgecolor="white",linewidth=0.2,zorder=3,label="wild-derived")
    ax.axvline(50,ls=(0,(3,2)),lw=0.7,color="#333",zorder=2)
axBt.legend(fontsize=5.4,frameon=False,loc="upper left",handlelength=1.0)
axBt.set_ylim(hi*0.94, hi*1.05); axBt.set_yticks([int(round(hi/100)*100)])
axBb.set_ylim(0, 135)
# break styling
axBt.spines['bottom'].set_visible(False); axBb.spines['top'].set_visible(False)
axBt.tick_params(bottom=False,labelbottom=False)
d=.015
kw=dict(transform=axBt.transAxes,color='k',clip_on=False,lw=0.9)
axBt.plot((-d,+d),(-d*3,+d*3),**kw); axBt.plot((1-d,1+d),(-d*3,+d*3),**kw)
kw=dict(transform=axBb.transAxes,color='k',clip_on=False,lw=0.9)
axBb.plot((-d,+d),(1-d*1.2,1+d*1.2),**kw); axBb.plot((1-d,1+d),(1-d*1.2,1+d*1.2),**kw)
axBb.annotate("disrupted\n<50%",xy=(25,70),fontsize=6.0,color="#333",ha="center")
retained=(frac.fraction_lifted>=0.5).mean()*100
axBb.set_xlabel("locus span retained (%)", fontsize=7.5)
axBb.set_ylabel("locus × strain count", fontsize=7.5); axBb.yaxis.set_label_coords(-0.18,0.72)
axBt.set_title(f"B   {retained:.1f}% retained\n(broken axis)", fontsize=8.3, fontweight="bold", loc="left", pad=4)

# ── C: per-strain QC — pangenome vs UCSC loci lifted (counts atop) ─────────────
xc=np.arange(len(order)); bw=0.40
base_col=[VERM if s in WILD else BLUE for s in order]    # hue = wild/classical
pale_col=[pale(c) for c in base_col]
axC.bar(xc-bw/2, qc.pan_pct, width=bw, color=base_col, edgecolor="white", linewidth=0.3, zorder=3)   # pangenome = SOLID
axC.bar(xc+bw/2, qc.ucsc_pct, width=bw, color=pale_col, edgecolor="white", linewidth=0.3, zorder=3)   # UCSC = PALE
for xi,n,pct in zip(xc-bw/2, qc.pan_lifted, qc.pan_pct):
    axC.text(xi, pct+1.2, str(int(n)), rotation=90, ha="center", va="bottom", fontsize=4.8, color="#222", zorder=5, fontweight="bold")
for xi,n,pct in zip(xc+bw/2, qc.ucsc_lifted, qc.ucsc_pct):
    axC.text(xi, pct+1.2, str(int(n)), rotation=90, ha="center", va="bottom", fontsize=4.8, color="#777", zorder=5, fontweight="bold")
axC.set_xticks(xc)
labs=axC.set_xticklabels(order, rotation=55, ha="right", fontsize=6.4)
for lab,s in zip(labs,order): lab.set_color("#C0392B" if s in WILD else "#222222")
axC.set_ylim(0,116); axC.set_yticks([0,25,50,75,100])
axC.set_ylabel(f"loci lifted (% of {N})", fontsize=8)
axC.yaxis.grid(True, lw=0.3, color="#eee", zorder=0); axC.set_axisbelow(True)
axC.legend(handles=[Patch(fc=BLUE,label="classical"),Patch(fc=VERM,label="wild-derived"),
                    Patch(fc="#555",label="pangenome (solid)"),Patch(fc=pale("#555"),label="UCSC (pale)")],
           fontsize=6.2, frameon=False, loc="lower right", bbox_to_anchor=(1.0,1.005), ncol=4, handlelength=1.0, columnspacing=1.0)
panA=100*qc.pan_lifted.sum()/(N*len(order)); ucA=100*qc.ucsc_lifted.sum()/(N*len(order))
axC.set_title(f"C   QC: per-strain loci lifted — pangenome {panA:.0f}% vs UCSC {ucA:.0f}% (n atop bar, of {N}); UCSC under-calls divergent strains",
              fontsize=7.8, fontweight="bold", loc="left", pad=16)

fig.suptitle("Pachytene piRNA loci are structurally retained across 16 strains (pangenome); UCSC liftOver under-calls divergent loci",
             fontsize=9.3, fontweight="bold")
fig.text(0.5,-0.02,
    "fraction_lifted = % of a locus's 100-bp windows projecting into a strain via cactus halLiftover (pangenome); Panel A bar = mean over 214 loci · "
    "UCSC liftOver (≥95% identity) scores sequence-diverged loci as absent (SPRET 52/214) → false 'not present'; the pangenome resolves divergence · "
    "pachytene loci = conserved position, divergent sequence (Yu 2021 PMID 33397987)",
    ha="center", fontsize=5.4, color="#666")
base=f"{CR}/Fig_zamore_retention_pangenome"
for ext in ("pdf","svg","png"): fig.savefig(f"{base}.{ext}", bbox_inches="tight")
print("wrote", base)
print(f"retained {retained:.1f}% | pangenome {panA:.1f}% vs UCSC {ucA:.1f}% | hi-bar={hi}")
