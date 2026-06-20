#!/usr/bin/env python3
"""Locus SCHEMATIC #2 (16-strain) — a strain-private IAP (LTR/ERVK endogenous retrovirus) insertion that
spawned a strain-specific piRNA cluster. Real example = SPRET_EiJ chr7:3,079,492-3,083,591, IAPLTR3-int,
80 strain-private piRNAs, ALL antisense to the IAP (silencing-competent). Top = the detailed SPRET locus
(IAP + its antisense piRNA cluster). Below = the cross-strain presence across ALL 16 strains (canonical
order): SPRET carries the insertion + cluster; the other 15 have the orthologous locus with no insertion and
no piRNAs (verified absent — Step-4 genome-anchored expression test + pangenome halLiftover; private 1/16)."""
import numpy as np
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
CHROM,S,E="chr7",3079492,3083591; TE="IAPLTR3-int (LTR/ERVK — IAP)"; NPI=80
OTHERS=[s for s in STRAIN_ORDER if s not in ("C57BL_6","SPRET_EiJ")]   # 15 other strains, canonical order
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(9.8,9.6),dpi=300); ax.set_xlim(0,10); ax.set_ylim(0,17.4); ax.axis("off")
L,R=2.4,8.8; ins0,ins1=4.6,6.0; ctr=(ins0+ins1)/2
# ---- title ----
ax.text(5,17.1,"A strain-private IAP (endogenous retrovirus) insertion creates a strain-specific piRNA SOURCE LOCUS (individual piRNAs; not a PICB-called cluster)",ha="center",fontsize=9.6,fontweight="bold")
ax.text(5,16.7,f"real example: SPRET/EiJ {CHROM}:{S:,}-{E:,}  (IAP = the most active mouse LTR retrotransposon)",ha="center",fontsize=8,color="#555")
# ---- detailed SPRET locus (top) ----
yD=14.7
ax.plot([L,R],[yD,yD],color="#333",lw=2.4)
ax.text(L-0.18,yD,"SPRET/EiJ",ha="right",va="center",fontsize=9,fontweight="bold",color="#C0392B")
ax.add_patch(Rectangle((ins0,yD-0.22),ins1-ins0,0.44,facecolor="#6a3d9a",edgecolor="#3d1f5c",lw=1.0,zorder=3,alpha=0.85))
ax.annotate("",xy=(ins1-0.12,yD),xytext=(ins0+0.12,yD),arrowprops=dict(arrowstyle="-|>",color="white",lw=1.4))
ax.text(ctr,yD-0.5,f"{TE}  —  SPRET-private insertion",ha="center",va="top",fontsize=7.5,color="#3d1f5c")
xs=np.linspace(ins0-0.05,ins1+0.05,60); peak=0.8*np.exp(-((xs-ctr)/0.55)**2)
ax.fill_between(xs,yD+0.36,yD+0.36+peak,color="#C0392B",alpha=0.28,zorder=2)
for xx in np.linspace(ins0+0.12,ins1-0.12,8): ax.annotate("",xy=(xx-0.30,yD+0.6),xytext=(xx,yD+0.6),arrowprops=dict(arrowstyle="-|>",color="#C0392B",lw=1.3))
ax.text(ctr,yD+1.5,f"piRNA cluster — {NPI} strain-private piRNAs (100% ANTISENSE to the IAP → silencing-competent)",ha="center",va="bottom",fontsize=7.6,color="#C0392B",fontweight="bold")
ax.text((L+ins0)/2,yD+0.28,"shared flank",ha="center",fontsize=6.3,color="#999"); ax.text((ins1+R)/2,yD+0.28,"shared flank",ha="center",fontsize=6.3,color="#999")
# ---- cross-strain presence across all 16 strains ----
ax.text(L-0.18,12.95,"Cross-strain presence — all 16 strains (canonical order):",ha="left",fontsize=8.2,fontweight="bold",color="#333")
rows=["SPRET_EiJ"]+OTHERS; ys=np.linspace(12.3,1.7,len(rows))
for nm,y in zip(rows,ys):
    has=(nm=="SPRET_EiJ")
    ax.text(L-0.18,y,nm.replace("_","/"),ha="right",va="center",fontsize=6.6,fontweight="bold" if has else "normal",color="#C0392B" if has else "#555")
    if has:
        ax.plot([L,R],[y,y],color="#333",lw=1.7)
        ax.add_patch(Rectangle((ins0,y-0.11),ins1-ins0,0.22,fc="#6a3d9a",ec="#3d1f5c",lw=0.6,zorder=3,alpha=0.85))
        ax.text(R+0.12,y,f"IAP + {NPI} antisense piRNAs (present)",ha="left",va="center",fontsize=6.0,color="#C0392B",fontweight="bold")
    else:
        ax.plot([L,ins0],[y,y],color="#333",lw=1.4); ax.plot([ins1,R],[y,y],color="#333",lw=1.4)
        ax.plot([ins0,ins1],[y,y],color="#cccccc",lw=0.9,ls=(0,(2,2)))
ax.text(R+0.12,ys[len(rows)//2],"no insertion · 0 cluster piRNAs (absent in all 15)",ha="left",va="center",fontsize=6.2,color="#888",style="italic")
for xx in (ins0,ins1): ax.plot([xx,xx],[1.4,12.6],color="#e6e6e6",lw=0.6,ls=(0,(1,3)),zorder=0)
ax.text(5,0.75,"Mechanism: an IAP endogenous retrovirus inserted in the SPRET lineage (absent at the orthologous locus in the 15 other strains) becomes a piRNA source — processed into piRNAs\n"
              "ENTIRELY ANTISENSE to the IAP, poised to silence the active element. A textbook TE-vs-piRNA arms-race event captured strain-specifically: TE-driven birth of a new piRNA cluster.\n"
              "Absence in the 15 other strains is data-derived (Step-4 genome-anchored expression test + pangenome halLiftover), not assumed.",ha="center",va="center",fontsize=6.6,color="#444")
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_locus_example_IAP.{e}",bbox_inches="tight")
print(f"wrote Fig_locus_example_IAP.{{png,pdf,svg}} (16-strain cross-strain) for {CHROM}:{S}-{E} {TE}")
