#!/usr/bin/env python3
"""Locus-specific SCHEMATIC #2 — a strain-private IAP (LTR/ERVK endogenous retrovirus) insertion that
spawned a strain-specific piRNA cluster. Real example = SPRET_EiJ chr7:3,079,492-3,083,591, IAPLTR3-int
(LTR/ERVK), 80 strain-private piRNAs, ALL antisense to the IAP (100% silencing-competent). Multi-strain
genome-track cartoon: SPRET carries the IAP + the antisense piRNA cluster; the 15 other strains have the
orthologous locus with no insertion and no piRNAs."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
CHROM,S,E="chr7",3079492,3083591; TE="IAPLTR3-int (LTR/ERVK — IAP)"; NPI=80
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(9.4,6.0),dpi=300); ax.set_xlim(0,10); ax.set_ylim(0,11.2); ax.axis("off")
L,R=1.2,8.8; ins0,ins1=4.1,5.9
def genome_line(y,gap=False):
    if gap:
        ax.plot([L,ins0],[y,y],color="#333",lw=2.2); ax.plot([ins1,R],[y,y],color="#333",lw=2.2)
        ax.plot([ins0,ins1],[y,y],color="#bbb",lw=1.0,ls=(0,(2,2)))
    else:
        ax.plot([L,R],[y,y],color="#333",lw=2.2)
for name,y,has in [("SPRET/EiJ",8.0,True),("C57BL/6NJ",6.2,False),("CAST/EiJ",4.8,False),("+ 13 other strains",3.4,False)]:
    genome_line(y,gap=not has)
    ax.text(L-0.15,y,name,ha="right",va="center",fontsize=9,fontweight="bold" if has else "normal",
            color="#C0392B" if has else "#555")
    if has:
        # IAP box (+ strand)
        ax.add_patch(Rectangle((ins0,y-0.22),ins1-ins0,0.44,facecolor="#6a3d9a",edgecolor="#3d1f5c",lw=1.0,zorder=3,alpha=0.85))
        ax.annotate("",xy=(ins1-0.12,y),xytext=(ins0+0.12,y),arrowprops=dict(arrowstyle="-|>",color="white",lw=1.4))
        ax.text((ins0+ins1)/2,y-0.45,f"{TE}  —  SPRET-private insertion",ha="center",va="top",fontsize=7.5,color="#3d1f5c")
        # piRNA cluster: ALL antisense (red, pointing opposite the + IAP) + peak
        xs=np.linspace(ins0-0.05,ins1+0.05,60); peak=0.7*np.exp(-((xs-(ins0+ins1)/2)/0.55)**2)
        ax.fill_between(xs,y+0.35,y+0.35+peak,color="#C0392B",alpha=0.28,zorder=2)
        for xx in np.linspace(ins0+0.12,ins1-0.12,8):
            ax.annotate("",xy=(xx-0.30,y+0.55),xytext=(xx,y+0.55),arrowprops=dict(arrowstyle="-|>",color="#C0392B",lw=1.3))
        ax.text((ins0+ins1)/2,y+1.28,f"piRNA cluster — {NPI} strain-private piRNAs\n(100% ANTISENSE to the IAP → silencing-competent)",
                ha="center",va="bottom",fontsize=7.6,color="#C0392B",fontweight="bold")
    else:
        ax.text((ins0+ins1)/2,y+0.28,"no insertion · no piRNAs",ha="center",va="bottom",fontsize=6.8,color="#888",style="italic")
for xx in (ins0,ins1): ax.plot([xx,xx],[3.2,8.2],color="#cccccc",lw=0.7,ls=(0,(1,3)),zorder=0)
ax.text((L+ins0)/2,8.45,"shared flank",ha="center",fontsize=6.5,color="#999"); ax.text((ins1+R)/2,8.45,"shared flank",ha="center",fontsize=6.5,color="#999")
ax.text(5,10.9,"A strain-private IAP (endogenous retrovirus) insertion creates a strain-specific piRNA cluster",ha="center",fontsize=10.5,fontweight="bold")
ax.text(5,10.5,f"real example: SPRET/EiJ {CHROM}:{S:,}-{E:,}  (IAP = the most active mouse LTR retrotransposon)",ha="center",fontsize=8,color="#555")
ax.text(5,1.9,"Mechanism: an IAP endogenous retrovirus inserted in the SPRET lineage (absent at the orthologous locus in the 15 other strains)\n"
              "becomes a piRNA source — processed into piRNAs that are ENTIRELY ANTISENSE to the IAP, poised to silence the active element.\n"
              "A textbook TE-vs-piRNA arms-race event, captured strain-specifically: TE-driven birth of a new piRNA cluster.",
        ha="center",va="top",fontsize=6.8,color="#444")
ax.annotate("",xy=(3.4,2.7),xytext=(3.8,2.7),arrowprops=dict(arrowstyle="-|>",color="#C0392B",lw=1.3))
ax.text(3.9,2.7,"antisense piRNA (silences the IAP)",va="center",fontsize=6.8,color="#C0392B")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_locus_example_IAP.{e}",bbox_inches="tight")
print(f"wrote Fig_locus_example_IAP.{{png,pdf,svg}} for {CHROM}:{S}-{E} {TE} ({NPI} all-antisense piRNAs)")
