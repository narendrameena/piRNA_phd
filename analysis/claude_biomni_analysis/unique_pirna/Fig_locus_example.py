#!/usr/bin/env python3
"""Locus-specific SCHEMATIC (biology, not numbers): a real strain-private TE insertion that spawned a
strain-specific piRNA cluster. Example = SPRET_EiJ chr2 ~150.45-150.52 Mb, L1MdN_I (LINE/L1), 368
strain-private piRNAs (mostly antisense -> silencing-competent). Drawn as a multi-strain genome-track
cartoon: SPRET carries the insertion + piRNA cluster; the 15 other strains have the orthologous locus
with NO insertion and NO piRNAs. Illustrates TE-driven piRNA-cluster birth."""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, FancyArrowPatch
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
ex=json.load(open(f"{PG}/example_locus.json"))
fam=ex["te_family"].split("|")[0]; cls=ex["te_family"].split("|")[-1]
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(9.4,6.0),dpi=300); ax.set_xlim(0,10); ax.set_ylim(0,11.2); ax.axis("off")
L,R=1.2,8.8; ins0,ins1=4.0,6.0   # flank coords + insertion span (schematic)
def genome_line(y,gap=False):
    if gap:
        ax.plot([L,ins0],[y,y],color="#333",lw=2.2); ax.plot([ins1,R],[y,y],color="#333",lw=2.2)
        ax.plot([ins0,ins1],[y,y],color="#bbb",lw=1.0,ls=(0,(2,2)))   # no insertion (ancestral)
    else:
        ax.plot([L,R],[y,y],color="#333",lw=2.2)
tracks=[("SPRET/EiJ",8.4,True),("C57BL/6NJ",6.2,False),("CAST/EiJ",4.8,False),("+ 13 other strains",3.4,False)]
for name,y,has in tracks:
    genome_line(y,gap=not has)
    ax.text(L-0.15,y,name,ha="right",va="center",fontsize=9,fontweight="bold" if has else "normal",
            color="#D55E00" if has else "#555")
    if has:
        # the L1 insertion box (- strand)
        ax.add_patch(Rectangle((ins0,y-0.22),ins1-ins0,0.44,facecolor="#E69F00",edgecolor="#7a5b00",lw=1.0,zorder=3))
        ax.annotate("",xy=(ins0+0.15,y),xytext=(ins1-0.15,y),arrowprops=dict(arrowstyle="-|>",color="#5a4400",lw=1.4))
        ax.text((ins0+ins1)/2,y-0.45,f"{fam}  ({cls})  — SPRET-private insertion",ha="center",va="top",fontsize=7.5,color="#7a5b00")
        # piRNA cluster above: filled peak + arrows (mostly antisense, red; some sense, grey)
        import numpy as np
        xs=np.linspace(ins0-0.1,ins1+0.1,60); peak=0.7*np.exp(-((xs-(ins0+ins1)/2)/0.7)**2)
        ax.fill_between(xs,y+0.35,y+0.35+peak,color="#C0392B",alpha=0.25,zorder=2)
        for i,xx in enumerate(np.linspace(ins0+0.1,ins1-0.1,9)):
            anti = i not in (3,7)   # mostly antisense
            col="#C0392B" if anti else "#888"; dirn=0.32 if anti else -0.32
            ax.annotate("",xy=(xx+dirn,y+0.55),xytext=(xx,y+0.55),arrowprops=dict(arrowstyle="-|>",color=col,lw=1.3))
        ax.text((ins0+ins1)/2,y+1.25,"piRNA cluster — 368 strain-private piRNAs\n(mostly ANTISENSE to the L1 → silencing-competent)",
                ha="center",va="bottom",fontsize=7.6,color="#C0392B",fontweight="bold")
    else:
        ax.text((ins0+ins1)/2,y+0.28,"no insertion · no piRNAs",ha="center",va="bottom",fontsize=6.8,color="#888",style="italic")
# homology connectors (flanks shared across strains)
for xx in (ins0,ins1):
    ax.plot([xx,xx],[3.2,8.6],color="#cccccc",lw=0.7,ls=(0,(1,3)),zorder=0)
ax.text((L+ins0)/2,8.95,"shared flank",ha="center",fontsize=6.5,color="#999")
ax.text((ins1+R)/2,8.95,"shared flank",ha="center",fontsize=6.5,color="#999")
# title + mechanism caption
ax.text(5,10.9,"A strain-private LINE/L1 insertion creates a strain-specific piRNA cluster",ha="center",fontsize=11,fontweight="bold")
ax.text(5,10.5,f"real example: SPRET/EiJ chr2:{ex['ins_start']:,}-{ex['ins_end']:,}",ha="center",fontsize=8,color="#555")
ax.text(5,1.9,"Mechanism: an L1 (LINE) element inserted in the SPRET lineage (absent at the orthologous locus in the 15 other strains) becomes a piRNA\n"
              "source — the locus is transcribed and processed into piRNAs that are predominantly ANTISENSE to the L1, poised to silence it.\n"
              "This is TE-driven birth of a strain-specific piRNA cluster (the minority, ~800-1200/divergent strain; most strain-private piRNAs are divergence-driven).",
        ha="center",va="top",fontsize=6.8,color="#444")
# legend
ax.annotate("",xy=(2.0,2.7),xytext=(1.6,2.7),arrowprops=dict(arrowstyle="-|>",color="#C0392B",lw=1.3)); ax.text(2.1,2.7,"antisense piRNA (silencing)",va="center",fontsize=6.5,color="#C0392B")
ax.annotate("",xy=(5.6,2.7),xytext=(6.0,2.7),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.3)); ax.text(5.7,2.7,"sense piRNA (ping-pong)",va="center",fontsize=6.5,color="#888")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_locus_example.{e}",bbox_inches="tight")
print(f"wrote Fig_locus_example.{{png,pdf,svg}} for {ex['chrom']}:{ex['ins_start']}-{ex['ins_end']} {ex['te_family']} ({ex['n_pirna']} piRNAs)")
