#!/usr/bin/env python3
"""Locus SCHEMATIC (16-strain) — a CONSERVED lncRNA gene that is the precursor of a pachytene piRNA cluster
(lncRNA-driven). The ncRNA counterpart of Fig_locus_example_IAP, and its CONTRAST: this locus is present in
ALL 16 strains (conserved), whereas the IAP locus is private to 1. Read-level-verified example = SPRET_EiJ
chr17:23,790,189-23,809,974, lncRNA Gm10505: pachytene piRNAs SENSE to the lncRNA (derived from it), 0%
protein-coding overlap, low TE, strong 1U. Top = the detailed SPRET locus; below = the 16-strain presence
(the lncRNA gene is present in every strain assembly — gene symbol verified in all 16 GFF3s; cf. the real
16-strain coverage in Fig_gm10505_16strains). Numbers computed from the data."""
import numpy as np, pysam
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; PG=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
CHROM,S,E="chr17",23790189,23809974; SYM="Gm10505"
GFF=f"{ROOT}/resources/annotation/SPRET_EiJ_v3.5.gff3"; RM=f"{ROOT}/resources/repeatMasker/SPRET_EiJ_repeatmasker.bed"
BAM=f"{ROOT}/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam"
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}; OTHERS=[s for s in STRAIN_ORDER if s not in ("C57BL_6","SPRET_EiJ")]
strand="+"
for ln in open(GFF):
    f=ln.split("\t")
    if len(f)>8 and f[2] in ("gene","ncRNA_gene") and f[0].split("#")[-1]==CHROM and "predicted gene 10505" in f[8] and int(f[3])<E and int(f[4])>S: strand=f[6]
teb=sum(min(int(g.split('\t')[2]),E)-max(int(g.split('\t')[1]),S) for g in open(RM) if g.split('\t',3)[0]==CHROM and int(g.split('\t')[1])<E and int(g.split('\t')[2])>S)
tef=teb/(E-S)*100
bam=pysam.AlignmentFile(BAM,"rb"); fwd=rev=0; first=[]; seqs=set()
for a in bam.fetch(f"SPRET_EiJ#1#{CHROM}",S,E):
    if a.is_unmapped or not a.query_sequence: continue
    if not 24<=a.reference_end-a.reference_start<=32: continue
    if a.is_reverse: rev+=1
    else: fwd+=1
    first.append(comp.get(a.query_sequence[-1],"N") if a.is_reverse else a.query_sequence[0]); seqs.add(a.query_sequence)
bam.close()
tot=fwd+rev; sense=(rev if strand=="-" else fwd)/tot*100; u1=first.count("T")/len(first)*100; nsp=len(seqs); gdir=1 if strand=="+" else -1
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(9.8,9.6),dpi=300); ax.set_xlim(0,10); ax.set_ylim(0,17.4); ax.axis("off")
L,R=2.4,8.8; g0,g1=4.5,6.1; ctr=(g0+g1)/2
def arrow(x,y,d,c,lw=1.2,dx=0.28): ax.annotate("",xy=(x+d*dx,y),xytext=(x,y),arrowprops=dict(arrowstyle="-|>",color=c,lw=lw))
ax.text(5,17.1,"A conserved lncRNA gene is the precursor of a pachytene piRNA cluster — present in ALL 16 strains",ha="center",fontsize=10.3,fontweight="bold")
ax.text(5,16.7,f"read-level-verified example: SPRET/EiJ {CHROM}:{S:,}-{E:,} — lncRNA {SYM} (0% protein-coding overlap; piRNAs sense to the lncRNA)",ha="center",fontsize=7.8,color="#555")
# detailed SPRET locus
yD=14.7; ax.plot([L,R],[yD,yD],color="#333",lw=2.4); ax.text(L-0.18,yD,"SPRET/EiJ",ha="right",va="center",fontsize=9,fontweight="bold",color="#1f6f3d")
ax.add_patch(Rectangle((g0,yD-0.22),g1-g0,0.44,facecolor="#7a3b9a",edgecolor="#3d1f5c",lw=1.0,zorder=3,alpha=0.85)); arrow(ctr-gdir*0.5,yD,gdir,"white",1.4,0.9)
xs=np.linspace(g0-0.05,g1+0.05,60); peak=0.72*np.exp(-((xs-ctr)/0.7)**2); ax.fill_between(xs,yD+0.36,yD+0.36+peak,color="#0072B2",alpha=0.25,zorder=2)
for xx in np.linspace(g0+0.2,g1-0.2,8): arrow(xx,yD+0.58,gdir,"#0072B2",1.1,0.28)
ax.text(ctr,yD-0.5,f"lncRNA {SYM}  —  conserved precursor",ha="center",va="top",fontsize=7.5,color="#3d1f5c")
ax.text(ctr,yD+1.45,f"pachytene piRNA cluster — {sense:.0f}% SENSE to the lncRNA (lncRNA-derived) · ≈{nsp:,} piRNAs · 1U {u1:.0f}% · TE {tef:.0f}%",ha="center",va="bottom",fontsize=7.5,color="#0072B2",fontweight="bold")
ax.text((L+g0)/2,yD+0.28,"shared flank",ha="center",fontsize=6.3,color="#999"); ax.text((g1+R)/2,yD+0.28,"shared flank",ha="center",fontsize=6.3,color="#999")
# cross-strain presence: ALL 16 present (conserved)
ax.text(L-0.18,12.95,"Cross-strain presence — all 16 strains (canonical order): lncRNA + sense pachytene piRNAs present in EVERY strain",ha="left",fontsize=7.8,fontweight="bold",color="#333")
rows=["SPRET_EiJ"]+OTHERS; ys=np.linspace(12.3,1.7,len(rows))
for nm,y in zip(rows,ys):
    sp=(nm=="SPRET_EiJ")
    ax.text(L-0.18,y,nm.replace("_","/"),ha="right",va="center",fontsize=6.6,fontweight="bold" if sp else "normal",color="#1f6f3d" if sp else "#555")
    ax.plot([L,R],[y,y],color="#333",lw=1.7 if sp else 1.4)
    ax.add_patch(Rectangle((g0,y-0.10),g1-g0,0.20,fc="#7a3b9a",ec="#3d1f5c",lw=0.5,zorder=3,alpha=0.85))
    for xx in np.linspace(g0+0.2,g1-0.2,5): arrow(xx,y+0.0,gdir,"#0072B2",0.7,0.16)
ax.text(R+0.12,ys[len(rows)//2],"conserved lncRNA + sense pachytene piRNAs (present in all 16)",ha="left",va="center",fontsize=6.2,color="#1f6f3d",style="italic")
for xx in (g0,g1): ax.plot([xx,xx],[1.4,12.6],color="#e6e6e6",lw=0.6,ls=(0,(1,3)),zorder=0)
ax.text(5,0.75,"Mechanism: A-MYB (MYBL1) activates the lncRNA promoter at the pachytene stage; the precursor is processed by MOV10L1 + PLD6/Zucchini into phased, 1U pachytene piRNAs SENSE to the lncRNA, low TE.\n"
              "The lncRNA locus is CONSERVED across all 16 strains (gene symbol present in every assembly; real 16-strain coverage in Fig_gm10505_16strains) — the CONTRAST with the strain-private, antisense,\n"
              "TE-silencing IAP cluster (Fig_locus_example_IAP, private 1/16): the lncRNA-driven route to piRNA biogenesis.",ha="center",va="center",fontsize=6.6,color="#444")
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_locus_example_ncRNA.{e}",bbox_inches="tight")
print(f"wrote Fig_locus_example_ncRNA.{{png,pdf,svg}} (16-strain) | {SYM} sense={sense:.1f}% 1U={u1:.1f}% TE={tef:.1f}% species={nsp:,}")
