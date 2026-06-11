#!/usr/bin/env python3
"""Comprehensive locus figure: EXAMPLE + ANATOMY + NUCLEOTIDE RESOLUTION (real read depth).
(1) cross-strain status strip (strain-private); (2) TE content; (3) strand-resolved piRNA coverage
(antisense-to-TE=silencing red, sense grey) + 1U inset; (4) a REAL piRNA pair at base resolution
(ping-pong: sense+antisense 10-nt 5' overlap with 1U/10A; or the top piRNA if no pair).
Usage: make_locus_full.py <tag> <chrom> <start> <end> <te_label> <te_strand>"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import numpy as np, pandas as pd, pysam
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
BAM="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam"
tag,CHROM,S,E,TELAB,TEST=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],sys.argv[6]
N=E-S; comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))
nb=200; fwd=np.zeros(nb); rev=np.zeros(nb); first=[]
sense=defaultdict(Counter); anti=defaultdict(Counter)
bam=pysam.AlignmentFile(BAM,"rb")
for a in bam.fetch(CHROM,S,E):
    if a.is_unmapped: continue
    L=a.reference_end-a.reference_start
    if not(24<=L<=32) or not a.query_sequence: continue
    b0=int((a.reference_start-S)/N*nb); b1=int((a.reference_end-S)/N*nb)
    for b in range(max(0,b0),min(nb,b1+1)): (rev if a.is_reverse else fwd)[b]+=1
    if a.is_reverse: anti[a.reference_end-1][rc(a.query_sequence)]+=1; first.append(comp.get(a.query_sequence[-1],"N"))
    else: sense[a.reference_start][a.query_sequence]+=1; first.append(a.query_sequence[0])
bam.close()
fc=pd.Series(first).value_counts(normalize=True)*100
antiC,senseC=(fwd,rev) if TEST=="-" else (rev,fwd)
# find a real ping-pong pair (sense@i, antisense@i+9), prefer 1U(anti)+10A(sense)
pair=None
for i in sense:
    if i+9 in anti:
        ss,sc=sense[i].most_common(1)[0]; aq,acn=anti[i+9].most_common(1)[0]
        if len(ss)>=10 and len(aq)>=10:
            sgn=(sc+acn)*(3 if(aq[0]=="T" and ss[9]=="A") else 1)
            if pair is None or sgn>pair[0]: pair=(sgn,i,ss,sc,aq,acn)
TE=pd.read_csv(f"{U}/sense_antisense/SPRET_EiJ.TE_stranded.bed",sep="\t",header=None,names=["c","s","e","fam","sc","st"])
TE=TE[(TE.c==CHROM)&(TE.s<E)&(TE.e>S)]
TECOL={"LINE/L1":"#E69F00","LTR/ERVK":"#6a3d9a","LTR/ERVL-MaLR":"#b15928","SINE/Alu":"#33a02c","SINE/B2":"#1f78b4","LTR/ERVL":"#a6cee3"}
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
plt.rcParams.update({"font.family":"DejaVu Sans"})
fig=plt.figure(figsize=(10.5,8.4),dpi=300)
gs=fig.add_gridspec(4,4,height_ratios=[0.9,0.4,1.9,1.5],hspace=0.5,wspace=0.3)
axS=fig.add_subplot(gs[0,:]); axT=fig.add_subplot(gs[1,:3]); axC=fig.add_subplot(gs[2,:3]); ax1=fig.add_subplot(gs[2,3]); axN=fig.add_subplot(gs[3,:])
# strip
axS.set_xlim(0,16); axS.set_ylim(0,1.7); axS.axis("off"); order=[s for s in STRAIN_ORDER if s!="C57BL_6"]
for i,s in enumerate(order):
    has=s=="SPRET_EiJ"; axS.add_patch(Rectangle((i+0.1,0.55),0.8,0.6,fc=("#C0392B" if has else "#eee"),ec="#999",lw=0.5))
    if has: axS.text(i+0.5,0.85,"●",ha="center",va="center",fontsize=8,color="white")
    axS.text(i+0.5,0.4,s.replace("_","/"),ha="right",va="top",fontsize=5,rotation=45)
axS.text(0,1.5,"Insertion + piRNA cluster present in SPRET/EiJ only (1/16) — strain-private",fontsize=8.5,fontweight="bold")
# TE
axT.set_xlim(S,E); axT.set_ylim(0,1); axT.axis("off"); axT.text(S,0.78,"TEs:",fontsize=6.5)
for _,r in TE.iterrows(): axT.add_patch(Rectangle((max(r.s,S),0.3),min(r.e,E)-max(r.s,S),0.4,fc=TECOL.get(r.fam.split("|")[-1],"#ccc"),ec="none"))
hx=0.12
for cl,col in list(TECOL.items())[:5]: axT.add_patch(Rectangle((S+N*hx,0.78),N*0.012,0.18,fc=col,ec="none")); axT.text(S+N*(hx+0.014),0.87,cl,fontsize=5,va="center"); hx+=0.16
# coverage
x=np.linspace(S,E,nb)
axC.fill_between(x,0,antiC,step="mid",color="#C0392B",alpha=0.85,label="antisense to TE (silencing)")
axC.fill_between(x,0,-senseC,step="mid",color="#999",alpha=0.85,label="sense to TE")
axC.axhline(0,color="#333",lw=0.8); axC.set_xlim(S,E); axC.set_ylabel("piRNA coverage",fontsize=8); axC.set_xlabel(f"{CHROM.split('#')[-1]} position",fontsize=8)
axC.ticklabel_format(axis="x",style="plain"); axC.tick_params(labelsize=6); axC.legend(fontsize=6,frameon=False,loc="upper right"); axC.spines[['top','right']].set_visible(False)
# 1U
vals=[fc.get(n,0) for n in ["A","C","G","T"]]; ax1.bar(["A","C","G","U"],vals,color=["#bbb","#bbb","#bbb","#C0392B"]); ax1.set_ylim(0,max(vals+[1])*1.25)
for i,v in enumerate(vals): ax1.text(i,v+1,f"{v:.0f}",ha="center",fontsize=6)
ax1.set_title("5′ nt (1U)",fontsize=7.4,fontweight="bold"); ax1.set_ylabel("%",fontsize=7); ax1.tick_params(labelsize=7); ax1.spines[['top','right']].set_visible(False)
# nucleotide panel
axN.axis("off"); bw=0.42
if pair:
    _,i,ss,sc,aq,acn=pair
    for k,b in enumerate(ss):
        axN.add_patch(Rectangle((k*bw,1.5),bw*0.9,0.5,fc=NT.get(b,"#999"),ec="white",lw=0.4)); axN.text(k*bw+bw*0.45,1.75,b,ha="center",va="center",fontsize=7,color="white",fontweight="bold")
    axN.text(-0.35,1.75,"5′",ha="right",fontsize=8); axN.text(len(ss)*bw+0.1,1.75,"3′ sense",ha="left",fontsize=7.5,color="#555")
    x0=9*bw
    for k,b in enumerate(aq):
        xx=x0-k*bw; axN.add_patch(Rectangle((xx+bw*0.1,0.7),bw*0.9,0.5,fc=NT.get(b,"#999"),ec="white",lw=0.4)); axN.text(xx+bw*0.55,0.95,b,ha="center",va="center",fontsize=7,color="white",fontweight="bold")
    axN.text(x0+bw+0.1,0.95,"5′",ha="left",fontsize=8); axN.text(x0-len(aq)*bw,0.95,"3′ antisense",ha="right",fontsize=7.5,color="#555")
    axN.add_patch(Rectangle((0,0.65),10*bw,1.45,fc="#fff3cd",ec="#e0a800",lw=1,zorder=0))
    for k in range(10): axN.plot([k*bw+bw*0.45]*2,[1.45,1.25],color="#888",lw=0.6)
    axN.annotate("1U",xy=(x0+bw*0.55,1.2),xytext=(x0+2,0.3),fontsize=8,color="#e31a1c",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#e31a1c",lw=1))
    axN.annotate("10A",xy=(9*bw+bw*0.45,2.0),xytext=(9*bw+1.6,2.5),fontsize=8,color="#1f78b4",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#1f78b4",lw=1))
    axN.text(5*bw,2.35,"nucleotide resolution: real ping-pong pair — 10-nt complementary 5′ overlap, 1U + 10A",ha="center",fontsize=8,fontweight="bold",color="#b8860b")
    # genomic coordinate ruler aligned to the nucleotides (genomic = i + x/bw; both strands share the frame)
    gmin=i+10-len(aq); gmax=i+len(ss)-1; yr=0.45
    axN.plot([(gmin-i)*bw,(gmax-i)*bw],[yr,yr],color="#555",lw=0.8)
    g=int(np.ceil(gmin/5.0)*5)
    while g<=gmax:
        xx=(g-i)*bw; axN.plot([xx,xx],[yr,yr-0.05],color="#555",lw=0.8)
        axN.text(xx,yr-0.08,f"{g:,}",ha="right",va="top",fontsize=5,rotation=90,color="#555"); g+=5
    axN.text(((gmin+gmax)/2-i)*bw,yr-0.66,f"{CHROM.split('#')[-1]} position",ha="center",va="top",fontsize=7,color="#333")
    axN.set_xlim((10-len(aq))*bw-2.0,len(ss)*bw+3); axN.set_ylim(-0.62,2.7)
if pair:   # zoom-in callout: tie the nucleotide panel to its true position on the coverage track (different scales)
    pc=(gmin+gmax)/2.0
    axC.axvspan(gmin,gmax,color="#e0a800",alpha=0.25,zorder=1)
    axC.text(pc,axC.get_ylim()[1]*0.98,f"ping-pong pair\n{CHROM.split('#')[-1]}:{gmin:,}\n(zoom ↓)",ha="center",va="top",fontsize=5.4,color="#b8860b",fontweight="bold",linespacing=0.95)
    yC=axC.get_ylim()[0]; xL=(10-len(aq))*bw-0.3; xR=(len(ss)-1)*bw+0.3; yNt=2.12
    for xa,xb in [(gmin,xL),(gmax,xR)]:
        fig.add_artist(ConnectionPatch(xyA=(xa,yC),coordsA=axC.transData,xyB=(xb,yNt),coordsB=axN.transData,color="#e0a800",lw=0.9,ls=(0,(4,2)),alpha=0.9,zorder=20))
fig.suptitle(f"Strain-private {TELAB} → piRNA cluster: example + anatomy + nucleotide resolution (SPRET/EiJ {CHROM.split('#')[-1]}:{S:,}-{E:,})",fontsize=9.6,fontweight="bold",y=0.995)
for ext in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_locus_full_{tag}.{ext}",bbox_inches="tight")
print(f"[{tag}] 1U={fc.get('T',0):.0f}% anti={antiC.sum():.0f} sense={senseC.sum():.0f} TEs={len(TE)} pair={'yes' if pair else 'no'}")
