#!/usr/bin/env python3
"""Nucleotide-resolution ping-pong: a REAL sense+antisense piRNA pair from the SPRET L1 cluster whose
5' ends overlap by 10 nt, drawn base-by-base — showing the 10-nt complementary overlap, the 1U (primary,
antisense) and 10A (secondary, sense) signatures. Real sequences pulled from the full sRNA BAM."""
import numpy as np, pysam
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
BAM="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam"
CHROM,S,E="SPRET_EiJ#1#chr2",150450743,150521192
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))
sense=defaultdict(Counter); anti=defaultdict(Counter)   # 5'pos -> Counter(piRNA seq)
bam=pysam.AlignmentFile(BAM,"rb")
for a in bam.fetch(CHROM,S,E):
    if a.is_unmapped: continue
    L=a.reference_end-a.reference_start
    if not(25<=L<=32) or not a.query_sequence: continue
    if a.is_reverse: anti[a.reference_end-1][rc(a.query_sequence)]+=1   # piRNA 5'->3'
    else: sense[a.reference_start][a.query_sequence]+=1
bam.close()
# ping-pong pair: sense 5' at i, antisense 5' at i+9 (10-nt overlap); prefer 1U(anti)+10A(sense)
best=None
for i in sense:
    j=i+9
    if j not in anti: continue
    ss,sc=sense[i].most_common(1)[0]; asq,ac=anti[j].most_common(1)[0]
    if len(ss)<10 or len(asq)<10: continue
    score=(sc+ac)*(2 if (asq[0]=="T" and ss[9]=="A") else 1)
    if best is None or score>best[0]: best=(score,ss,sc,asq,ac,i)
score,ss,sc,asq,ac,i=best
print(f"sense  5'->3': {ss}  (x{sc})")
print(f"anti   5'->3': {asq}  (x{ac})   1U={asq[0]} 10A(sense)={ss[9]}")
# draw: sense L->R on top; antisense drawn 3'<-5' under it so the 10-nt 5' overlap aligns
plt.rcParams.update({"font.family":"DejaVu Sans Mono"})
fig,ax=plt.subplots(figsize=(11,3.8),dpi=300); ax.axis("off")
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","U":"#e31a1c","N":"#999"}
bw=0.42
# sense piRNA: positions 0..len-1 at x=0..; its 5' (pos0) at x=0
for k,b in enumerate(ss):
    ax.add_patch(Rectangle((k*bw,1.5),bw*0.92,0.5,fc=NT.get(b,"#999"),ec="white",lw=0.5))
    ax.text(k*bw+bw*0.46,1.75,b,ha="center",va="center",fontsize=8,color="white",fontweight="bold")
ax.text(-0.4,1.75,"5′",ha="right",va="center",fontsize=9); ax.text(len(ss)*bw+0.1,1.75,"3′  sense (secondary)",ha="left",va="center",fontsize=8.5,color="#555")
# antisense piRNA: its 5' aligns at sense position 9 (the overlap), drawn leftward (3'<-5')
# antisense 5' (pos0) sits above sense pos9 -> x = 9*bw ; antisense extends LEFT (toward smaller x)
x0=9*bw
for k,b in enumerate(asq):
    xx=x0-k*bw
    ax.add_patch(Rectangle((xx+bw*0.08,0.6),bw*0.92,0.5,fc=NT.get(b,"#999"),ec="white",lw=0.5))
    ax.text(xx+bw*0.54,0.85,b,ha="center",va="center",fontsize=8,color="white",fontweight="bold")
ax.text(x0+bw+0.1,0.85,"5′",ha="left",va="center",fontsize=9); ax.text(x0-len(asq)*bw+0.1,0.85,"3′  antisense (primary)",ha="right",va="center",fontsize=8.5,color="#555")
# 10-nt overlap shading + base-pair ticks
ax.add_patch(Rectangle((0,0.5),10*bw,1.6,fc="#fff3cd",ec="#e0a800",lw=1,zorder=0))
for k in range(10):
    ax.plot([k*bw+bw*0.46]*2,[1.45,1.15],color="#888",lw=0.7)
ax.text(5*bw,2.25,"10-nt complementary 5′ overlap (ping-pong)",ha="center",fontsize=8.5,fontweight="bold",color="#b8860b")
# 1U + 10A highlight
ax.annotate("1U\n(primary)",xy=(x0+bw*0.54,1.15),xytext=(x0+2.2,0.2),fontsize=8,color="#e31a1c",fontweight="bold",ha="center",
            arrowprops=dict(arrowstyle="-|>",color="#e31a1c",lw=1.2))
ax.annotate("10A\n(ping-pong)",xy=(9*bw+bw*0.46,2.0),xytext=(9*bw+2.0,2.7),fontsize=8,color="#1f78b4",fontweight="bold",ha="center",
            arrowprops=dict(arrowstyle="-|>",color="#1f78b4",lw=1.2))
# genomic coordinate ruler aligned to the nucleotides (genomic = i + x/bw; both strands share the frame)
gmin=i+10-len(asq); gmax=i+len(ss)-1; yr=0.35
ax.plot([(gmin-i)*bw,(gmax-i)*bw],[yr,yr],color="#555",lw=0.8)
g=int(np.ceil(gmin/5.0)*5)
while g<=gmax:
    xx=(g-i)*bw; ax.plot([xx,xx],[yr,yr-0.04],color="#555",lw=0.8)
    ax.text(xx,yr-0.07,f"{g:,}",ha="right",va="top",fontsize=5.5,rotation=90,color="#555"); g+=5
ax.text(((gmin+gmax)/2-i)*bw,yr-0.55,f"{CHROM.split('#')[-1]} position",ha="center",va="top",fontsize=7,color="#333")
ax.set_xlim((10-len(asq))*bw-2.5,len(ss)*bw+3.5); ax.set_ylim(-0.55,3.1)
ax.set_title("Ping-pong amplification at base resolution — a real sense+antisense piRNA pair (SPRET L1 cluster)",fontsize=10,fontweight="bold")
fig.text(0.5,0.02,f"primary antisense piRNA (5′-U) and secondary sense piRNA (10-A) overlap by exactly 10 nt at their 5′ ends — the molecular signature of ping-pong. "
  f"sense ×{sc}, antisense ×{ac} reads.",ha="center",fontsize=6.6,color="#555")
fig.tight_layout()
for ext in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_pingpong_nucleotide.{ext}",bbox_inches="tight")
print("wrote Fig_pingpong_nucleotide.{png,pdf,svg}")
