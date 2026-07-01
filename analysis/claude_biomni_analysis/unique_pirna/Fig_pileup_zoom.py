#!/usr/bin/env python3
"""COORDINATE-ANCHORED nucleotide pileup at the SPRET L1 cluster (real reads, IGV-style). Reads are
drawn at their ACTUAL genomic positions over the reference sequence, so the biogenesis signatures are
read directly off the coordinates: ping-pong (a sense + antisense read whose 5' ends are 10 nt apart),
phasing (consecutive same-strand reads ~27 nt apart, head-to-tail), and 1U at the 5' ends. Not a floating
schematic — every base is at its genomic coordinate."""
import numpy as np, pysam
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
BAM="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam"
FA="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/SPRET_EiJ_chromosomes_MT.fasta"
CHROM,RS,RE="SPRET_EiJ#1#chr2",150450743,150521192
# collect read 5' (genomic start for +, genomic start for - too; we track strand) with counts, and full read span
plus=defaultdict(int); minus=defaultdict(int); span={}
bam=pysam.AlignmentFile(BAM,"rb")
for a in bam.fetch(CHROM,RS,RE):
    if a.is_unmapped: continue
    L=a.reference_end-a.reference_start
    if not(25<=L<=32): continue
    key=(a.reference_start,a.reference_end,a.is_reverse)
    span[key]=span.get(key,0)+1
    if a.is_reverse: minus[a.reference_end-1]+=1     # 5' (genomic) of - read = right end
    else: plus[a.reference_start]+=1                  # 5' of + read = left end
bam.close()
# find a window: ping-pong (sense 5'=i, antisense 5'=i+9) AND a phased sense read near i+27
best=None
for i in list(plus):
    if plus[i]<50: continue
    if (i+9) in minus and minus[i+9]>=5:
        ph=max((plus.get(i+d,0) for d in range(24,31)),default=0)   # phased sense ~27 downstream
        sc=plus[i]+minus[i+9]+ph*3
        if best is None or sc>best[0]: best=(sc,i)
i=best[1]; W0,W1=i-10,i+45
ref=pysam.FastaFile(FA).fetch(CHROM,W0,W1).upper()
# pick reads to display: pin-pong pair + phased sense + a couple abundant, all overlapping window
cand=[]
for (s,e,rev),c in span.items():
    if e<=W0 or s>=W1: continue
    cand.append((c,s,e,rev))
cand.sort(reverse=True)
# ensure the ping-pong pair present: sense starting at i, antisense ending at i+9
disp=[]; seen=set()
def add(s,e,rev):
    k=(s,e,rev)
    if k in span and k not in seen: disp.append((s,e,rev,span[k])); seen.add(k)
# sense@i (longest), antisense@i+9 (longest), phased sense@~i+27
for (s,e,rev),c in sorted(span.items(),key=lambda kv:-kv[1]):
    if not rev and s==i and len(disp)<1: add(s,e,rev)
for (s,e,rev),c in sorted(span.items(),key=lambda kv:-kv[1]):
    if rev and e-1==i+9: add(s,e,rev); break
for (s,e,rev),c in sorted(span.items(),key=lambda kv:-kv[1]):
    if not rev and abs(s-(i+27))<=2: add(s,e,rev); break
for c,s,e,rev in cand:
    if len(disp)>=7: break
    add(s,e,rev)
print(f"window {CHROM}:{W0}-{W1} ({W1-W0}bp), ping-pong sense@{i} antisense5'@{i+9}, {len(disp)} reads")
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
plt.rcParams.update({"font.family":"DejaVu Sans Mono"})
fig,ax=plt.subplots(figsize=(10,5.0),dpi=170); ax.set_xlim(W0-1,W1+1)
nrow=len(disp); ax.set_ylim(-1.5,nrow+1.5); ax.axis("off")
def cell(x,y,b,h=0.8,fc=None,txt=True,ec="white"):
    ax.add_patch(Rectangle((x,y),1,h,fc=fc or NT.get(b,"#999"),ec=ec,lw=0.4,zorder=3))
    if txt: ax.text(x+0.5,y+h/2,b,ha="center",va="center",fontsize=6.5,color="white",fontweight="bold")
# reference row (y=0)
for k,b in enumerate(ref): cell(W0+k,0,b)
ax.text(W0-0.6,0.4,"ref",ha="right",va="center",fontsize=8,fontweight="bold")
# reads, sorted by start, stacked
disp.sort(key=lambda r:r[0])
for ri,(s,e,rev,c) in enumerate(disp):
    y=ri+1
    seq=pysam.FastaFile(FA).fetch(CHROM,s,e).upper()    # genomic-forward bases (aligned to ref)
    for k in range(s,e): cell(k,y,seq[k-s])
    # 5' marker: + read left, - read right
    if rev:
        ax.annotate("",xy=(e+0.9,y+0.4),xytext=(e+0.1,y+0.4),arrowprops=dict(arrowstyle="-|>",color="#C0392B",lw=1.6))
        ax.text(e+1.1,y+0.4,f"− 5′ (×{c})",ha="left",va="center",fontsize=6.5,color="#C0392B")
    else:
        ax.annotate("",xy=(s-0.9,y+0.4),xytext=(s-0.1,y+0.4),arrowprops=dict(arrowstyle="-|>",color="#0072B2",lw=1.6))
        ax.text(s-1.1,y+0.4,f"+ 5′ (×{c})",ha="right",va="center",fontsize=6.5,color="#0072B2")
# annotate ping-pong (sense 5'=i, antisense 5'=i+9: 10-nt overlap i..i+9)
ax.add_patch(Rectangle((i,-0.1),10,nrow+1.2,fc="#fff3cd",ec="#e0a800",lw=1,zorder=0))
ax.text(i+5,nrow+1.0,"10-nt 5′ overlap → ping-pong",ha="center",fontsize=8,fontweight="bold",color="#b8860b")
ax.annotate("",xy=(i+9.5,nrow+0.4),xytext=(i+0.5,nrow+0.4),arrowprops=dict(arrowstyle="<->",color="#b8860b",lw=1))
# phasing arrow (i -> i+27)
ax.annotate("",xy=(i+27,-0.9),xytext=(i,-0.9),arrowprops=dict(arrowstyle="<->",color="#0072B2",lw=1))
ax.text(i+13,-1.25,"~27 nt → phasing (head-to-tail)",ha="center",fontsize=7.5,color="#0072B2")
ax.set_title(f"Coordinate-anchored piRNA pileup at the SPRET L1 cluster — ping-pong, phasing & 1U read off the genome ({CHROM.split('#')[-1]}:{W0:,}-{W1:,})",fontsize=9.8,fontweight="bold")
# coordinate ticks
for xt in range(W0,W1+1,10): ax.text(xt+0.5,-1.7,f"{xt:,}",ha="center",fontsize=6,color="#666")
for ext in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_pileup_zoom.{ext}",bbox_inches="tight")
print("wrote Fig_pileup_zoom.{png,pdf,svg}")
