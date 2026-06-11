#!/usr/bin/env python3
"""Concept figure: the FOUR routes by which a piRNA can look 'strain-specific', and which are genuinely
unique. Each panel = the same locus across 3 strains (SPRET/CAST/C57); the piRNA shown is a REAL
representative SPRET sequence for that route (base resolution, from the Step-4 classification), and the
cross-strain pattern (expressed-in-all / SNP-allele / silent-elsewhere / private-locus) shows why only
routes 3 & 4 are genuinely unique. Counts = real SPRET_EiJ Step-4 classification."""
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
d=pd.read_csv(f"{U}/step4/SPRET_EiJ.step4_classified.csv.gz")
def rep_seq(klass):
    sub=d[d.klass==klass]; s=sub[sub.sequence.str[0]=="T"].sequence   # prefer 1U
    s=(s if len(s) else sub.sequence).tolist()
    comp=lambda q: max(q.count(b) for b in "ACGT")/len(q)             # least homopolymer-dominated = representative
    return min(s,key=comp)
SEQ={1:rep_seq("expressed elsewhere (exact)"),2:rep_seq("SNP-variant of expressed (1-3mm)"),
     3:rep_seq("unique: conserved-but-silent"),4:rep_seq("unique: strain-private locus")}
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,axs=plt.subplots(2,2,figsize=(12,8.2),dpi=300)
STR=["SPRET/EiJ","CAST/EiJ","C57BL/6NJ"]; YS=[0.80,0.68,0.56]
def mini(ax,seq,xc,y,w=0.26,faded=False,ticks=0):       # real piRNA as a base-coloured mini-strip on a strain line
    bw=w/len(seq); x0=xc-w/2
    for k,b in enumerate(seq):
        ax.add_patch(Rectangle((x0+k*bw,y),bw*0.92,0.034,fc=NT.get(b,"#999"),ec="none",alpha=0.32 if faded else 1.0))
    for t in range(ticks): ax.plot([x0+(len(seq)//3+t*4)*bw]*2,[y+0.038,y+0.06],color="black",lw=0.8)
def seqrow(ax,seq,y,label):                              # the real sequence at readable base resolution
    x0=0.06; w=0.88; bw=w/len(seq)
    for k,b in enumerate(seq):
        ax.add_patch(Rectangle((x0+k*bw,y),bw*0.9,0.072,fc=NT.get(b,"#999"),ec="white",lw=0.2))
        ax.text(x0+k*bw+bw*0.45,y+0.036,b,ha="center",va="center",fontsize=4.3,color="white",fontweight="bold")
    ax.text(x0,y+0.095,label,fontsize=5.6,color="#333",va="bottom")
def gline(ax,y,L,R,gap=None):
    if gap: ax.plot([L,gap[0]],[y,y],color="#333",lw=1.6); ax.plot([gap[1],R],[y,y],color="#333",lw=1.6); ax.plot(gap,[y,y],color="#ccc",lw=0.8,ls=(0,(2,2)))
    else: ax.plot([L,R],[y,y],color="#333",lw=1.6)
def panel(ax,title,verdict,vcol,note,count,seq,draw):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.01,0.01),0.98,0.98,boxstyle="round,pad=0.01",fc="none",ec="#ddd",lw=1))
    ax.text(0.04,0.96,title,fontsize=9.0,fontweight="bold",va="top")
    ax.text(0.965,0.965,verdict,fontsize=8.0,fontweight="bold",color="white",ha="right",va="top",bbox=dict(boxstyle="round,pad=0.25",fc=vcol,ec="none"))
    L,R=0.20,0.80; ins0,ins1=0.44,0.56
    for s,y in zip(STR,YS):
        ax.text(0.17,y,s,ha="right",va="center",fontsize=6.8,color="#333")
        draw(ax,s,y,L,R,ins0,ins1)
    ax.text(0.04,0.42,note,fontsize=6.6,color="#444",va="top")
    ax.text(0.965,0.42,f"SPRET: {count:,}",fontsize=7.2,color=vcol,ha="right",va="top",fontweight="bold")
    seqrow(ax,seq,0.12,"real SPRET piRNA (this route):")
def d1(ax,s,y,L,R,i0,i1): gline(ax,y,L,R); mini(ax,SEQ[1],(i0+i1)/2,y+0.010)               # exact in all strains
panel(axs[0,0],"1 · Exact sequence expressed in ALL strains","NOT unique","#9e9e9e",
      "Identical piRNA made by every strain (0-mm match to an expressed\nsequence elsewhere) — looked strain-specific only from a threshold.",24124,SEQ[1],d1)
def d2(ax,s,y,L,R,i0,i1):
    gline(ax,y,L,R)
    if s=="SPRET/EiJ": mini(ax,SEQ[2],(i0+i1)/2,y+0.010,ticks=2)
    else: mini(ax,SEQ[2],(i0+i1)/2,y+0.010,faded=True)
panel(axs[0,1],"2 · SNP-variant of a conserved, expressed piRNA","NOT unique","#E69F00",
      "Same locus in all strains; SPRET's copy carries 1–3 SNPs (ticks) — a\nstrain ALLELE of a conserved piRNA, not a new one. (≤3 mm = Poisson cutoff.)",100513,SEQ[2],d2)
def d3(ax,s,y,L,R,i0,i1):
    gline(ax,y,L,R)
    if s=="SPRET/EiJ": mini(ax,SEQ[3],(i0+i1)/2,y+0.010)
    else: ax.text((i0+i1)/2,y+0.028,"silent",ha="center",fontsize=5.6,color="#bbb",style="italic")
panel(axs[1,0],"3 · Conserved locus, expressed only in one strain","UNIQUE (expression)","#0072B2",
      "Locus EXISTS in all strains but is transcribed into piRNAs only in\nSPRET — regulatory / expression divergence. The LARGEST unique class.",140969,SEQ[3],d3)
def d4(ax,s,y,L,R,i0,i1):
    if s=="SPRET/EiJ": gline(ax,y,L,R); ax.add_patch(Rectangle((i0,y-0.017),i1-i0,0.034,fc="#C0392B",ec="#7a1f15",lw=0.5)); mini(ax,SEQ[4],(i0+i1)/2,y+0.020)
    else: gline(ax,y,L,R,gap=(i0,i1))
panel(axs[1,1],"4 · Strain-private locus (e.g. TE insertion)","UNIQUE (locus gain)","#C0392B",
      "A TE inserted only in SPRET (absent in others) becomes a new piRNA\nsource — the L1 (chr2) & IAP (chr7) examples. TE-driven; the minority.",61208,SEQ[4],d4)
fig.suptitle("Four routes to a 'strain-specific' piRNA — only routes 3 & 4 are genuinely unique (expression is the criterion); each piRNA shown is a real SPRET Step-4 sequence",
             fontsize=10.2,fontweight="bold",y=1.004)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_concept_four_routes.{e}",bbox_inches="tight")
print("wrote Fig_concept_four_routes.{png,pdf,svg} | route seqs:",{k:SEQ[k] for k in SEQ})
