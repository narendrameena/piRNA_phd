#!/usr/bin/env python3
"""Concept figure — the FOUR routes by which a piRNA can look 'strain-specific', and which are genuinely
unique. Three columns per route: (L) verdict + real Step-4 count; (M) the INTUITION — the same locus across
3 strains (SPRET/CAST/C57) showing the route's defining cross-strain pattern (expressed-in-all / SNP-allele
/ silent-elsewhere / private); (R) the DATA — the real representative SPRET piRNA at its true locus, zoomed
from real P20.5 coverage to single-base + genomic-coordinate resolution. Counts = real SPRET_EiJ Step-4
classification; representatives = uniquely-mapping 1U candidates (step4_classified ∩ cand_self BAM), loci verified 2026-06-11."""
import numpy as np, pysam
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, ConnectionPatch
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"
BAM="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam"
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
R=[(1,"Exact sequence in ALL strains","NOT unique","#9e9e9e",24124,"TAAATGCATCTGAAGCCTTGGACGGTCC","chr12",2096248,2096276,"+",
    "Identical piRNA made by every strain (0-mm match elsewhere)."),
   (2,"SNP-variant of a conserved piRNA","NOT unique","#E69F00",100513,"TAAAGGTCACTCTGAATCCTGCGAGGCT","chr19",2790047,2790075,"+",
    "SPRET ≤3-SNP allele of a conserved piRNA the others express."),
   (3,"Conserved locus, expressed only in SPRET","UNIQUE (expression)","#0072B2",140969,"TAACGGTATCAGGTAGGTAGCACCTCTC","chr10",75578163,75578191,"-",
    "Locus in all strains, transcribed only in SPRET — expression divergence."),
   (4,"Strain-private locus (e.g. TE insertion)","UNIQUE (locus gain)","#C0392B",61208,"TAAACAGTGTCAGGGCAGTCTGTACCTC","chr14",56538402,56538430,"-",
    "Locus only in SPRET (absent in others) — TE-driven birth; the minority.")]
W=110; bam=pysam.AlignmentFile(BAM,"rb")
def cov(chrom,s,e):
    a,b=s-W,e+W; c=np.zeros(b-a)
    for r in bam.fetch(f"SPRET_EiJ#1#{chrom}",a,b):
        if r.is_unmapped: continue
        L=r.reference_end-r.reference_start
        if 24<=L<=32:
            for p in range(max(a,r.reference_start),min(b,r.reference_end)): c[p-a]+=1
    return a,b,c
def xstrain(ax,ridx,seq):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    ax.text(0.5,0.99,"the same locus across strains",ha="center",va="top",fontsize=6.4,color="#777",style="italic")
    L,Rg,i0,i1=0.34,0.96,0.55,0.74; ctr=(i0+i1)/2
    def line(y,gap=False):
        if gap: ax.plot([L,i0],[y,y],color="#333",lw=1.5); ax.plot([i1,Rg],[y,y],color="#333",lw=1.5); ax.plot([i0,i1],[y,y],color="#ccc",lw=0.8,ls=(0,(2,2)))
        else: ax.plot([L,Rg],[y,y],color="#333",lw=1.5)
    def strip(y,faded=False,ticks=0):
        w=0.17; bw=w/len(seq); x0=ctr-w/2
        for k,bse in enumerate(seq): ax.add_patch(Rectangle((x0+k*bw,y+0.012),bw*0.92,0.045,fc=NT.get(bse,"#999"),ec="none",alpha=0.3 if faded else 1))
        for t in range(ticks): ax.plot([x0+(len(seq)//3+t*5)*bw]*2,[y+0.058,y+0.092],color="black",lw=0.8)
    for i,(nm,y) in enumerate([("SPRET/EiJ",0.72),("CAST/EiJ",0.47),("C57BL/6NJ",0.22)]):
        sp=(i==0); ax.text(L-0.03,y,nm,ha="right",va="center",fontsize=5.9,color=("#C0392B" if sp else "#555"),fontweight="bold" if sp else "normal")
        if ridx==1: line(y); strip(y)
        elif ridx==2: line(y); strip(y,ticks=2 if sp else 0,faded=not sp)
        elif ridx==3:
            line(y); (strip(y) if sp else ax.text(ctr,y+0.045,"silent",ha="center",fontsize=5.2,color="#bbb",style="italic"))
        else:
            if sp: line(y); ax.add_patch(Rectangle((i0,y-0.018),i1-i0,0.036,fc="#C0392B",ec="#7a1f15",lw=0.4,alpha=0.45)); strip(y)
            else: line(y,gap=True)
    for xx in (i0,i1): ax.plot([xx,xx],[0.14,0.80],color="#ddd",lw=0.6,ls=(0,(1,3)),zorder=0)
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(15.5,13),dpi=300)
gs=fig.add_gridspec(8,3,width_ratios=[1.0,1.55,3.0],height_ratios=[1,1.05]*4,hspace=0.5,wspace=0.1,left=0.012,right=0.985,top=0.94,bottom=0.025)
for i,(n,title,verdict,col,count,seq,chrom,st,en,strand,note) in enumerate(R):
    axv=fig.add_subplot(gs[2*i:2*i+2,0]); axx=fig.add_subplot(gs[2*i:2*i+2,1]); axs=fig.add_subplot(gs[2*i,2]); axn=fig.add_subplot(gs[2*i+1,2])
    axv.set_xlim(0,1); axv.set_ylim(0,1); axv.axis("off")
    axv.add_patch(FancyBboxPatch((0.03,0.05),0.94,0.9,boxstyle="round,pad=0.01",fc="none",ec="#ddd",lw=1))
    axv.text(0.07,0.9,f"{n} · {title}",fontsize=7.6,fontweight="bold",va="top",wrap=True)
    axv.text(0.07,0.52,verdict,fontsize=7.6,fontweight="bold",color="white",va="center",bbox=dict(boxstyle="round,pad=0.3",fc=col,ec="none"))
    axv.text(0.07,0.37,f"SPRET: {count:,}",fontsize=7.0,color=col,fontweight="bold",va="center")
    axv.text(0.07,0.28,note,fontsize=6.0,color="#444",va="top",wrap=True)
    xstrain(axx,n,seq)
    a,b,c=cov(chrom,st,en); x=np.arange(a,b)
    axs.fill_between(x,0,c,step="mid",color=col,alpha=0.55,lw=0); axs.axvspan(st,en,color="#e0a800",alpha=0.4)
    axs.set_xlim(a,b); axs.set_ylim(0,max(c.max(),1)*1.18)
    axs.set_title(f"zoom-out · SPRET/EiJ {chrom}:{st:,} — real P20.5 piRNA coverage (24–32 nt)",fontsize=6.8,fontweight="bold",loc="left")
    axs.tick_params(labelsize=5.3); axs.ticklabel_format(axis="x",style="plain"); axs.spines[['top','right']].set_visible(False); axs.set_ylabel("cov",fontsize=6)
    axn.axis("off")
    for k,bse in enumerate(seq):
        axn.add_patch(Rectangle((k+0.04,1.0),0.92,0.5,fc=NT.get(bse,"#999"),ec="white",lw=0.4)); axn.text(k+0.5,1.25,bse,ha="center",va="center",color="white",fontsize=6.3,fontweight="bold")
    axn.text(-0.5,1.25,"5′",ha="right",va="center",fontsize=8); axn.text(len(seq)+0.2,1.25,f"3′  real SPRET piRNA ({strand} strand)",ha="left",va="center",fontsize=6.8,color="#444")
    yr=0.78; axn.plot([0.5,len(seq)-0.5],[yr,yr],color="#555",lw=0.8)
    for k in range(0,len(seq),5):
        g=st+k if strand=="+" else en-1-k
        axn.plot([k+0.5,k+0.5],[yr,yr-0.07],color="#555",lw=0.8); axn.text(k+0.5,yr-0.1,f"{g:,}",ha="right",va="top",fontsize=5,rotation=90,color="#555")
    axn.text(len(seq)/2,yr-0.78,f"{chrom} position",ha="center",va="top",fontsize=6.4,color="#333"); axn.set_xlim(-3,len(seq)+11); axn.set_ylim(yr-0.98,1.7)
    for xa,xb in [(st,0.5),(en,len(seq)-0.5)]:
        fig.add_artist(ConnectionPatch(xyA=(xa,0),coordsA=axs.transData,xyB=(xb,1.55),coordsB=axn.transData,color="#e0a800",lw=0.8,ls=(0,(3,2)),alpha=0.85,zorder=20))
bam.close()
fig.suptitle("Four routes to a 'strain-specific' piRNA — only routes 3 & 4 are genuinely unique (expression is the criterion).\nMiddle = the cross-strain intuition (same locus in 3 strains); right = the real SPRET piRNA at its true locus, zoomed to single-base + coordinate resolution.",
             fontsize=9.2,fontweight="bold",y=0.99,linespacing=1.4)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_concept_four_routes.{e}",bbox_inches="tight")
print("wrote Fig_concept_four_routes.{png,pdf,svg} (3-col: verdict | cross-strain intuition | coverage+nucleotide)")
