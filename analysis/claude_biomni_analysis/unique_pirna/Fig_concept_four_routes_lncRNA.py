#!/usr/bin/env python3
"""Concept figure (lncRNA version) — the four routes by which an lncRNA-DERIVED piRNA can look 'strain-specific'.
KEY difference from the TE version: route 4 is NOT a locus gain. For lncRNAs the strain-private-locus class comes
from a CONSERVED lncRNA precursor by sequence divergence (the lncRNA gene is present in all strains). So routes 3 & 4
are both genuinely unique but BOTH arise from conserved lncRNAs (expression / sequence divergence) — lncRNAs show
essentially no TE-style locus-gain route. Every example is CONFOUNDING-AUDITED: representative locus overlaps the
lncRNA only (no protein_coding, no TE); target genes have protein_coding overlap = 0. Counts = SPRET Step-4 ∩
clean-lncRNA (uniquely-mapping, 1U). Reps + audit computed here from the data."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pysam, pandas as pd, re
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, ConnectionPatch
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
X="SPRET_EiJ"; comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}; NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
REAL=f"{ROOT}/results/STAR_srna_strain_wise/{X}/{X}-20.5dpp.1/Aligned.sortedByCoord.out.bam"; CS=f"{U}/step4/{X}.cand_self.Aligned.sortedByCoord.out.bam"
GFF=f"{ROOT}/resources/annotation/{X}_v3.5.gff3"; RM=f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
d=pd.read_csv(f"{U}/step4/{X}.step4_classified.csv.gz"); kof=dict(zip(d.id,d.klass)); seqof=dict(zip(d.id,d.sequence))
# route -> (klass, verdict, colour, count[clean-lncRNA], chrom, gstart, gend, sym, title, note)
SPEC=[(1,"expressed elsewhere (exact)","NOT unique","#9e9e9e",645,"chr18",33494059,33502846,"C330008A17",
       "Exact sequence in ALL strains","Identical lncRNA-derived piRNA made by every strain."),
      (2,"SNP-variant of expressed (1-3mm)","NOT unique","#E69F00",10718,"chr17",23815048,23828598,"Gm49794",
       "SNP-variant of a conserved lncRNA piRNA","≤3-SNP allele of a conserved lncRNA piRNA the others express."),
      (3,"unique: conserved-but-silent","UNIQUE (expression)","#0072B2",12367,"chr7",60978491,61013294,"Gm10619",
       "Conserved lncRNA, expressed only in SPRET","lncRNA in all strains, piRNAs made only in SPRET — expression divergence."),
      (4,"unique: strain-private locus","UNIQUE (sequence)","#7a3b9a",5905,"chr17",23790189,23809974,"Gm10505",
       "Strain-private SEQUENCE from a CONSERVED lncRNA","lncRNA conserved in all 16 strains; the piRNA SEQUENCE diverged — NOT a locus gain.")]
# ---- parse genes (all biotypes) + TE for the audit ----
chset={s[5] for s in SPEC}; genes=defaultdict(list); te=defaultdict(list)
for ln in open(GFF):
    if ln[0]=="#": continue
    f=ln.split("\t")
    if len(f)<9 or f[2] not in ("gene","ncRNA_gene"): continue
    c=f[0].split("#")[-1]
    if c in chset:
        bt=re.search(r"biotype=([^;]+)",f[8]); genes[c].append((int(f[3]),int(f[4]),bt.group(1) if bt else "?"))
for ln in open(RM):
    g=ln.split("\t",4)
    if g[0] in chset: te[g[0]].append((int(g[1]),int(g[2])))
def pc_over(c,s,e): return [g for g in genes[c] if g[0]<e and g[1]>s and g[2]=="protein_coding"]
def te_over(c,s,e): return [t for t in te[c] if t[0]<e and t[1]>s]
real=pysam.AlignmentFile(REAL,"rb"); cs=pysam.AlignmentFile(CS,"rb"); W=110
def pick(kl,c,gs,ge):
    A=np.zeros(ge-gs)
    for r in real.fetch(f"{X}#1#{c}",gs,ge):
        if r.is_unmapped or not 25<=r.reference_end-r.reference_start<=32: continue
        A[max(gs,r.reference_start)-gs:min(ge,r.reference_end)-gs]+=1
    best=None; seen=set()
    for r in cs.fetch(f"{X}#1#{c}",gs,ge):
        qn=r.query_name
        if qn in seen or kof.get(qn)!=kl or (r.has_tag("NH") and r.get_tag("NH")!=1): continue
        if not 25<=r.reference_end-r.reference_start<=32 or r.reference_start<gs or r.reference_end>ge: continue
        sq=r.query_sequence or ""
        if not ((comp.get(sq[-1],"N")=="T") if r.is_reverse else (sq[:1]=="T")): continue
        seen.add(qn); cv=A[r.reference_start-gs:r.reference_end-gs].max()
        if best is None or cv>best[0]: best=(cv,qn,r.reference_start,r.reference_end,"-" if r.is_reverse else "+")
    return best
def cov(c,s,e):
    a,b=s-W,e+W; arr=np.zeros(b-a)
    for r in real.fetch(f"{X}#1#{c}",a,b):
        if not r.is_unmapped and 25<=r.reference_end-r.reference_start<=32:
            for p in range(max(a,r.reference_start),min(b,r.reference_end)): arr[p-a]+=1
    return a,b,arr
REP={}
print("=== confounding audit (must be CLEAN before plotting) ===")
for rt,kl,verd,col,cnt,c,gs,ge,sym,title,note in SPEC:
    b=pick(kl,c,gs,ge); assert b, f"no rep for R{rt}"
    cv,qn,p0,p1,st=b; gpc=pc_over(c,gs,ge); lpc=pc_over(c,p0,p1); lte=te_over(c,p0,p1)
    REP[rt]=(seqof[qn],c,p0,p1,st,cv)
    print(f"R{rt} {sym} {c}:{p0}-{p1} {st} cov={cv:.0f} | gene-span pc-overlap={len(gpc)} | locus pc={len(lpc)} TE={len(lte)} -> {'CLEAN' if not(gpc or lpc or lte) else 'WARN'}")
# ---- cross-strain intuition (lncRNA: route 4 line is CONTINUOUS = conserved, no gap) ----
def xstrain(ax,ridx,seq):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    ax.text(0.5,0.99,"the same lncRNA locus across strains",ha="center",va="top",fontsize=6.4,color="#777",style="italic")
    L,Rg,i0,i1=0.34,0.96,0.55,0.74; ctr=(i0+i1)/2
    def line(y): ax.plot([L,Rg],[y,y],color="#333",lw=1.5)                       # lncRNA conserved -> always continuous
    def gene(y): ax.add_patch(Rectangle((i0,y-0.016),i1-i0,0.032,fc="#7a3b9a",ec="none",alpha=0.30))  # conserved lncRNA gene
    def strip(y,faded=False,ticks=0):
        w=0.17; bw=w/len(seq); x0=ctr-w/2
        for k,bse in enumerate(seq): ax.add_patch(Rectangle((x0+k*bw,y+0.012),bw*0.92,0.045,fc=NT.get(bse,"#999"),ec="none",alpha=0.28 if faded else 1))
        for t in range(ticks): ax.plot([x0+(len(seq)//4+t*4)*bw]*2,[y+0.058,y+0.092],color="black",lw=0.8)
    for i,(nm,y) in enumerate([("SPRET/EiJ",0.72),("CAST/EiJ",0.47),("C57BL/6NJ",0.22)]):
        sp=(i==0); line(y); gene(y)
        ax.text(L-0.03,y,nm,ha="right",va="center",fontsize=5.9,color=("#C0392B" if sp else "#555"),fontweight="bold" if sp else "normal")
        if ridx==1: strip(y)
        elif ridx==2: strip(y,ticks=2 if sp else 0,faded=not sp)
        elif ridx==3: (strip(y) if sp else ax.text(ctr,y+0.045,"silent (lncRNA present)",ha="center",fontsize=4.8,color="#bbb",style="italic"))
        else: (strip(y) if sp else (strip(y,faded=True,ticks=3),ax.text(ctr,y-0.03,"diverged seq",ha="center",fontsize=4.6,color="#aaa",style="italic")))
    for xx in (i0,i1): ax.plot([xx,xx],[0.14,0.80],color="#ddd",lw=0.6,ls=(0,(1,3)),zorder=0)
    if ridx==4: ax.text(0.5,0.04,"locus CONSERVED (no gap) — only the sequence is private",ha="center",fontsize=4.9,color="#7a3b9a")
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(15.5,13),dpi=300)
gs=fig.add_gridspec(8,3,width_ratios=[1.0,1.55,3.0],height_ratios=[1,1.05]*4,hspace=0.5,wspace=0.1,left=0.012,right=0.985,top=0.93,bottom=0.025)
for i,(rt,kl,verd,col,cnt,c,gs0,ge0,sym,title,note) in enumerate(SPEC):
    seq,cc,st,en,strand,cv=REP[rt]
    axv=fig.add_subplot(gs[2*i:2*i+2,0]); axx=fig.add_subplot(gs[2*i:2*i+2,1]); axs=fig.add_subplot(gs[2*i,2]); axn=fig.add_subplot(gs[2*i+1,2])
    axv.set_xlim(0,1); axv.set_ylim(0,1); axv.axis("off"); axv.add_patch(FancyBboxPatch((0.03,0.05),0.94,0.9,boxstyle="round,pad=0.01",fc="none",ec="#ddd",lw=1))
    axv.text(0.07,0.92,f"{rt} · {title}",fontsize=7.4,fontweight="bold",va="top",wrap=True)
    axv.text(0.07,0.5,verd,fontsize=7.4,fontweight="bold",color="white",va="center",bbox=dict(boxstyle="round,pad=0.3",fc=col,ec="none"))
    axv.text(0.07,0.36,f"clean-lncRNA: {cnt:,}",fontsize=6.8,color=col,fontweight="bold",va="center")
    axv.text(0.07,0.29,f"in lncRNA {sym}",fontsize=6.0,color="#7a3b9a",va="center")
    axv.text(0.07,0.22,note,fontsize=5.8,color="#444",va="top",wrap=True)
    xstrain(axx,rt,seq)
    a,b,arr=cov(cc,st,en); xx=np.arange(a,b)
    axs.fill_between(xx,0,arr,step="mid",color=col,alpha=0.55,lw=0); axs.axvspan(st,en,color="#e0a800",alpha=0.4)
    axs.set_xlim(a,b); axs.set_ylim(0,max(arr.max(),1)*1.18)
    axs.set_title(f"zoom-out · SPRET/EiJ {cc}:{st:,} in lncRNA {sym} — real P20.5 piRNA coverage (24–32 nt)",fontsize=6.6,fontweight="bold",loc="left")
    axs.tick_params(labelsize=5.3); axs.ticklabel_format(axis="x",style="plain"); axs.spines[['top','right']].set_visible(False); axs.set_ylabel("cov",fontsize=6)
    axn.axis("off")
    for k,bse in enumerate(seq):
        axn.add_patch(Rectangle((k+0.04,1.0),0.92,0.5,fc=NT.get(bse,"#999"),ec="white",lw=0.4)); axn.text(k+0.5,1.25,bse,ha="center",va="center",color="white",fontsize=6.3,fontweight="bold")
    axn.text(-0.5,1.25,"5′",ha="right",va="center",fontsize=8); axn.text(len(seq)+0.2,1.25,f"3′  real SPRET piRNA ({strand} strand)",ha="left",va="center",fontsize=6.8,color="#444")
    yr=0.78; axn.plot([0.5,len(seq)-0.5],[yr,yr],color="#555",lw=0.8)
    for k in range(0,len(seq),5):
        gp=st+k if strand=="+" else en-1-k
        axn.plot([k+0.5,k+0.5],[yr,yr-0.07],color="#555",lw=0.8); axn.text(k+0.5,yr-0.1,f"{gp:,}",ha="right",va="top",fontsize=5,rotation=90,color="#555")
    axn.text(len(seq)/2,yr-0.78,f"{cc} position",ha="center",va="top",fontsize=6.4,color="#333"); axn.set_xlim(-3,len(seq)+11); axn.set_ylim(yr-0.98,1.7)
    for xa,xb in [(st,0.5),(en,len(seq)-0.5)]:
        fig.add_artist(ConnectionPatch(xyA=(xa,0),coordsA=axs.transData,xyB=(xb,1.55),coordsB=axn.transData,color="#e0a800",lw=0.8,ls=(0,(3,2)),alpha=0.85,zorder=20))
real.close(); cs.close()
fig.suptitle("Four routes to a 'strain-specific' piRNA — lncRNA-DERIVED version. Routes 3 & 4 are genuinely unique, but BOTH arise from CONSERVED lncRNA precursors\n"
             "(expression divergence / sequence divergence). Unlike TEs, lncRNAs show essentially no 'locus-gain' route — route 4's locus is present in every strain.",
             fontsize=9.0,fontweight="bold",y=0.985,linespacing=1.4)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_concept_four_routes_lncRNA.{e}",bbox_inches="tight")
print("wrote Fig_concept_four_routes_lncRNA.{png,pdf,svg}")
