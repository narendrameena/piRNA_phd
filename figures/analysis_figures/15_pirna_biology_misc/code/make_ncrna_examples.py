#!/usr/bin/env python3
"""ncRNA-DRIVEN (pachytene) piRNA cluster examples — real SPRET P20.5 clusters that coincide with an
annotated lncRNA gene. Each example uses a DEDICATED gene-model track BELOW the coverage (no text overlap):
 top  = strand-resolved piRNA coverage (+ up / - down)
 below= source lncRNA gene model (exons + intron line) on its own axis
 right= 5' nt (1U)
Includes 4 unidirectional clusters AND 1 divergent BIDIRECTIONAL cluster (RIKEN 4932441J04: + arm left, - arm
right). Data-driven: clusters from PICB scan, lncRNA + exons from the v3.5 GFF3, coverage from the real BAM,
TE% from RepeatMasker. Architecture (uni/bidirectional) computed from the +/- read balance per cluster."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pysam, re
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"
BAM="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam"
GFF="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/SPRET_EiJ_v3.5.gff3"
RM="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/SPRET_EiJ_repeatmasker.bed"
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
# low-TE lncRNA-overlapping SPRET P20.5 clusters; last is a DIVERGENT BIDIRECTIONAL cluster
EX=[("chr2",91078926,91151480,"Gm13817"),
    ("chr17",23789396,23833400,"Gm10505"),
    ("chr5",144808901,144855100,"Gm21221"),
    ("chr7",60307311,60365235,"RIKEN 4933435G04"),
    ("chr5",54121761,54264726,"RIKEN 4932441J04")]
chroms={c for c,_,_,_ in EX}
tx=defaultdict(list); exons=defaultdict(list)
for ln in open(GFF):
    if ln[0]=="#": continue
    f=ln.rstrip("\n").split("\t")
    if len(f)<9: continue
    c=f[0].split("#")[-1]
    if c not in chroms: continue
    s,e=int(f[3]),int(f[4])
    if f[2]=="lnc_RNA":
        tid=(re.search(r"ID=transcript:([^;]+)",f[8]) or re.search(r"ID=([^;]+)",f[8])).group(1)
        tx[c].append((s,e,tid,f[6]))
    elif f[2]=="exon":
        p=re.search(r"Parent=transcript:([^;]+)",f[8])
        if p: exons[p.group(1)].append((s,e))
rm=defaultdict(list)
for ln in open(RM):
    g=ln.split("\t",3)
    if g[0] in chroms: rm[g[0]].append((int(g[1]),int(g[2])))
bam=pysam.AlignmentFile(BAM,"rb")
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(12,14),dpi=300)
outer=fig.add_gridspec(5,1,hspace=0.62)
src_rows=[]
for i,(c,S,E,lab) in enumerate(EX):
    inner=outer[i].subgridspec(2,4,height_ratios=[1,0.22],hspace=0.05,wspace=0.34)
    axC=fig.add_subplot(inner[0,:3]); axG=fig.add_subplot(inner[1,:3],sharex=axC); ax1=fig.add_subplot(inner[:,3])
    N=E-S; nb=240; fwd=np.zeros(nb); rev=np.zeros(nb); first=[]; chrom=f"SPRET_EiJ#1#{c}"
    for a in bam.fetch(chrom,S,E):
        if a.is_unmapped or not a.query_sequence: continue
        L=a.reference_end-a.reference_start
        if not(25<=L<=32): continue
        b0=int((a.reference_start-S)/N*nb); b1=int((a.reference_end-S)/N*nb)
        for b in range(max(0,b0),min(nb,b1+1)): (rev if a.is_reverse else fwd)[b]+=1
        first.append(comp.get(a.query_sequence[-1],"N") if a.is_reverse else a.query_sequence[0])
    fc=pd.Series(first).value_counts(normalize=True)*100 if first else pd.Series(dtype=float)
    teb=sum(min(b,E)-max(a,S) for a,b in rm[c] if a<E and b>S); tef=teb/N
    bal=min(fwd.sum(),rev.sum())/max(fwd.sum(),rev.sum(),1)
    arch="bidirectional (divergent)" if bal>=0.3 else "unidirectional"
    x=np.linspace(S,E,nb)
    axC.fill_between(x,0,fwd,step="mid",color="#0072B2",alpha=0.85,label="+ strand")
    axC.fill_between(x,0,-rev,step="mid",color="#56B4E9",alpha=0.9,label="− strand")
    axC.axhline(0,color="#333",lw=0.7); axC.set_xlim(S,E); axC.tick_params(labelsize=6,labelbottom=False)
    axC.set_ylabel("piRNA cov\n(+↑ / −↓)",fontsize=7); axC.spines[['top','right']].set_visible(False)
    axC.set_title(f"SPRET/EiJ {c}:{S:,}-{E:,} ({N//1000} kb) — source lncRNA: {lab}  |  TE {tef*100:.0f}%  |  {arch} (bal {bal:.2f})",fontsize=7.6,fontweight="bold",loc="left")
    axC.legend(fontsize=6,frameon=False,loc="upper right",ncol=2)
    # dedicated gene-model track (own axis, no overlap)
    lncs=sorted([t for t in tx[c] if t[0]<E and t[1]>S],key=lambda t:-(min(t[1],E)-max(t[0],S)))
    axG.set_xlim(S,E); axG.set_ylim(0,1); axG.set_yticks([]); axG.spines[['top','right','left']].set_visible(False)
    axG.tick_params(labelsize=6); axG.ticklabel_format(axis="x",style="plain"); axG.set_xlabel(f"{c} position",fontsize=7)
    if lncs:
        ts,te,tid,strand=lncs[0]
        axG.plot([max(ts,S),min(te,E)],[0.5,0.5],color="#333",lw=0.9,zorder=1)
        for (es,ee) in exons.get(tid,[]):
            if es<E and ee>S: axG.add_patch(Rectangle((max(es,S),0.18),min(ee,E)-max(es,S),0.64,fc="#7a3b9a",ec="none",zorder=2))
        axG.annotate("",xy=(min(te,E) if strand=="+" else max(ts,S),0.5),xytext=((min(te,E)+max(ts,S))/2,0.5),
                     arrowprops=dict(arrowstyle="-|>",color="#333",lw=0.8),zorder=3)
        axG.text(S,1.18,f"lncRNA precursor: {lab} ({strand})",fontsize=6,color="#7a3b9a",va="bottom")
    vals=[fc.get(n,0) for n in ["A","C","G","T"]]; ax1.bar(["A","C","G","U"],vals,color=["#bbb","#bbb","#bbb","#0072B2"]); ax1.set_ylim(0,max(vals+[1])*1.28)
    for j,v in enumerate(vals): ax1.text(j,v+1,f"{v:.0f}",ha="center",fontsize=5.5)
    ax1.set_title("5′ nt (1U)",fontsize=6.8,fontweight="bold"); ax1.set_ylabel("%",fontsize=6); ax1.tick_params(labelsize=6); ax1.spines[['top','right']].set_visible(False)
    src_rows.append(dict(locus=f"{c}:{S}-{E}",lncRNA=lab,kb=N//1000,TE_pct=round(tef*100,1),fwd_cov=int(fwd.sum()*N/nb),rev_cov=int(rev.sum()*N/nb),bal=round(bal,3),architecture=arch,U1_pct=round(fc.get("T",0),1)))
bam.close()
pd.DataFrame(src_rows).to_csv(f"{PG}/SourceData_ncrna_examples.csv",index=False)
fig.suptitle("ncRNA-DRIVEN (pachytene) piRNA clusters — real SPRET P20.5 clusters over annotated lncRNA genes (strong 1U; 4 unidirectional + 1 divergent bidirectional)",fontsize=9.6,fontweight="bold",y=0.998)
fig.text(0.5,0.004,"Source lncRNA = the annotated lncRNA the cluster overlaps (purple gene model, own track). Pachytene piRNAs are processed from lncRNA precursors (A-MYB-driven). "
  "Bottom panel = a DIVERGENT bidirectional cluster (+ arm left / − arm right). Overlap is annotation-based evidence the lncRNA is the precursor; architecture computed from +/- balance.",ha="center",fontsize=6,color="#555")
for ext in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_ncrna_examples.{ext}",bbox_inches="tight")
print("wrote Fig_ncrna_examples.{png,pdf,svg} + SourceData_ncrna_examples.csv");
print(pd.DataFrame(src_rows).to_string())
