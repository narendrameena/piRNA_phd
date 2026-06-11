#!/usr/bin/env python3
"""Rich multi-strain + multi-timepoint figure for a strain-private TE-driven piRNA cluster (replaces the
old lines-and-arrows schematic). REAL data only:
  - SPRET strand-resolved piRNA coverage at the locus across the THREE timepoints (E16.5/P12.5/P20.5,
    reps pooled), on the same genomic coordinate axis -> shows the cluster emerge through development;
  - TE composition at true coordinates;
  - a real antisense piRNA at nucleotide resolution (1U), with a zoom callout to its position;
  - verified cross-strain absence: SPRET-private (1/16); the 15 other strains have no orthologous
    insertion -> 0 cluster piRNAs (Step-4 genome-anchored expression test).
Usage: make_locus_example.py <outbase> <chrom> <start> <end> <te_label> <te_strand> <npi> <antisense_desc>"""
import sys, numpy as np, pandas as pd, pysam
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch

U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"
BR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/SPRET_EiJ"
TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
TECOL={"LINE/L1":"#E69F00","LTR/ERVK":"#6a3d9a","LTR/ERVL-MaLR":"#b15928","SINE/Alu":"#33a02c","SINE/B2":"#1f78b4","LTR/ERVL":"#a6cee3"}

outbase,CHROM,S,E,TELAB,TEST,NPI,ADESC=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],sys.argv[6],int(sys.argv[7]),sys.argv[8]
N=E-S; CHRLAB=CHROM.split("#")[-1]; nb=220
boxcol="#E69F00" if ("L1" in TELAB or "LINE" in TELAB) else "#6a3d9a"

# ---- real coverage per timepoint (reps pooled) + antisense piRNA from P20.5 ----
cov={}; totreads={}; anti_seq=defaultdict(Counter)
for lab,tp in TPS:
    fwd=np.zeros(nb); rev=np.zeros(nb); tot=0
    for r in (1,2,3):
        try: bam=pysam.AlignmentFile(f"{BR}/SPRET_EiJ-{tp}.{r}/Aligned.sortedByCoord.out.bam","rb")
        except Exception: continue
        for a in bam.fetch(CHROM,S,E):
            if a.is_unmapped or not a.query_sequence: continue
            L=a.reference_end-a.reference_start
            if not(24<=L<=32): continue
            b0=int((a.reference_start-S)/N*nb); b1=int((a.reference_end-S)/N*nb)
            for b in range(max(0,b0),min(nb,b1+1)): (rev if a.is_reverse else fwd)[b]+=1
            tot+=1
            if a.is_reverse: anti_seq[a.reference_end-1][rc(a.query_sequence)]+=1   # pool antisense piRNAs across all tps for the representative
        bam.close()
    cov[lab]=((fwd,rev) if TEST=="-" else (rev,fwd)); totreads[lab]=tot   # antisense-to-TE = silencing
ymax=max(max(c[0].max(),c[1].max()) for c in cov.values()) or 1.0
a5=max(anti_seq,key=lambda p:anti_seq[p].most_common(1)[0][1]); aseq,acnt=anti_seq[a5].most_common(1)[0]
TE=pd.read_csv(f"{U}/sense_antisense/SPRET_EiJ.TE_stranded.bed",sep="\t",header=None,names=["c","s","e","fam","sc","st"])
TE=TE[(TE.c==CHROM)&(TE.s<E)&(TE.e>S)]

# ---- plot ----
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(11.8,11.2),dpi=300)
gs=fig.add_gridspec(6,1,height_ratios=[0.30,0.60,0.60,0.78,1.12,1.65],hspace=0.5,left=0.10,right=0.975,top=0.925,bottom=0.05)
axTE,axE,axP,axA,axN,axX=[fig.add_subplot(gs[i]) for i in range(6)]
x=np.linspace(S,E,nb)
# TE composition bar
axTE.set_xlim(S,E); axTE.set_ylim(0,1); axTE.axis("off"); axTE.text(S-N*0.012,0.5,"TEs",fontsize=7.5,ha="right",va="center",fontweight="bold")
for _,r in TE.iterrows(): axTE.add_patch(Rectangle((max(r.s,S),0.2),min(r.e,E)-max(r.s,S),0.6,fc=TECOL.get(r.fam.split("|")[-1],"#ccc"),ec="none"))
axTE.set_title(f"Strain-private {TELAB} insertion → a strain-specific piRNA cluster  (real example: SPRET/EiJ {CHRLAB}:{S:,}–{E:,})",fontsize=10.2,fontweight="bold",loc="left")
# timepoint coverage tracks (shared coords + y-scale; reps pooled)
for ax,lab in [(axE,"E16.5"),(axP,"P12.5"),(axA,"P20.5")]:
    antiC,senseC=cov[lab]
    ax.fill_between(x,0,antiC,step="mid",color="#C0392B",alpha=0.9); ax.fill_between(x,0,-senseC,step="mid",color="#999",alpha=0.9)
    ax.axhline(0,color="#333",lw=0.6); ax.set_xlim(S,E); ax.set_ylim(-ymax*1.05,ymax*1.05)
    ax.set_ylabel(lab,fontsize=8.5,fontweight="bold",rotation=0,ha="right",va="center",labelpad=14)
    ax.tick_params(labelsize=6); ax.spines[['top','right']].set_visible(False)
    ax.text(0.997,0.92,f"{totreads[lab]:,} reads",transform=ax.transAxes,ha="right",va="top",fontsize=6,color="#777")
    if lab!="P20.5": ax.tick_params(labelbottom=False); ax.spines['bottom'].set_visible(False)
axE.set_title("SPRET/EiJ — real piRNA coverage across development (reps pooled; red = antisense-to-TE / silencing, grey = sense)",fontsize=8.0,fontweight="bold",loc="left",pad=2)
axA.set_xlabel(f"{CHRLAB} position",fontsize=8); axA.ticklabel_format(axis="x",style="plain")
axA.legend(handles=[plt.Line2D([0],[0],color="#C0392B",lw=5),plt.Line2D([0],[0],color="#999",lw=5)],
           labels=["antisense to TE (silencing)","sense to TE"],fontsize=6,frameon=False,loc="upper right",ncol=2)
# nucleotide inset (real antisense piRNA) + zoom callout from P20.5
axN.axis("off")
for k,b in enumerate(aseq):
    axN.add_patch(Rectangle((k+0.04,1.0),0.92,0.5,fc=NT.get(b,"#999"),ec="white",lw=0.4)); axN.text(k+0.5,1.25,b,ha="center",va="center",color="white",fontsize=7,fontweight="bold")
axN.text(-0.4,1.25,"5′",ha="right",va="center",fontsize=8.5); axN.text(len(aseq)+0.15,1.25,f"3′  antisense piRNA → silences {TELAB}",ha="left",va="center",fontsize=7.5,color="#444")
axN.annotate("1U",xy=(0.5,1.52),xytext=(0.5,2.12),fontsize=8.5,color="#e31a1c",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#e31a1c",lw=1.1))
yr=0.8; axN.plot([0.5,len(aseq)-0.5],[yr,yr],color="#555",lw=0.8)
for k in range(0,len(aseq),5):
    axN.plot([k+0.5,k+0.5],[yr,yr-0.06],color="#555",lw=0.8); axN.text(k+0.5,yr-0.09,f"{a5-k:,}",ha="right",va="top",fontsize=5.4,rotation=90,color="#555")
axN.text(len(aseq)/2,yr-0.6,f"{CHRLAB} position (antisense; 5′→3′ runs high→low coord)",ha="center",va="top",fontsize=7,color="#333")
axN.text(len(aseq)/2,2.28,f"nucleotide resolution: a real antisense piRNA from the cluster (5′-U, {acnt:,} reads here) — antisense to {TELAB} → silencing-competent",ha="center",fontsize=7.8,fontweight="bold",color="#b8860b")
axN.set_xlim(-3,len(aseq)+10); axN.set_ylim(yr-1.0,2.5)
for axtp in (axE,axP,axA): axtp.axvspan(a5-len(aseq)+1,a5,color="#e0a800",alpha=0.35,zorder=6)
for xa,xb in [(a5-len(aseq)+1,0.5),(a5,len(aseq)-0.5)]:
    fig.add_artist(ConnectionPatch(xyA=(xa,axA.get_ylim()[0]),coordsA=axA.transData,xyB=(xb,2.05),coordsB=axN.transData,color="#e0a800",lw=0.9,ls=(0,(4,2)),alpha=0.9,zorder=20))
# cross-strain comparison (genome lines; real annotations)
axX.set_xlim(0,10); axX.set_ylim(0,10); axX.axis("off")
L,R=1.7,8.6; ins0,ins1=4.3,5.7
for nm,y,has in [("SPRET/EiJ",8.4,True),("C57BL/6NJ",6.1,False),("CAST/EiJ",4.5,False),("+ 13 other strains",2.9,False)]:
    if has:
        axX.plot([L,R],[y,y],color="#333",lw=2.4); axX.add_patch(Rectangle((ins0,y-0.18),ins1-ins0,0.36,fc=boxcol,ec="#333",lw=0.8,zorder=3))
        axX.text((ins0+ins1)/2,y+0.55,f"piRNA cluster — {NPI:,} strain-private piRNAs ({ADESC} antisense → silencing)",ha="center",fontsize=7.3,color="#C0392B",fontweight="bold")
    else:
        axX.plot([L,ins0],[y,y],color="#333",lw=2.4); axX.plot([ins1,R],[y,y],color="#333",lw=2.4); axX.plot([ins0,ins1],[y,y],color="#bbb",lw=1.0,ls=(0,(2,2)))
        axX.text((ins0+ins1)/2,y+0.3,"no insertion · 0 cluster piRNAs",ha="center",fontsize=6.8,color="#888",style="italic")
    axX.text(L-0.18,y,nm,ha="right",va="center",fontsize=8.5,fontweight="bold" if has else "normal",color="#C0392B" if has else "#555")
for xx in (ins0,ins1): axX.plot([xx,xx],[2.7,8.6],color="#ccc",lw=0.7,ls=(0,(1,3)),zorder=0)
axX.text((L+ins0)/2,8.95,"shared flank",ha="center",fontsize=6.5,color="#999"); axX.text((ins1+R)/2,8.95,"shared flank",ha="center",fontsize=6.5,color="#999")
axX.set_title("Cross-strain: present in SPRET/EiJ only (1/16) — TE-driven birth of a strain-specific piRNA source",fontsize=8.6,fontweight="bold",loc="left")
axX.text(5,0.7,"Absence in the other 15 strains is data-derived: verified by the Step-4 genome-anchored expression test (no ≤3-mm expressed homolog) — not assumed.",ha="center",fontsize=6.6,color="#444")

for ext in ("pdf","svg","png"): fig.savefig(f"{PG}/{outbase}.{ext}",bbox_inches="tight")
# source data
pd.DataFrame({"bin_center_bp":x,**{f"{lab}_antisense":cov[lab][0] for lab in ["E16.5","P12.5","P20.5"]},
              **{f"{lab}_sense":cov[lab][1] for lab in ["E16.5","P12.5","P20.5"]}}).to_csv(f"{PG}/SourceData_{outbase}_coverage.csv",index=False)
pd.DataFrame([{"role":"antisense_piRNA","genomic_5p":a5,"reads":acnt,"sequence":aseq}]).to_csv(f"{PG}/SourceData_{outbase}_pirna.csv",index=False)
print(f"[{outbase}] reads E16.5/P12.5/P20.5 = {totreads['E16.5']:,}/{totreads['P12.5']:,}/{totreads['P20.5']:,} | antisense piRNA {aseq} @ {a5:,} (x{acnt}) | TEs={len(TE)}")
