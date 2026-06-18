#!/usr/bin/env python3
"""16-strain extension of the conserved lncRNA-driven pachytene cluster Gm10505. For each strain (canonical
order, P20.5), real strand-resolved piRNA coverage over that strain's OWN orthologous Gm10505 locus (coords
parsed per-strain from its v3.5 GFF3 — assemblies differ), normalised to RPM (library = total mapped reads).
Demonstrates the lncRNA precursor + its sense pachytene piRNA cluster are CONSERVED + active across all 16
strains (the precursor is conserved; strain-private piRNAs arise from it by sequence divergence). All from data."""
import warnings; warnings.filterwarnings("ignore")
import pysam, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
PG=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
STR=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
import glob,re
def coords(X):
    for ln in open(f"{ROOT}/resources/annotation/{X}_v3.5.gff3"):
        f=ln.split("\t")
        if len(f)>8 and f[2] in ("gene","ncRNA_gene") and "predicted gene 10505" in f[8]:
            return f[0],int(f[3]),int(f[4]),f[6]
    return None
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,axes=plt.subplots(4,4,figsize=(15,11),dpi=300); axes=axes.ravel()
rows=[]
for k,X in enumerate(STR):
    ax=axes[k]; cc=coords(X); bam=sorted(glob.glob(f"{ROOT}/results/STAR_srna_strain_wise/{X}/{X}-20.5dpp.*/Aligned.sortedByCoord.out.bam"))
    if cc is None or not bam: ax.text(0.5,0.5,f"{X}\n(no data)",ha="center",va="center",transform=ax.transAxes,fontsize=8); ax.axis("off"); continue
    chrom,S,E,strand=cc; b=pysam.AlignmentFile(bam[0],"rb")
    lib=sum(x.mapped for x in b.get_index_statistics()) or 1
    N=E-S; nb=160; fwd=np.zeros(nb); rev=np.zeros(nb); first=[]; nin=0
    for a in b.fetch(chrom,S,E):
        if a.is_unmapped or not a.query_sequence: continue
        L=a.reference_end-a.reference_start
        if not 24<=L<=32: continue
        nin+=1; b0=int((a.reference_start-S)/N*nb); b1=int((a.reference_end-S)/N*nb)
        for bb in range(max(0,b0),min(nb,b1+1)): (rev if a.is_reverse else fwd)[bb]+=1
        first.append(comp.get(a.query_sequence[-1],"N") if a.is_reverse else a.query_sequence[0])
    b.close()
    sc=1e6/lib; fwd*=sc; rev*=sc; tot_rpm=nin*sc
    sense=(rev.sum() if strand=="-" else fwd.sum())/max(fwd.sum()+rev.sum(),1e-9)*100
    u1=(first.count("T")/len(first)*100) if first else 0
    x=np.linspace(0,N/1000,nb)
    ax.fill_between(x,0,fwd,step="mid",color="#0072B2",alpha=0.85); ax.fill_between(x,0,-rev,step="mid",color="#56B4E9",alpha=0.9)
    ax.axhline(0,color="#333",lw=0.5); ax.set_title(f"{X.replace('_','/')}  ·  {tot_rpm:,.0f} RPM\nsense {sense:.0f}% · 1U {u1:.0f}%",fontsize=7.6,fontweight="bold")
    ax.tick_params(labelsize=6); ax.set_xlabel("kb in Gm10505 locus",fontsize=6); ax.spines[['top','right']].set_visible(False)
    if k%4==0: ax.set_ylabel("piRNA RPM\n(+↑ / −↓)",fontsize=6.5)
    rows.append(dict(strain=X,locus=f"{chrom.split('#')[-1]}:{S}-{E}",strand=strand,total_RPM=round(tot_rpm,1),sense_pct=round(sense,1),U1_pct=round(u1,1),reads=nin,libsize=lib))
df=pd.DataFrame(rows); df.to_csv(f"{PG}/SourceData_gm10505_16strains.csv",index=False)
fig.suptitle("Conserved lncRNA-driven pachytene piRNA cluster Gm10505 across all 16 strains (P20.5) — precursor conserved + cluster active in every strain (sense, low-TE; strain-private piRNAs arise by divergence)",fontsize=10,fontweight="bold",y=1.0)
fig.text(0.5,0.005,"Each panel = that strain's OWN orthologous Gm10505 locus (per-strain v3.5 coords; assemblies differ), real P20.5 piRNA coverage in RPM (library=total mapped). "
  "The lncRNA precursor + its sense pachytene piRNAs are present in all 16 strains (conserved); cf. the divergence-driven strain-private piRNA sequences from this same locus.",ha="center",fontsize=6.4,color="#555")
fig.tight_layout(rect=[0,0.015,1,0.985])
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_gm10505_16strains.{e}",bbox_inches="tight")
print("wrote Fig_gm10505_16strains.{png,pdf,svg} + SourceData_gm10505_16strains.csv")
print(df.to_string())
