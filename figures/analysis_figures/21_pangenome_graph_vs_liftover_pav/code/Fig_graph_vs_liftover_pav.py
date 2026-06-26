#!/usr/bin/env python3
"""THEME 21 Fig 1 — pangenome GRAPH (odgi sequence-PAV) vs HAL-LIFTOVER (cluster-PAV) at 42,384 master piRNA-cluster
loci. The liftover 'private/dispensable' classification is overwhelmingly REFERENCE BIAS / epigenetic: 99% of
liftover-private clusters are sequence-SHARED in the graph, 52-57% of strain-locus events are SILENCING (sequence
present, NO cluster) vs only 3-8% genetic loss. The graph SEPARATES epigenetic silencing from genetic loss; liftover
(one linear reference) conflates them. Robust across PAV-ratio thresholds 0.5-0.95."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
T="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav"; D=f"{T}/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
lift=pd.read_csv(f"{D}/liftover_pav_matrix.tsv",sep="\t",dtype={"chrom":str})
g=pd.read_csv(f"{D}/graph_pav_matrix.tsv",sep="\t").rename(columns={"name":"locus_id"})
m=lift.merge(g[["locus_id"]+S].rename(columns={s:"g_"+s for s in S}),on="locus_id")
lp=m[S].values.astype(int); gpr=m[["g_"+s for s in S]].values
comp=pd.read_csv(f"{D}/graph_vs_liftover_comparison.tsv",sep="\t")
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(13,9.7),dpi=300)
ORD=["core","dispensable","private"]
# A: confusion heatmap (liftover class rows x graph class cols)
conf=pd.crosstab(comp.liftover_class,comp.graph_class).reindex(index=ORD,columns=["core","dispensable","private","absent"],fill_value=0)
im=axA.imshow(conf.values,cmap="Blues",aspect="auto",norm=matplotlib.colors.LogNorm(vmin=1,vmax=conf.values.max()))
axA.set_xticks(range(4)); axA.set_xticklabels(["core","disp.","private","absent"],fontsize=9); axA.set_yticks(range(3)); axA.set_yticklabels(["core","dispensable","private"],fontsize=9)
for i in range(3):
    for j in range(4):
        v=conf.values[i,j]; axA.text(j,i,f"{v:,}",ha="center",va="center",fontsize=8.5,fontweight=("bold" if (i>0 and j==0) else "normal"),color=("#B2182B" if (i>0 and j==0) else ("white" if v>conf.values.max()*0.2 else "#333")))
axA.set_xlabel("GRAPH class (sequence PAV)",fontsize=9.5); axA.set_ylabel("HAL-LIFTOVER class (cluster PAV)",fontsize=9.5)
axA.text(0,2.62,"liftover private/dispensable\n-> graph CORE = reference bias",ha="center",fontsize=7.2,color="#B2182B",style="italic")
axA.set_title("A  Liftover under-calls sharing: most 'dispensable' & 'private'\nclusters are sequence-CORE in the graph",fontsize=9.6,fontweight="bold",loc="left")
# B: strain-locus event breakdown across thresholds (robustness)
ths=[0.5,0.8,0.9,0.95]; conc=[];sil=[];loss=[]
for th in ths:
    gp=(gpr>=th).astype(int); tot=lp.size
    conc.append(100*((lp==1)&(gp==1)).sum()/tot); sil.append(100*((lp==0)&(gp==1)).sum()/tot); loss.append(100*((lp==0)&(gp==0)).sum()/tot)
x=np.arange(len(ths));
axB.bar(x,conc,color="#4393C3",label="concordant (cluster + sequence)",edgecolor="white")
axB.bar(x,sil,bottom=conc,color="#E8852B",label="SILENCING (sequence+, cluster−)",edgecolor="white")
axB.bar(x,loss,bottom=np.array(conc)+np.array(sil),color="#B2182B",label="genetic loss (no sequence)",edgecolor="white")
for i in range(len(ths)): axB.text(i,conc[i]/2,f"{conc[i]:.0f}",ha="center",va="center",fontsize=7.5,color="white",fontweight="bold"); axB.text(i,conc[i]+sil[i]/2,f"{sil[i]:.0f}",ha="center",va="center",fontsize=7.5,color="white",fontweight="bold"); axB.text(i,conc[i]+sil[i]+loss[i]/2+1.5,f"{loss[i]:.0f}",ha="center",va="center",fontsize=7,color="#B2182B",fontweight="bold")
axB.set_xticks(x); axB.set_xticklabels([f"≥{t}" for t in ths],fontsize=9); axB.set_xlabel("graph PAV-ratio 'present' threshold",fontsize=9.2); axB.set_ylabel("% of strain×locus events",fontsize=9.5); axB.set_ylim(0,105)
axB.legend(fontsize=7.3,frameon=False,loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=1); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Silencing dominates, genetic loss is rare\n— robust across thresholds",fontsize=9.6,fontweight="bold",loc="left",y=1.18)
# C: per liftover class, fraction the graph rescues to shared (n>=2)
gn=(gpr>=0.5).astype(int).sum(1); resc=[]
for c in ORD:
    sub=comp.liftover_class==c; gns=gn[sub.values]; resc.append([100*(gns>=2).mean(),100*((gns==1)|(gns==0)).mean()])
resc=np.array(resc)
axC.bar(range(3),resc[:,0],color="#1B7837",label="sequence-SHARED in graph (≥2 strains)",edgecolor="white")
axC.bar(range(3),resc[:,1],bottom=resc[:,0],color="#cccccc",label="graph private/absent",edgecolor="white")
for i,c in enumerate(ORD): axC.text(i,resc[i,0]-5,f"{resc[i,0]:.0f}%",ha="center",va="top",fontsize=10,color="white",fontweight="bold")
axC.set_xticks(range(3)); axC.set_xticklabels([f"liftover\n{c}" for c in ORD],fontsize=8.5); axC.set_ylabel("% of loci",fontsize=9.5); axC.set_ylim(0,108)
axC.legend(fontsize=7.4,frameon=False,loc="lower center"); axC.spines[["top","right"]].set_visible(False)
axC.text(2,50,"99% of liftover-PRIVATE\nloci are sequence-shared",ha="center",fontsize=7.6,color="#1B7837",fontweight="bold")
axC.set_title("C  The graph rescues liftover-'private' clusters\nas genetically SHARED (epigenetically silenced)",fontsize=9.6,fontweight="bold",loc="left")
# D: model schematic
axD.axis("off"); axD.set_xlim(0,1); axD.set_ylim(0,1)
axD.text(0.5,0.95,"Two questions, two answers",ha="center",fontsize=10,fontweight="bold")
axD.text(0.5,0.84,"HAL-liftover  →  is there a piRNA CLUSTER here?  (regulatory)",ha="center",fontsize=8.2,color="#E8852B",fontweight="bold")
axD.text(0.5,0.77,"odgi graph  →  is the SEQUENCE here?  (genetic)",ha="center",fontsize=8.2,color="#1B7837",fontweight="bold")
axD.plot([0.08,0.92],[0.70,0.70],color="#ddd",lw=1)
axD.text(0.5,0.60,"median strains:  graph SEQUENCE in 16   vs   liftover CLUSTER in 4",ha="center",fontsize=8.4,color="#222",fontweight="bold")
axD.text(0.5,0.50,"→ cross-strain piRNA-cluster variation is mostly EPIGENETIC SILENCING\n(sequence conserved, cluster not expressed), NOT genetic loss (only 3–8%).",ha="center",fontsize=7.9,color="#444")
axD.text(0.5,0.37,"Liftover (one linear reference) CONFLATES silencing & loss;\nthe pangenome graph SEPARATES them.",ha="center",fontsize=8.1,color="#B2182B",fontweight="bold")
axD.text(0.5,0.22,"the genetically-NOVEL minority (clusters absent from GRCm39) =\nyoung L1 / IAP TE insertions  →  Fig 2",ha="center",fontsize=7.6,color="#555",style="italic")
axD.text(0.5,0.06,f"42,384 master piRNA-cluster loci · 16 strains + GRCm39 · odgi pav vs halLiftover",ha="center",fontsize=6.8,color="#888")
axD.set_title("D  Biology: silencing vs loss",fontsize=9.6,fontweight="bold",loc="left")
fig.suptitle("Pangenome graph vs reference liftover: strain-variable piRNA clusters are mostly EPIGENETICALLY SILENCED, not genetically lost",fontsize=11.6,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{T}/figures/Fig_graph_vs_liftover_pav"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out)
