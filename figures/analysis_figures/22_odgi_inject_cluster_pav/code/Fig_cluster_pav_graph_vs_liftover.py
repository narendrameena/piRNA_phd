#!/usr/bin/env python3
"""THEME 22 figure — graph-inject cluster-PAV (odgi inject + odgi pav, reference-free cluster placement) vs HAL-liftover
cluster-PAV, at 42,384 master loci. Both ask 'which strains have a piRNA CLUSTER here', via two independent methods.
Result: 99% per-strain-locus agreement -> the cluster-PAV is method-robust (cross-validated); the graph is slightly
more conservative (inject drops 0.4% boundary clusters + coverage threshold), never spuriously adds (graph-only=0)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
T21=f"{B}/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav/data"
T=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav"; D=f"{T}/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
lift=pd.read_csv(f"{T21}/liftover_pav_matrix.tsv",sep="\t",dtype={"chrom":str})
g=pd.read_csv(f"{D}/graph_cluster_pav_matrix.tsv",sep="\t").rename(columns={"name":"locus_id"})
m=lift.merge(g[["locus_id"]+S].rename(columns={s:"g_"+s for s in S}),on="locus_id")
lp=m[S].values.astype(int); gr=m[["g_"+s for s in S]].values; gp=(gr>=0.1).astype(int)
comp=pd.read_csv(f"{D}/cluster_pav_comparison.tsv",sep="\t")
ORD=["core","dispensable","private"]
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(13,9.6),dpi=300)
# A: confusion heatmap
conf=pd.crosstab(comp.liftover_class,comp.graph_class).reindex(index=ORD,columns=ORD+["absent"],fill_value=0)
axA.imshow(conf.values,cmap="Greens",aspect="auto",norm=matplotlib.colors.LogNorm(vmin=1,vmax=conf.values.max()))
axA.set_xticks(range(4)); axA.set_xticklabels(["core","disp.","private","absent"],fontsize=9); axA.set_yticks(range(3)); axA.set_yticklabels(ORD,fontsize=9)
for i in range(3):
    for j in range(4):
        v=conf.values[i,j]; axA.text(j,i,f"{v:,}",ha="center",va="center",fontsize=9,fontweight=("bold" if i==j else "normal"),color=("white" if v>conf.values.max()*0.2 else "#333"))
axA.set_xlabel("GRAPH-inject cluster class",fontsize=9.5); axA.set_ylabel("HAL-liftover cluster class",fontsize=9.5)
axA.set_title("A  The two cluster-PAV methods agree (strong diagonal):\ngraph-inject cross-validates HAL-liftover",fontsize=9.6,fontweight="bold",loc="left")
# B: per-strain agreement
ag=[100*(lp[:,i]==gp[:,i]).mean() for i in range(len(S))]; idx=np.argsort(ag)
axB.barh(range(len(S)),[ag[i] for i in idx],color=["#C0392B" if S[i] in WILD else "#4393C3" for i in idx],edgecolor="white")
axB.set_yticks(range(len(S))); axB.set_yticklabels([S[i].replace("_","/") for i in idx],fontsize=6.8)
axB.set_xlim(90,100.5); axB.set_xlabel("% strain×locus agreement (graph vs liftover)",fontsize=9); axB.spines[["top","right"]].set_visible(False)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(fc="#C0392B",label="wild-derived"),Patch(fc="#4393C3",label="classical")],fontsize=7,frameon=False,loc="lower left")
axB.set_title("B  Agreement is 99–100% for EVERY strain,\nincluding the divergent wild strains",fontsize=9.6,fontweight="bold",loc="left")
# C: threshold robustness — class counts graph vs liftover
ths=[0.05,0.1,0.2]; gcore=[];gdisp=[];gpriv=[]
for th in ths:
    gpx=(gr>=th).astype(int).sum(1); gcore.append((gpx>=16).sum()); gdisp.append(((gpx>=2)&(gpx<16)).sum()); gpriv.append((gpx==1).sum())
x=np.arange(len(ths)); w=0.6
lc=[(lift[S].values.sum(1)>=16).sum(),((lift[S].values.sum(1)>=2)&(lift[S].values.sum(1)<16)).sum(),(lift[S].values.sum(1)==1).sum()]
axC.bar(x,gcore,w,color="#1B7837",label="core",edgecolor="white"); axC.bar(x,gdisp,w,bottom=gcore,color="#80b1d3",label="dispensable",edgecolor="white"); axC.bar(x,gpriv,w,bottom=np.array(gcore)+np.array(gdisp),color="#E8852B",label="private",edgecolor="white")
axC.axhline(lc[0],color="#1B7837",ls="--",lw=1); axC.text(2.4,lc[0],f"liftover core {lc[0]}",fontsize=6.5,va="center",color="#1B7837")
axC.set_xticks(x); axC.set_xticklabels([f"≥{t}" for t in ths],fontsize=9); axC.set_xlabel("graph cluster-coverage threshold",fontsize=9.2); axC.set_ylabel("# loci",fontsize=9.3)
axC.legend(fontsize=7.5,frameon=False,loc="upper right"); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  Graph class counts ≈ liftover (5983/22866/13535),\nrobust to the coverage threshold",fontsize=9.6,fontweight="bold",loc="left")
# D: interpretation
axD.axis("off"); axD.set_xlim(0,1); axD.set_ylim(0,1)
axD.text(0.5,0.95,"Two methods, one answer — cross-validation",ha="center",fontsize=10,fontweight="bold")
axD.text(0.5,0.83,"HAL-liftover  →  project clusters to GRCm39 (coordinate overlap)",ha="center",fontsize=8.1,color="#E8852B")
axD.text(0.5,0.76,"odgi inject  →  place clusters on each strain's own graph path",ha="center",fontsize=8.1,color="#1B7837")
axD.text(0.5,0.74,"             then odgi pav node-coverage at the same loci",ha="center",fontsize=8.1,color="#1B7837")
axD.plot([0.08,0.92],[0.67,0.67],color="#ddd",lw=1)
axD.text(0.5,0.575,"99% per-strain-locus AGREEMENT  →  the cluster conservation\nclassification (core/dispensable/private) is METHOD-ROBUST.",ha="center",fontsize=8.4,color="#222",fontweight="bold")
axD.text(0.5,0.45,"graph is slightly MORE conservative (never adds spurious clusters,\ngraph-only=0; misses ~0.4% that span fragment boundaries).",ha="center",fontsize=7.6,color="#555")
axD.text(0.5,0.31,"Unlike theme 21 (graph SEQUENCE vs liftover CLUSTER = silencing),\nhere BOTH measure cluster presence → they concur.",ha="center",fontsize=7.6,color="#444",style="italic")
axD.text(0.5,0.13,"Scope: compared at the lifted-cluster master loci. Reference-FREE locus\ndefinition (non-ref clusters) needs odgi overlap/untangle, which hit graph-scale limits at 365k paths.",ha="center",fontsize=6.7,color="#888")
axD.set_title("D  Interpretation",fontsize=9.6,fontweight="bold",loc="left")
fig.suptitle("Graph-inject vs HAL-liftover cluster-PAV: two independent methods cross-validate the piRNA-cluster conservation classification (99% agreement)",fontsize=10.8,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{T}/figures/Fig_cluster_pav_graph_vs_liftover"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out)
