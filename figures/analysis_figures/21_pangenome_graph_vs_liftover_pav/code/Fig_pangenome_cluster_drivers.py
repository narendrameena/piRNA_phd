#!/usr/bin/env python3
"""THEME 21 Fig 2 — biology of the strain-variable piRNA clusters: (A) the genetically-variable loci are driven by
young TE insertions (L1MdTf/L1MdA + IAP LTR); (B) developmental class x gene context — pachytene clusters are
intergenic (A-MYB), pre-pachytene/hybrid are 3'UTR-enriched (mRNA-derived); (C) antisense pseudogene fragments in
pachytene clusters (candidate mRNA-targeting piRNA sources); (D) piC-DoG transcriptional readthrough from whole-testis
bulk RNA (~65% of downstream clusters). Independent of the graph PAV; uses the precomputed VCF/RM, GFF3, bulk RNA."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
T="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav"; D=f"{T}/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(13.5,9.7),dpi=300)
CLS_COL={"L1":"#E69F00","LTR":"#6a3d9a","SINE":"#33a02c","LINE":"#E69F00"}
# A: TE drivers of the variable loci
te=pd.read_csv(f"{D}/sv_te_hits.tsv",sep="\t",header=None); te["cls"]=te[1].str.split(":").str[0]
top=te[1].value_counts().head(9)[::-1]
cols=[CLS_COL.get(f.split(":")[0],"#888") for f in top.index]
axA.barh(range(len(top)),top.values,color=cols,edgecolor="white")
axA.set_yticks(range(len(top))); axA.set_yticklabels([f.replace(":","\n") for f in top.index],fontsize=7.2)
for i,v in enumerate(top.values): axA.text(v+max(top.values)*0.01,i,f"{v:,}",va="center",fontsize=7.5,fontweight="bold")
axA.set_xlabel("TE-insertion SVs at variable piRNA-cluster loci",fontsize=9.3); axA.set_xlim(0,max(top.values)*1.16); axA.spines[["top","right"]].set_visible(False)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(fc=CLS_COL[k],label=k) for k in ["L1","LTR","SINE"]],fontsize=7.5,frameon=False,loc="lower right",title="TE class",title_fontsize=7.5)
axA.set_title("A  Strain-variable clusters arise at young TE insertions\n(L1MdTf/L1MdA LINE-1 + IAP LTR dominate)",fontsize=9.5,fontweight="bold",loc="left")
# B: dev class x gene context (row-normalised stacked)
dv=pd.read_csv(f"{D}/master_devclass.tsv",sep="\t",header=None,names=["locus_id","devclass","pattern"],dtype={"pattern":str})
gc=pd.read_csv(f"{D}/master_genecontext.tsv",sep="\t",header=None,names=["locus_id","context"])
mm=dv.merge(gc,on="locus_id"); CTX=["3UTR","CDS","gene_body","downstream","intergenic"]
CTXCOL={"3UTR":"#4393C3","CDS":"#80b1d3","gene_body":"#bbbbbb","downstream":"#E8852B","intergenic":"#B2182B"}
DEV=["pachytene","pre_pachytene","hybrid","P12.5_only","fetal"]
ct=pd.crosstab(mm.devclass,mm.context,normalize="index").reindex(index=DEV,columns=CTX,fill_value=0)*100
bot=np.zeros(len(DEV))
for c in CTX:
    axB.bar(range(len(DEV)),ct[c].values,bottom=bot,color=CTXCOL[c],label=c,edgecolor="white",width=0.72); bot+=ct[c].values
axB.set_xticks(range(len(DEV))); axB.set_xticklabels([d.replace("_","\n") for d in DEV],fontsize=8); axB.set_ylabel("% of loci (gene context)",fontsize=9.3); axB.set_ylim(0,100)
axB.legend(fontsize=6.8,frameon=False,ncol=5,loc="lower center",bbox_to_anchor=(0.5,1.0)); axB.spines[["top","right"]].set_visible(False)
axB.text(0,100-ct.loc["pachytene","intergenic"]/2,f"{ct.loc['pachytene','intergenic']:.0f}%\nintergenic",ha="center",va="center",fontsize=6.6,color="white",fontweight="bold")
axB.set_title("B  Pachytene clusters = intergenic (A-MYB); pre-pachytene\n/hybrid = 3'UTR-enriched (mRNA-derived)",fontsize=9.5,fontweight="bold",loc="left",y=1.1)
# C: pseudogene fragments (antisense = mRNA-targeting)
pgf=pd.read_csv(f"{D}/pgf_hits.tsv",sep="\t",header=None); ori=pgf.iloc[:,-1].value_counts()
axC.bar([0,1],[ori.get("antisense",0),ori.get("sense",0)],color=["#C0392B","#bbbbbb"],width=0.6,edgecolor="white")
for i,v in enumerate([ori.get("antisense",0),ori.get("sense",0)]): axC.text(i,v+ori.max()*0.02,f"{v:,}",ha="center",fontsize=10,fontweight="bold")
axC.set_xticks([0,1]); axC.set_xticklabels(["ANTISENSE\n(mRNA-targeting)","sense"],fontsize=8.5); axC.set_ylabel("pseudogene-fragment hits",fontsize=9.3); axC.set_ylim(0,ori.max()*1.2); axC.spines[["top","right"]].set_visible(False)
axC.text(0.5,0.86,f"{ori.get('antisense',0):,} antisense PGFs in pachytene clusters\n= candidate trans-acting mRNA-silencing piRNA sources",transform=axC.transAxes,ha="center",fontsize=7.3,color="#444",style="italic")
axC.set_title("C  Pachytene clusters carry ANTISENSE pseudogene\nfragments (mRNA-targeting; Loubalova 2025)",fontsize=9.5,fontweight="bold",loc="left")
# D: piC-DoG readthrough per strain
rt=[]
for s in S:
    try:
        df=pd.read_csv(f"{D}/readthrough_{s}.tsv",sep="\t"); t=len(df); y=(df.readthrough=="YES").sum(); rt.append((s,100*y/t if t else 0))
    except Exception: rt.append((s,0))
rt=sorted(rt,key=lambda x:x[1]); xs=[r[1] for r in rt]
axD.barh(range(len(rt)),xs,color=["#C0392B" if r[0] in WILD else "#4393C3" for r in rt],edgecolor="white")
axD.set_yticks(range(len(rt))); axD.set_yticklabels([r[0].replace("_","/") for r in rt],fontsize=6.6)
ov=100*323715/496813
axD.axvline(ov,color="#222",lw=1.2,ls="--"); axD.text(ov+1,0.5,f"overall {ov:.0f}%",fontsize=7.5,color="#222",fontweight="bold")
axD.set_xlabel("% of downstream clusters with piC-DoG readthrough",fontsize=9.2); axD.set_xlim(0,100); axD.spines[["top","right"]].set_visible(False)
axD.legend(handles=[Patch(fc="#C0392B",label="wild-derived"),Patch(fc="#4393C3",label="classical")],fontsize=7,frameon=False,loc="lower right")
axD.set_title("D  piC-DoG readthrough: whole-testis RNA reads through\ngene 3' ends into ~65% of downstream clusters",fontsize=9.5,fontweight="bold",loc="left")
fig.suptitle("Drivers of strain-variable piRNA clusters: TE insertions (genetic novelty), gene context, antisense pseudogene fragments, and piC-DoG readthrough",fontsize=11,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{T}/figures/Fig_pangenome_cluster_drivers"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out)
