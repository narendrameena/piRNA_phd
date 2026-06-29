#!/usr/bin/env python3
"""THEME 22 figure 2 — NON-REFERENCE piRNA clusters, CONFOUNDER-CHECKED. (A) per-strain counts (1,393; 93% TE);
(B) multimapping-CORRECTED expression: all-primary FPM is inflated (non-ref TE-rich, 25% multimap) but on UNIQUE reads
non-ref ~ reference -> genuinely expressed, not 'higher'; (C) RIGOROUS evolution test on the DIRECT VCF TE-insertion
burden: non-ref cluster count tracks genome-wide non-ref TE insertions (Spearman), robust to multimapping (unique-read
share rho=0.51) and total-output (partial r=0.50); (D) young LTR/ERVK + LINE-1, BioMNI-verified arms race + confounders."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
T=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav"; D=f"{T}/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
te=pd.read_csv(f"{D}/nonref_te_summary.csv"); fam=pd.read_csv(f"{D}/nonref_te_families.csv")
mm=pd.read_csv(f"{D}/multimapping_test.csv"); burden=json.load(open(f"{D}/te_insertion_burden.json"))
ev=mm.merge(te[["strain","n_nonref"]],on="strain"); ev["te_burden"]=ev.strain.map(burden); ev["wild"]=ev.strain.isin(WILD)
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(13.2,9.9),dpi=300)
# A: per-strain non-ref counts
a=te.merge(ev[["strain","wild"]],on="strain").sort_values("n_nonref")
axA.barh(range(len(a)),a.n_nonref,color=["#C0392B" if w else "#4393C3" for w in a.wild],edgecolor="white")
axA.set_yticks(range(len(a))); axA.set_yticklabels([s.replace("_","/") for s in a.strain],fontsize=6.8)
axA.set_xlabel("non-reference piRNA clusters (absent from GRCm39)",fontsize=9.2); axA.spines[["top","right"]].set_visible(False)
axA.legend(handles=[Patch(fc="#C0392B",label="wild-derived"),Patch(fc="#4393C3",label="classical")],fontsize=7,frameon=False,loc="lower right")
axA.text(0.97,0.05,"1,393 total · 93% overlap a TE",transform=axA.transAxes,ha="right",fontsize=7.6,color="#444",style="italic",fontweight="bold")
axA.set_title("A  Strain-specific clusters the reference lacks\n(0.4% of all clusters; TE-insertion-driven)",fontsize=9.4,fontweight="bold",loc="left")
# B: multimapping-corrected expression (per-cluster medians, all vs unique)
al_nr=[];al_rf=[];uq_nr=[];uq_rf=[];mm_nr=[];mm_rf=[]
for X in S:
    f=pd.read_csv(f"{CP}/{X}.clusters_fpm.bed",sep="\t",header=None,names=["c","s","e","allF","uniqF","st","tp"],dtype={"c":str})
    nr=pd.read_csv(f"{D}/nonref/{X}.nonref.bed",sep="\t",header=None,names=["c","s","e","id"],dtype={"c":str}); nrset=set(zip(nr.c.astype(str),nr.s,nr.e))
    f["nr"]=[(c,s,e) in nrset for c,s,e in zip(f.c.astype(str),f.s,f.e)]
    cl=f.groupby(["c","s","e","nr"],as_index=False).agg(allF=("allF","sum"),uniqF=("uniqF","sum")); cl["mm"]=1-cl.uniqF/cl.allF.replace(0,np.nan)
    al_nr+=list(cl[cl.nr].allF);al_rf+=list(cl[~cl.nr].allF);uq_nr+=list(cl[cl.nr].uniqF);uq_rf+=list(cl[~cl.nr].uniqF);mm_nr+=list(cl[cl.nr].mm.dropna());mm_rf+=list(cl[~cl.nr].mm.dropna())
med=[np.median(al_rf),np.median(al_nr),np.median(uq_rf),np.median(uq_nr)]
x=[0,0.9,2.3,3.2]; cols=["#cccccc","#1B7837","#cccccc","#1B7837"]
axB.bar(x,med,width=0.8,color=cols,edgecolor="white")
for xi,v in zip(x,med): axB.text(xi,v+0.3,f"{v:.1f}",ha="center",fontsize=8,fontweight="bold")
axB.set_xticks([0.45,2.75]); axB.set_xticklabels(["all-primary FPM\n(multimapping incl.)","UNIQUE-read FPM\n(multimapping removed)"],fontsize=8.3)
axB.set_ylabel("median cluster expression (FPM)",fontsize=9.2); axB.spines[["top","right"]].set_visible(False); axB.set_ylim(0,15)
axB.legend(handles=[Patch(fc="#cccccc",label="reference"),Patch(fc="#1B7837",label="non-reference")],fontsize=7.4,frameon=False,loc="upper right")
axB.text(0.5,0.62,f"non-ref are 25% multimapping (vs 0.2% ref): all-primary\nlooks higher, but on UNIQUE reads non-ref ≈ reference\n(genuinely expressed, NOT inflated 'higher')",transform=axB.transAxes,ha="center",fontsize=7.0,color="#B2182B",fontweight="bold")
axB.set_title("B  Genuinely expressed — once multimapping is removed",fontsize=9.4,fontweight="bold",loc="left")
# C: evolution on DIRECT VCF TE-insertion burden
xc=ev.te_burden/1000; yc=ev.n_nonref; cc=["#C0392B" if w else "#4393C3" for w in ev.wild]
axC.scatter(xc,yc,c=cc,s=48,edgecolor="white",lw=0.5,zorder=3)
z=np.polyfit(xc,yc,1); xx=np.linspace(xc.min(),xc.max(),20); axC.plot(xx,np.polyval(z,xx),color="#888",lw=1.2,ls="--",zorder=2)
rho,p=spearmanr(ev.te_burden,ev.n_nonref)
# robustness, computed live so the annotations are reproducible: unique-read-share rho + Pearson partial r vs total output
_ru,_pu=spearmanr(ev.te_burden,mm.set_index("strain").loc[ev.strain,"nr_pct_uniq"].values)
_cf=pd.read_csv(f"{D}/confounding.csv"); _cf["burden"]=_cf.strain.map(burden)
_rn=_cf.n_nonref.values-np.polyval(np.polyfit(_cf.total_fpm.values,_cf.n_nonref.values.astype(float),1),_cf.total_fpm.values)
_rb=_cf.burden.values.astype(float)-np.polyval(np.polyfit(_cf.total_fpm.values,_cf.burden.values.astype(float),1),_cf.total_fpm.values)
_partial=float(np.corrcoef(_rn,_rb)[0,1])
for s,xv,yv in zip(ev.strain,xc,yc):
    if s in WILD or xv>20: axC.annotate(s.split("_")[0],(xv,yv),fontsize=5.7,color="#555",xytext=(3,2),textcoords="offset points")
axC.set_xlabel("genome-wide non-reference TE insertions (×1000, from VCF)",fontsize=8.8); axC.set_ylabel("non-reference piRNA clusters",fontsize=9.2); axC.spines[["top","right"]].set_visible(False)
axC.text(0.50,0.20,f"Spearman ρ = {rho:+.2f}, p = {p:.3f}\nrobust: unique-read share ρ={_ru:.2f} (p={_pu:.3f})\n+ partial r={_partial:.2f} vs total output",transform=axC.transAxes,va="top",fontsize=7.6,color="#B2182B",fontweight="bold")
axC.set_title("C  piRNA-cluster novelty tracks the DIRECT TE-insertion\nburden across strains (multimapping-independent)",fontsize=9.4,fontweight="bold",loc="left")
# D: TE families + biology + confounders
fl=fam[~fam.TE_family.isin(["Simple_repeat","Low_complexity","Unknown","Satellite","tRNA","rRNA","snRNA"])].head(7)[::-1]
TC={"LTR":"#6a3d9a","LINE":"#E69F00","SINE":"#33a02c","DNA":"#b15928"}
axD.barh(range(len(fl)),fl.n,color=[TC.get(f.split("/")[0],"#888") for f in fl.TE_family],edgecolor="white")
axD.set_yticks(range(len(fl))); axD.set_yticklabels(fl.TE_family,fontsize=6.8)
axD.set_xlabel("annotations at non-reference piRNA clusters",fontsize=8.8); axD.spines[["top","right"]].set_visible(False)
axD.legend(handles=[Patch(fc=TC[k],label=k) for k in ["LTR","LINE","SINE"]],fontsize=6.6,frameon=False,loc="lower right",title="TE class",title_fontsize=6.6)
axD.set_title("D  Young LTR/ERVK + LINE-1 — the piRNA–TE arms race",fontsize=9.4,fontweight="bold",loc="left")
axD.text(0.5,-0.28,"A young active TE (ERVK/L1) inserts AFTER strains diverge → host mounts a piRNA response → a NEW, strain-specific,\nwell-expressed piRNA cluster absent from the reference (BioMNI-verified: Aravin 2007; Girard 2010; Gainetdinov 2015).\nCONFOUNDERS CHECKED: multimapping (corrected ✓) · total output (partial-controlled ✓) · lift-artifact (not divergence-differential ✓).\nCAVEAT: effect is wild-vs-classical driven (within-classical n.s.).  GRAPH-VALIDATED: the ~8% cross-strain-shared tail is GENUINELY absent\nfrom GRCm39 (odgi pav — GRCm39 covers frame+lifted controls ≈1.0, these loci 0.0); the non-reference catch is genuine, not a liftover artifact.",transform=axD.transAxes,ha="center",va="top",fontsize=6.2,color="#444",style="italic")
fig.suptitle("Non-reference piRNA clusters: the young, TE-insertion-driven leading edge of piRNA-cluster evolution (confounder-checked across 16 strains)",fontsize=10.4,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0.075,1,0.96])
# --- per-figure SourceData: every plotted value (A counts, B expression, C evolution, D TE families) ---
sd=[]
for _,r in a.iterrows(): sd.append(dict(panel="A_counts",key=r.strain,group=("wild" if r.wild else "classical"),metric="n_nonref_clusters",value=int(r.n_nonref)))
sd.append(dict(panel="A_counts",key="(total)",group="",metric="n_nonref_total",value=int(a.n_nonref.sum())))
sd.append(dict(panel="A_counts",key="(total)",group="",metric="pct_overlap_TE",value=round(100*te.n_TE_overlap.sum()/te.n_nonref.sum(),1)))
for lab,v in zip(["reference_all","nonref_all","reference_unique","nonref_unique"],med): sd.append(dict(panel="B_expression",key=lab,group="",metric="median_FPM",value=round(float(v),2)))
sd.append(dict(panel="B_expression",key="nonref",group="",metric="median_multimap_frac",value=round(float(np.median(mm_nr)),3)))
sd.append(dict(panel="B_expression",key="reference",group="",metric="median_multimap_frac",value=round(float(np.median(mm_rf)),4)))
for _,r in ev.iterrows(): sd.append(dict(panel="C_evolution",key=r.strain,group=("wild" if r.wild else "classical"),metric="te_insertion_burden_VCF",value=int(r.te_burden)))
for _,r in ev.iterrows(): sd.append(dict(panel="C_evolution",key=r.strain,group=("wild" if r.wild else "classical"),metric="n_nonref_clusters",value=int(r.n_nonref)))
sd.append(dict(panel="C_evolution",key="(spearman_main)",group=f"p={p:.3f}",metric="rho_burden_vs_count",value=round(float(rho),3)))
sd.append(dict(panel="C_evolution",key="(robustness_uniqueread_share)",group=f"p={_pu:.3f}",metric="rho_burden_vs_count",value=round(_ru,3)))
sd.append(dict(panel="C_evolution",key="(robustness_partial_vs_total_output)",group="",metric="partial_r",value=round(_partial,3)))
for _,r in fl.iterrows(): sd.append(dict(panel="D_TEfamily",key=str(r.TE_family),group="",metric="n_annotations",value=int(r.n)))
pd.DataFrame(sd).to_csv(f"{D}/source_data/SourceData_Fig_nonreference_clusters.csv",index=False)
out=f"{T}/figures/Fig_nonreference_clusters"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out,"| total",int(a.n_nonref.sum()),"| spearman",round(float(rho),3))
