#!/usr/bin/env python3
"""THEME 18 — biological difference between the two kinds of within-tp genuinely-unique piRNAs, and a CONTROL for
whether the heterochronic signal is real or a multimapping artifact.
  HETERO     = 'stage-shifted (heterochronic)': unique in its home stage but the SAME sequence is expressed in
               another strain at a DIFFERENT stage (re-deployed across development).
  RESTRICTED = 'conserved-but-silent' + 'strain-private': unique within stage AND not expressed at any other stage.
Findings (each TESTED): heterochronic piRNAs are MORE TE-derived AND MORE multimapping -- and these are the SAME
phenomenon (TEs are multicopy: TE-fraction rises monotonically with NH). BUT it is NOT a pure multimapping artifact:
79% of heterochronic piRNAs map to a SINGLE locus (NH=1) and stay TE-enriched -> genuine single-locus re-expression,
with a TE-multicopy minority. Heterochronic are also more strongly stage-peaked (higher log2FC).
Data: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz + cand_self16/*.bam (NH) + sense_antisense percand."""
import warnings, os, subprocess; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
ST="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/samtools"; WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
HET="unique: stage-shifted (heterochronic)"; RES=["unique: conserved-but-silent","unique: strain-private locus"]
CH="#009E73"; CR="#9aa0a6"
d=pd.read_csv(f"{U}/deseq16_lenfilt/deseq_stagepeak_classified.csv.gz")
d=d[d.klass.isin([HET]+RES)].copy(); d["grp"]=np.where(d.klass==HET,"HETERO","RESTRICTED")
cache=f"{T}/data/source_data/SourceData_Fig_heterochronic_vs_restricted.csv.gz"
if os.path.exists(cache):
    d=pd.read_csv(cache)
else:
    need=set(d.cand_id); nh={}
    for X in sorted(d.strain.unique()):
        out=subprocess.run([ST,"view",f"{U}/cand_self16/{X}.cand_self16.bam"],capture_output=True,text=True,timeout=900).stdout
        for ln in out.splitlines():
            f=ln.split("\t")
            if len(f)<11 or int(f[1]) not in (0,16) or f[0] not in need: continue
            for t in f[11:]:
                if t.startswith("NH:i:"): nh[f[0]]=int(t[5:]); break
    d["nh"]=d.cand_id.map(nh); d=d.dropna(subset=["nh"]); d["nh"]=d.nh.astype(int)
    oc=pd.read_csv(f"{U}/sense_antisense/SourceData_sense_antisense16_percand.csv.gz")
    d=d.merge(oc[["id","family"]].drop_duplicates("id"),left_on="cand_id",right_on="id",how="left")
    d["te"]=d.family.notna(); d["classical"]=~d.strain.isin(WILD)
    d[["cand_id","grp","nh","te","classical","log2FC"]].to_csv(cache,index=False)
G=["HETERO","RESTRICTED"]; H=d[d.grp=="HETERO"]; R=d[d.grp=="RESTRICTED"]; d["multi"]=d.nh>1
pTE=fisher_exact(pd.crosstab(d.grp,d.te).values)[1]; pMM=fisher_exact(pd.crosstab(d.grp,d.multi).values)[1]
pFC=mannwhitneyu(H.log2FC,R.log2FC)[1]
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(12.8,9.6),dpi=300)
# A: both TE-derived and multimapping elevated in heterochronic
xb=np.arange(2); bw=0.36
teV=[H.te.mean()*100,R.te.mean()*100]; mmV=[(H.nh>1).mean()*100,(R.nh>1).mean()*100]
axA.bar(xb-bw/2,teV,bw,color="#1b9e77",label="TE-derived"); axA.bar(xb+bw/2,mmV,bw,color="#d95f02",label="multimapping (NH>1)")
for i in range(2):
    axA.text(i-bw/2,teV[i]+0.8,f"{teV[i]:.0f}%",ha="center",fontsize=8,fontweight="bold")
    axA.text(i+bw/2,mmV[i]+0.8,f"{mmV[i]:.0f}%",ha="center",fontsize=8,fontweight="bold")
axA.set_xticks(xb); axA.set_xticklabels(["heterochronic","stage-restricted"],fontsize=9); axA.set_ylabel("% of unique piRNAs",fontsize=9.5)
axA.set_ylim(0,max(teV)*1.3); axA.legend(fontsize=8,frameon=False,loc="upper right"); axA.spines[["top","right"]].set_visible(False)
axA.set_title(f"A  Heterochronic piRNAs are MORE TE-derived AND MORE multimapping\n(TE Fisher p={pTE:.0e}; multimap Fisher p={pMM:.0e})",fontsize=9.5,fontweight="bold",loc="left")
# B: TE and multimapping are the same phenomenon (TEs are multicopy)
d["nhbin"]=pd.cut(d.nh,[0,1,2,10,100,1e9],labels=["1","2","3-10","11-100",">100"])
teb=d.groupby("nhbin").te.mean()*100; nb=d.groupby("nhbin").size()
axB.bar(range(len(teb)),teb.values,color=plt.cm.YlOrRd(np.linspace(0.3,0.95,len(teb))),edgecolor="white")
for i,(v,n) in enumerate(zip(teb.values,nb.values)): axB.text(i,v+1.5,f"{v:.0f}%\n(n={n:,})",ha="center",va="bottom",fontsize=7.5)
axB.set_xticks(range(len(teb))); axB.set_xticklabels(teb.index); axB.set_xlabel("number of genomic loci a piRNA maps to (NH)",fontsize=9.5)
axB.set_ylabel("% TE-derived",fontsize=9.5); axB.set_ylim(0,112); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Multimapping IS TE: TE-fraction rises 31% -> 100% with NH\n(TEs are multicopy, so the two are inseparable)",fontsize=9.5,fontweight="bold",loc="left")
# C: disentangle -- not a pure multimapping artifact
u=H[H.nh==1]; mm=H[H.nh>1]; fu=100*len(u)/len(H); fm=100*len(mm)/len(H)
teU_h=H[H.nh==1].te.mean()*100; teU_r=R[R.nh==1].te.mean()*100; pU=fisher_exact(pd.crosstab(d[d.nh==1].grp,d[d.nh==1].te).values)[1]
axC.bar([0],[fu],0.6,color=CH,label="single-locus (NH=1): genuine re-expression")
axC.bar([0],[fm],0.6,bottom=[fu],color="#d95f02",label="multicopy (NH>1): TE-driven")
axC.text(0,fu/2,f"{fu:.0f}%",ha="center",va="center",color="white",fontsize=11,fontweight="bold")
axC.text(0,fu+fm/2,f"{fm:.0f}%",ha="center",va="center",color="white",fontsize=9,fontweight="bold")
axC.set_xlim(-0.8,1.6); axC.set_xticks([0]); axC.set_xticklabels(["heterochronic\npiRNAs"],fontsize=9); axC.set_ylabel("% of heterochronic piRNAs",fontsize=9.5)
axC.legend(fontsize=7.6,frameon=False,loc="upper right",bbox_to_anchor=(1.02,1.0)); axC.set_ylim(0,108); axC.spines[["top","right"]].set_visible(False)
axC.text(0.5,0.40,f"among SINGLE-locus piRNAs,\nheterochronic still TE-enriched\n({teU_h:.0f}% vs {teU_r:.0f}%, p={pU:.0e})\n-> real, not a mapping artifact",transform=axC.transAxes,fontsize=8,color="#333",fontweight="bold",va="center")
axC.set_title("C  79% map to ONE locus -> mostly genuine heterochrony\n(only 21% are the TE-multicopy multimappers)",fontsize=9.5,fontweight="bold",loc="left")
# D: orthogonal biology -- more strongly stage-peaked
parts=axD.violinplot([H.log2FC,R.log2FC],showmedians=True,widths=0.85)
for b,c in zip(parts["bodies"],[CH,CR]): b.set_facecolor(c); b.set_alpha(0.6)
for k in ("cmedians","cbars","cmins","cmaxes"): parts[k].set_color("#333")
axD.set_xticks([1,2]); axD.set_xticklabels([f"heterochronic\n(n={len(H):,})",f"stage-restricted\n(n={len(R):,})"],fontsize=9)
axD.set_ylabel("stage-peak strength  log2FC (DESeq2)",fontsize=9.5); axD.spines[["top","right"]].set_visible(False)
axD.text(0.5,0.96,f"medians {H.log2FC.median():.1f} vs {R.log2FC.median():.1f}; MWU p={pFC:.0e}",transform=axD.transAxes,ha="center",va="top",fontsize=8.5,fontweight="bold")
axD.set_title("D  Beyond TE/mapping: heterochronic are also\nmore strongly stage-peaked (higher fold-induction)",fontsize=9.5,fontweight="bold",loc="left")
fig.suptitle("Heterochronic vs stage-restricted unique piRNAs — TE-derived multicopy is involved, but 79% are single-locus genuine re-expression",fontsize=12,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{T}/figures/Fig_heterochronic_vs_restricted"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print(f"wrote {out}")
print(f"TE% {teV[0]:.0f}/{teV[1]:.0f} p={pTE:.0e}; multimap% {mmV[0]:.0f}/{mmV[1]:.0f} p={pMM:.0e}; NH=1 share heterochronic {fu:.0f}%; log2FC {H.log2FC.median():.2f}/{R.log2FC.median():.2f} p={pFC:.0e}")
