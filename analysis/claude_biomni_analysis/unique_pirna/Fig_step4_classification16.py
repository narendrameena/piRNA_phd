#!/usr/bin/env python3
"""16-strain unique-piRNA classification — CANONICAL pangenome-syntenic 4-class (incl. SNP-variant), the
'how genuinely unique are they' view. Per strain (canonical order): counts + composition of expressed-elsewhere
(exact) / SNP-variant (1-3 mm allele expressed at the orthologous locus elsewhere) = NOT unique, vs unique
conserved-but-silent / unique strain-private-locus = GENUINELY unique by expression; bottom panel annotates the
genuinely-unique % per strain. Input = unique16/final_classified_clean.csv.gz (klass5): mm0-clean strain-private; orthologous-locus
presence via the pangenome (halLiftover, one-time projection); SNP-variant = the other strain's allele at the
SAME locus is 1-3 mm from the piRNA AND in its >=2-read pool. (Validated > the STAR pairwise route, which
under-calls SNP-variant by missing expressed syntenic copies.)"""
import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
d=pd.read_csv(f"{U}/unique16/final_classified_clean.csv.gz")   # CANONICAL 5-class (klass5): strain-private filtered to mm0-mappable own-genome loci
CANON=[s for s in STRAIN_ORDER if s in set(d.strain)]
KL=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
LAB=["expressed-elsewhere (not unique)","SNP-variant (allelic — not unique)","low-quality (no mm0 own locus)","unique: conserved-but-silent","unique: strain-private locus (clean)"]
COL=["#9e9e9e","#E69F00","#cdb892","#0072B2","#7a3b9a"]
ct=pd.crosstab(d.strain,d.klass5).reindex(CANON)[KL]; ct.to_csv(f"{PG}/SourceData_step4_classification16.csv")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,(ax1,ax2,ax3,ax4)=plt.subplots(4,1,figsize=(13,19),dpi=300); x=np.arange(len(CANON)); bw=0.16
for i,(k,lab) in enumerate(zip(KL,LAB)):
    ax1.bar(x+(i-2)*bw,ct[k],bw,color=COL[i],label=lab)
    for xi,v in zip(x+(i-2)*bw,ct[k]): ax1.text(xi,v*1.06,f"{int(v):,}",ha="center",va="bottom",fontsize=3.8,rotation=90,color=COL[i])
ax1.set_yscale("log"); ax1.set_ylim(20,ct.values.max()*2.4); ax1.set_xticks(x); ax1.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8)
ax1.set_ylabel("strain-specific piRNA candidates (count, log)",fontsize=9); ax1.legend(fontsize=7,frameon=False,ncol=3,loc="upper left")
ax1.set_title("A  16-strain pangenome-syntenic unique-piRNA classification (mm0-clean strain-private) — counts per strain, canonical order",fontsize=10.4,fontweight="bold",pad=16)
ax1.spines[['top','right']].set_visible(False)
prop=ct.div(ct.sum(1),axis=0); bottom=np.zeros(len(CANON))
for i,(k,lab) in enumerate(zip(KL,LAB)):
    ax2.bar(x,prop[k],0.74,bottom=bottom,color=COL[i],label=lab); bottom+=prop[k].values
gu=ct[["unique: conserved-but-silent","unique: strain-private locus"]].sum(1)/ct.sum(1)*100
ax2.set_xticks(x); ax2.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8); ax2.set_ylim(0,1)
ax2.set_ylabel("class composition (proportion)",fontsize=9); ax2.set_title("B  Class composition per strain (text = % genuinely unique by expression)",fontsize=10,fontweight="bold",pad=14)
for xi,s in zip(x,CANON): ax2.text(xi,1.01,f"{gu[s]:.0f}%",ha="center",va="bottom",fontsize=5.6,color="#333")
ax2.spines[['top','right']].set_visible(False)
# C: developmental timepoint origin of GENUINELY-UNIQUE piRNAs, per strain (strain-wise × timepoint-wise — NO pooling across strains or timepoints)
TPC={"16.5dpc":"#0072B2","12.5dpp":"#009E73","20.5dpp":"#D55E00"}; TPL={"16.5dpc":"E16.5 (prepachytene)","12.5dpp":"P12.5","20.5dpp":"P20.5 (pachytene)"}
gud=d[d.klass5.str.startswith("unique")]
tpc=pd.crosstab(gud.strain,gud.timepoint).reindex(index=CANON).fillna(0)
for t in TPC:
    if t not in tpc.columns: tpc[t]=0.0
tpp=tpc[list(TPC)].div(tpc[list(TPC)].sum(1).replace(0,1),axis=0); bot=np.zeros(len(CANON))
for t in TPC:
    ax3.bar(x,tpp[t].values,0.74,bottom=bot,color=TPC[t],label=TPL[t]); bot+=tpp[t].values
ax3.set_xticks(x); ax3.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8); ax3.set_ylim(0,1)
ax3.set_ylabel("timepoint origin\n(proportion of genuinely-unique)",fontsize=9); ax3.legend(fontsize=7,frameon=False,ncol=3,loc="upper left")
ax3.set_title("C  Developmental timepoint origin of genuinely-unique piRNAs — per strain (strain-wise × timepoint-wise, NO pooling); text = n",fontsize=10,fontweight="bold",pad=14)
for xi,s in zip(x,CANON): ax3.text(xi,1.01,f"{int(tpc.loc[s,list(TPC)].sum()):,}",ha="center",va="bottom",fontsize=4.8,color="#777")
ax3.spines[['top','right']].set_visible(False)
tpc.assign(total=tpc[list(TPC)].sum(1)).to_csv(f"{PG}/SourceData_step4_classification16_tp_origin.csv")
# D: developmental timepoint origin BY CLASS (each class's E16.5/P12.5/P20.5 composition, pooled over 16 strains)
KLAB={"expressed elsewhere (exact)":"expressed-elsewhere","SNP-variant (1-3mm)":"SNP-variant","low-quality: no mm0 own-genome locus":"low-quality","unique: conserved-but-silent":"unique: CBS","unique: strain-private locus":"unique: strain-private"}
tpx=pd.crosstab(d.klass5,d.timepoint).reindex(index=KL)
for t in TPC:
    if t not in tpx.columns: tpx[t]=0.0
tpxp=tpx[list(TPC)].div(tpx[list(TPC)].sum(1).replace(0,1),axis=0); xb4=np.arange(len(KL)); bot4=np.zeros(len(KL))
for t in TPC:
    ax4.bar(xb4,tpxp[t].values,0.6,bottom=bot4,color=TPC[t],label=TPL[t]); bot4+=tpxp[t].values
ax4.set_xticks(xb4); ax4.set_xticklabels([KLAB[k] for k in KL],rotation=18,ha="right",fontsize=8.5); ax4.set_ylim(0,1)
ax4.set_ylabel("timepoint origin\n(proportion)",fontsize=9); ax4.legend(fontsize=7,frameon=False,ncol=3,loc="upper left")
ax4.set_title("D  Developmental timepoint origin BY CLASS (pooled over 16 strains); text = n",fontsize=10,fontweight="bold",pad=14)
for xi,k in zip(xb4,KL): ax4.text(xi,1.012,f"{int(tpx.loc[k,list(TPC)].sum()):,}",ha="center",va="bottom",fontsize=6.5,color="#777")
ax4.spines[['top','right']].set_visible(False)
tpx.assign(total=tpx[list(TPC)].sum(1)).to_csv(f"{PG}/SourceData_step4_classification16_tp_by_class.csv")
fig.text(0.5,0.005,"Pangenome-syntenic expression test: a candidate is NOT unique if another strain expresses the identical (exact) or a ≤3-mm SNP-variant allele at the SAME orthologous locus; GENUINELY unique = orthologous locus present elsewhere but truly silent there (conserved-but-silent) or no orthologous locus at all (strain-private locus). "
  "Orthology via the pangenome (halLiftover, one-time projection); coordinate-based, not sequence-containment. Strain-private = mm0-clean own-genome loci (low-quality tan = mm1-3 error/het-SNP, filtered). Source: unique16/final_classified_clean.csv.gz (klass5).",ha="center",fontsize=6.0,color="#555")
fig.tight_layout(rect=[0,0.02,1,1])
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_step4_classification16.{e}",bbox_inches="tight")
print("wrote Fig_step4_classification16.{png,pdf,svg} + source data"); print(ct.assign(total=ct.sum(1)).to_string())
