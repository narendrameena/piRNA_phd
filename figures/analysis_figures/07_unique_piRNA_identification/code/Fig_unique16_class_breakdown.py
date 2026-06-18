#!/usr/bin/env python3
"""16-strain unique-piRNA classification (pangenome locus-based, final_classified_4class.csv.gz) — the 16-strain
version of the pilot class-count figure. Per strain (CANONICAL thesis order), counts + composition of the FOUR
classes: expressed-elsewhere (not unique), SNP-variant (allelic 1-3 mm variant of a piRNA expressed at the
orthologous locus in another strain -> NOT novel), unique conserved-but-silent (locus present elsewhere but
truly silent there = expression divergence), unique strain-private-locus (novel locus). Wild-derived strains
(WSB/CAST/PWK/SPRET) carry far more unique piRNAs — the divergence signal. Source data (incl. timepoint) saved."""
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")   # ADOPTED ≥2-read absence (5-class klass5; expressed-elsewhere trimmed, genuinely-unique identical)
CANON=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
KL=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
LAB=["expressed-elsewhere (not unique)","SNP-variant (allelic — not novel)","low-quality (no mm0 own-genome locus)","unique: conserved-but-silent (expression)","unique: strain-private locus (clean)"]
COL=["#9e9e9e","#E69F00","#cdb892","#0072B2","#7a3b9a"]   # grey/orange/tan = NOT genuinely unique; blue+purple = genuinely unique (clean loci)
ct=pd.crosstab(d.strain,d.klass5).reindex(CANON)[KL]
pd.crosstab([d.strain,d.timepoint],d.klass5).reindex(KL,axis=1).to_csv(f"{PG}/SourceData_unique16_class_breakdown.csv")
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,(ax1,ax2,ax3)=plt.subplots(3,1,figsize=(13,13.8),dpi=300); x=np.arange(len(CANON)); bw=0.16
for i,k in enumerate(KL):
    ax1.bar(x+(i-2)*bw,ct[k],bw,color=COL[i],label=LAB[i])
    for xi,v in zip(x+(i-2)*bw,ct[k]): ax1.text(xi,v*1.06,f"{int(v):,}",ha="center",va="bottom",fontsize=3.8,rotation=90,color=COL[i])
ax1.set_yscale("log"); ax1.set_ylim(30,ct.values.max()*2.2); ax1.set_xticks(x); ax1.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8)
ax1.set_ylabel("piRNA sequences (count, log scale)",fontsize=9); ax1.legend(fontsize=7,frameon=False,ncol=3,loc="upper left")
ax1.set_title("16-strain unique-piRNA classification (pangenome, mm0-clean loci) — counts per strain, canonical order",fontsize=10.5,fontweight="bold")
ax1.spines[['top','right']].set_visible(False)
prop=ct.div(ct.sum(1),axis=0); bottom=np.zeros(len(CANON))
for i,k in enumerate(KL):
    ax2.bar(x,prop[k],0.74,bottom=bottom,color=COL[i],label=LAB[i]); bottom+=prop[k].values
ax2.set_xticks(x); ax2.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8)
ax2.set_ylim(0,1); ax2.set_ylabel("class composition (proportion)",fontsize=9); ax2.set_title("B  class composition per strain",fontsize=10,fontweight="bold")
ax2.spines[['top','right']].set_visible(False)
# Panel C: developmental TIMEPOINT origin of EACH class (pooled across strains) — when does each class arise?
TPO=["E16.5","P12.5","P20.5"]; TPM={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TCOL={"E16.5":"#0072B2","P12.5":"#009E73","P20.5":"#D55E00"}
d["tpL"]=d.timepoint.map(TPM); tpc=pd.crosstab(d.klass5,d.tpL).reindex(index=KL,columns=TPO).fillna(0)
tprop=tpc.div(tpc.sum(1).replace(0,np.nan),axis=0); xc=np.arange(len(KL)); bottomc=np.zeros(len(KL))
for tp in TPO:
    ax3.bar(xc,tprop[tp].values,0.6,bottom=bottomc,color=TCOL[tp],edgecolor="white",linewidth=0.4,label=tp); bottomc=bottomc+np.nan_to_num(tprop[tp].values)
ax3.set_xticks(xc); ax3.set_xticklabels([l.split(" (")[0] for l in LAB],fontsize=7); ax3.set_ylim(0,1.08)
ax3.set_ylabel("timepoint composition",fontsize=9); ax3.legend(fontsize=8,frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(0.5,1.04),title="developmental window")
ax3.set_title("C  Developmental timepoint origin per class (E16.5 prepachytene → P20.5 pachytene; pooled over 16 strains)",fontsize=10,fontweight="bold",loc="left")
for xi,k in zip(xc,KL): ax3.text(xi,1.015,f"n={int(tpc.loc[k].sum()):,}",ha="center",va="bottom",fontsize=6,color="#555")
ax3.spines[['top','right']].set_visible(False)
fig.text(0.5,0.005,"Genuinely unique = conserved-but-silent (ortholog present but truly silent elsewhere) + strain-private-locus (CLEAN = maps mm0 to a real own-genome locus). SNP-variant (orange) = allelic 1-3 mm variant EXPRESSED at the orthologous "
  "locus elsewhere -> NOT novel. low-quality (tan) = called strain-private only because the read maps mm1-3 (sequencing-error/het-SNP) to its OWN assembly and so finds no mm0 ortholog -> FILTERED OUT (60,857->20,846 clean strain-private). "
  "Wild strains (WSB/CAST/PWK/SPRET) still dominate the clean signal. Source: unique16/final_classified_clean.csv.gz (klass5); counts summed over E16.5/P12.5/P20.5.",ha="center",fontsize=6.0,color="#555")
fig.tight_layout(rect=[0,0.02,1,1])
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_unique16_class_breakdown.{e}",bbox_inches="tight")
print("wrote Fig_unique16_class_breakdown.{png,pdf,svg} + source data")
print(ct.assign(total=ct.sum(1)).to_string())
