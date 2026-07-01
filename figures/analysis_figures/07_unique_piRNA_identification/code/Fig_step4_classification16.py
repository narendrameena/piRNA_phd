#!/usr/bin/env python3
"""16-strain unique-piRNA classification — CANONICAL pangenome-syntenic 5-class (klass5), the 'how genuinely
unique are they' view. REORGANISED to a compact 2×2 (was a 4-row stack): A counts per strain (log) + D
developmental timepoint origin BY CLASS on top; B class composition per strain spanning the bottom. The old
per-strain timepoint-origin panel was REMOVED — it is covered by Fig_unique_pirna_timepoint (both genuinely-
unique subcategories, per strain × timepoint, side-by-side). expressed-elsewhere (exact) / SNP-variant
(1-3 mm allele expressed at the orthologous locus elsewhere) = NOT unique; conserved-but-silent / strain-private
locus = GENUINELY unique by expression. Input = unique16/final_classified_clean_2read.csv.gz (klass5):
mm0-clean strain-private; orthologous-locus presence via the pangenome (halLiftover); SNP-variant = the other
strain's allele at the SAME locus is 1-3 mm from the piRNA AND in its >=2-read pool."""
import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion, classical_wild_legend_handles
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")   # ADOPTED ≥2-read absence (5-class klass5)
CANON=[s for s in STRAIN_ORDER if s in set(d.strain)]
WPOS=[i for i,s in enumerate(CANON) if s in WILD]   # wild-derived block (contiguous at the end of canonical order)
def _wildmark(ax,span=True):
    for lab,s in zip(ax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
    if span and WPOS: ax.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
KL=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
LAB=["expressed-elsewhere (not unique)","SNP-variant (allelic — not unique)","low-quality (no mm0 own locus)","unique: conserved-but-silent","unique: strain-private locus (clean)"]
COL=["#9e9e9e","#E69F00","#cdb892","#0072B2","#7a3b9a"]
TPC={"16.5dpc":"#0072B2","12.5dpp":"#009E73","20.5dpp":"#D55E00"}; TPL={"16.5dpc":"E16.5 (prepachytene)","12.5dpp":"P12.5","20.5dpp":"P20.5 (pachytene)"}
ct=pd.crosstab(d.strain,d.klass5).reindex(CANON)[KL]; ct.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data/SourceData_step4_classification16.csv")
x=np.arange(len(CANON)); bw=0.16
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(15,10.8),dpi=300)
gs=fig.add_gridspec(2,2,height_ratios=[1.0,0.94],width_ratios=[1.95,1.0],hspace=0.46,wspace=0.17,left=0.065,right=0.985,top=0.875,bottom=0.225)
ax1=fig.add_subplot(gs[0,0]); ax4=fig.add_subplot(gs[0,1]); ax2=fig.add_subplot(gs[1,:])
# ---- A: 5-class counts per strain (grouped, log) ----
for i,(k,lab) in enumerate(zip(KL,LAB)):
    ax1.bar(x+(i-2)*bw,ct[k],bw,color=COL[i],label=lab)
    for xi,v in zip(x+(i-2)*bw,ct[k]): ax1.text(xi,v*1.08,f"{int(v):,}",ha="center",va="bottom",fontsize=5.5,rotation=90,color=COL[i],fontweight="bold")
ax1.set_yscale("log"); ax1.set_ylim(20,ct.values.max()*3.6); ax1.set_xticks(x); ax1.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.3)
ax1.set_ylabel("candidates (count, log)",fontsize=9)
ax1.set_title("A  Unique-piRNA 5-class counts per strain (log)",fontsize=10.2,fontweight="bold",loc="left",pad=6)
_wildmark(ax1); ax1.text(np.mean(WPOS),ct.values.max()*2.1,"wild-derived",ha="center",va="top",fontsize=7.2,fontweight="bold",color="#C0392B")
ax1.spines[['top','right']].set_visible(False)
# ---- D: developmental timepoint origin BY CLASS (pooled over 16 strains) ----
KLAB={"expressed elsewhere (exact)":"expressed-\nelsewhere","SNP-variant (1-3mm)":"SNP-\nvariant","low-quality: no mm0 own-genome locus":"low-\nquality","unique: conserved-but-silent":"unique:\nCBS","unique: strain-private locus":"unique:\nstrain-private"}
tpx=pd.crosstab(d.klass5,d.timepoint).reindex(index=KL)
for t in TPC:
    if t not in tpx.columns: tpx[t]=0.0
tpxp=tpx[list(TPC)].div(tpx[list(TPC)].sum(1).replace(0,1),axis=0); xb4=np.arange(len(KL)); bot4=np.zeros(len(KL))
for t in TPC:
    ax4.bar(xb4,tpxp[t].values,0.62,bottom=bot4,color=TPC[t],label=TPL[t],edgecolor="white",linewidth=0.4); bot4+=tpxp[t].values
ax4.set_xticks(xb4); ax4.set_xticklabels([KLAB[k] for k in KL],fontsize=7.4); ax4.set_ylim(0,1)
ax4.set_ylabel("timepoint origin (proportion)",fontsize=9)
ax4.legend(fontsize=7.2,frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(0.5,-0.20),columnspacing=1.6,handlelength=1.4)
ax4.set_title("D  Timepoint origin BY CLASS (pooled; text = n)",fontsize=10.2,fontweight="bold",loc="left",pad=14)
for xi,k in zip(xb4,KL): ax4.text(xi,1.012,f"{int(tpx.loc[k,list(TPC)].sum()):,}",ha="center",va="bottom",fontsize=6,color="#777")
ax4.spines[['top','right']].set_visible(False)
tpx.assign(total=tpx[list(TPC)].sum(1)).to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data/SourceData_step4_classification16_tp_by_class.csv")
# ---- B: class composition per strain (stacked) ----
prop=ct.div(ct.sum(1),axis=0); bottom=np.zeros(len(CANON))
for i,(k,lab) in enumerate(zip(KL,LAB)):
    ax2.bar(x,prop[k],0.74,bottom=bottom,color=COL[i]); bottom+=prop[k].values
gu=ct[["unique: conserved-but-silent","unique: strain-private locus"]].sum(1)/ct.sum(1)*100
ax2.set_xticks(x); ax2.set_xticklabels([]); ax2.set_ylim(0,1)
ax2.set_ylabel("class composition (proportion)",fontsize=9)
ax2.set_title("B  Class composition per strain (text = % genuinely unique by expression)",fontsize=10.2,fontweight="bold",loc="left",pad=8)
for xi,s in zip(x,CANON): ax2.text(xi,1.01,f"{gu[s]:.0f}%",ha="center",va="bottom",fontsize=5.8,color="#333")
_wildmark(ax2); ax2.spines[['top','right']].set_visible(False)
# classical(blue)/wild(orange) total-count companion below the composition (subspecies colour scheme, per request)
_cax=add_classical_wild_companion(fig, ax2, CANON, ct.sum(1).reindex(CANON).values, gap=0.02, height_frac=0.16, ylabel="total\ncands (log)")
_cax.set_xticks(x); _cax.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.5)
for lab,s in zip(_cax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — total candidates per strain",fontsize=7.5,fontweight="bold",loc="left")
# shared 5-class legend across the top (covers A and B)
fig.legend(handles=[Patch(facecolor=COL[i],label=LAB[i]) for i in range(len(KL))],ncol=5,loc="upper center",
           bbox_to_anchor=(0.5,0.945),frameon=False,fontsize=8.4,columnspacing=1.4,handlelength=1.4)
fig.suptitle("16-strain pangenome-syntenic unique-piRNA classification (klass5, mm0-clean strain-private) — counts, composition, and per-class developmental origin",fontsize=11,fontweight="bold",y=0.99)
fig.text(0.5,0.012,"Pangenome-syntenic expression test: a candidate is NOT unique if another strain expresses the identical (exact) or a ≤3-mm SNP-variant allele at the SAME orthologous locus; GENUINELY unique = orthologous locus present elsewhere but truly silent there (conserved-but-silent) or no orthologous locus at all (strain-private locus). "
  "Orthology via the pangenome (halLiftover); coordinate-based, not sequence-containment. Strain-private = mm0-clean own-genome loci (low-quality tan = mm1-3 error/het-SNP, filtered). The per-strain timepoint-origin panel was moved to Fig_unique_pirna_timepoint. Source: unique16/final_classified_clean_2read.csv.gz (klass5).",ha="center",fontsize=6.0,color="#555")
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_step4_classification16.{e}",bbox_inches="tight")
print("wrote Fig_step4_classification16 (reorganised 2×2: A counts, D class×timepoint, B composition; dropped per-strain timepoint panel)")
print(ct.assign(total=ct.sum(1)).to_string())
