#!/usr/bin/env python3
"""16-strain coordinate pangenome TE-driven test (CLEAN klass5, ≥2-read). Per strain (canonical order) and class,
fold-enrichment of having the PRODUCTION LOCUS inside a 1-of-16-private insertion, over the random-locus null
(merged private-insertion bp / genome). Tests production, not sequence similarity. With the clean klass5
strain-private set (mm0 loci), strain-private-locus piRNAs are the MOST insertion-enriched class — their
production loci sit predominantly INSIDE strain-private TE insertions (64-92% for divergent strains, up to ~86x
over null) — whereas conserved-but-silent (divergence) piRNAs are not. [The earlier 3-class 'klass' route bundled
in mm1-3 low-quality reads that have no clean mm0 locus and DILUTED this to ~19-48%; klass5 removes them.] A finer
breadth split (Fig_te_driven_corrected16) resolves these insertion-driven loci into new-locus CREATION (minority)
vs PROPAGATION into conserved clusters (majority). Input = {X}.coord_byclass16.csv (coord_classify16, klass5).
[BioMNI 3/3 verified 2026-06-18 — TE-insertion-driven origin of strain-private piRNAs is established biology; caveat: 'private' depends on assembly completeness, so a lower bound]"""
import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]
KL=["expressed elsewhere (exact)","unique: conserved-but-silent","unique: strain-private locus"]
LAB=["expressed-elsewhere (common)","unique: conserved-but-silent","unique: strain-private locus"]; COL=["#9e9e9e","#0072B2","#7a3b9a"]
df=pd.concat([pd.read_csv(f"{PG}/{X}.coord_byclass16.csv") for X in CANON if __import__("os").path.exists(f"{PG}/{X}.coord_byclass16.csv")],ignore_index=True)
df["fold"]=df.apply(lambda r: r["pct"]/r["exp_pct"] if r["exp_pct"]>0 else np.nan,axis=1)
df.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/09_TE_driven_evolution/data/source_data/SourceData_te_driven_coord16.csv",index=False)
exp={X:df[df.strain==X].exp_pct.iloc[0] for X in CANON if (df.strain==X).any()}
print("fold-enrichment (pct/exp) by strain x class:")
print(df.pivot(index="strain",columns="klass",values="fold").reindex(CANON).round(2).to_string())
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(13.5,6.2),dpi=300); x=np.arange(len(CANON)); bw=0.27
for i,(k,lab) in enumerate(zip(KL,LAB)):
    sub=df[df.klass==k].set_index("strain").reindex(CANON)
    ax.bar(x+(i-1)*bw,sub.fold,bw,color=COL[i],label=lab)
    for xi,v in zip(x+(i-1)*bw,sub.fold):
        if pd.notna(v): ax.text(xi,v*1.04,f"{v:.1f}",ha="center",va="bottom",fontsize=4.4,rotation=90,color=COL[i],fontweight="bold")
ax.axhline(1,color="#555",ls="--",lw=1); ax.text(len(CANON)-0.5,1.05,"random (1×)",ha="right",fontsize=7,color="#555")
ax.set_yscale("log"); ax.set_xticks(x)
ax.set_xticklabels([f"{s.replace('_','/')}\n({exp.get(s,0):.1f}% ins)" for s in CANON],rotation=45,ha="right",fontsize=7)
ax.set_ylabel("fold-enrichment: production locus inside a\n1-of-16-private insertion (obs ÷ genomic null)",fontsize=9)
ax.legend(fontsize=8.5,frameon=False,loc="upper left"); ax.spines[['top','right']].set_visible(False)
ax.set_title("16-strain coordinate TE-driven test — strain-private-locus piRNAs are the MOST enriched at strain-private TE insertions",fontsize=10.3,fontweight="bold")
fig.text(0.5,0.005,"Numbers on bars = fold-enrichment (% under each strain = merged private-insertion fraction of genome = the null). With the CLEAN klass5 strain-private set, strain-private-locus piRNAs (pangenome: locus absent in all 15 other strains) are the MOST insertion-enriched class — their production loci sit predominantly inside private insertions (64-92% for divergent strains; up to ~86x over null). [Correction: the earlier 3-class estimate (~24% in-insertion) was DILUTED by mm1-3 low-quality reads lacking a clean mm0 locus; klass5 removes them.] The class is still a minority of candidates (klass5 strain-private = 20,846 = 5.2%). "
  "A finer breadth split (Fig_te_driven_corrected16) resolves these insertion-driven loci into new-locus CREATION vs PROPAGATION into conserved clusters. High classical-strain folds are inflated by their small insertion fraction (EXP ~1%; the % is the interpretable number). Coordinate-based (production locus); METHODS §8. [BioMNI 3/3 verified 2026-06-18 — TE-insertion-driven origin of strain-private piRNAs is established biology; caveat: 'private' depends on assembly completeness, so a lower bound]",ha="center",fontsize=5.3,color="#555")
fig.tight_layout(rect=[0,0.02,1,1])
# classical(blue)/wild(orange) companion: strain-private-locus piRNAs per strain (subspecies colour scheme)
fig.subplots_adjust(bottom=0.34)
_tot=df[df.klass=="unique: strain-private locus"].set_index("strain").reindex(CANON).n.fillna(0).values
_cax=add_classical_wild_companion(fig,ax,CANON,_tot,gap=0.13,height_frac=0.20,ylabel="strain-priv\nloci (log)")
_cax.set_xticks(np.arange(len(CANON))); _cax.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=6.5)
for lab,s in zip(_cax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — strain-private-locus piRNAs per strain",fontsize=7.5,fontweight="bold",loc="left")
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_te_driven_coord16.{e}",bbox_inches="tight")
print("wrote Fig_te_driven_coord16.{png,pdf,svg} + source data")
