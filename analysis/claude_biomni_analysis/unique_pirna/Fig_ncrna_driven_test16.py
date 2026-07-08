#!/usr/bin/env python3
"""16-strain ncRNA-driven test — the 16-strain version of the pilot Fig_ncrna_driven_test. For each strain
(canonical order) and each UNIQUE class (conserved-but-silent, strain-private-locus), fold-enrichment of
uniquely-mapping piRNA loci inside CLEAN lncRNA genes (lncRNA MINUS protein_coding — confounding-free), over the
random-locus expectation. Loci = unique16/loci/{X}.cand_loci.ens.bed (novel candidates mapped to own genome,
Ensembl coords); class via join to final_classified on strain|tp|seq. (expressed-elsewhere has no novel loci, so
the two unique classes are shown.)"""
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np, bisect
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import WILD, add_classical_wild_companion
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"; ANN=f"{ROOT}/resources/annotation"
CANON=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
MAIN={str(i) for i in range(1,20)}|{"X"}
CLS=["unique: conserved-but-silent","unique: strain-private locus"]; LAB=["conserved-but-silent","strain-private locus"]; COL=["#0072B2","#7a3b9a"]
def stripc(c): return c.split("#")[-1].replace("chr","")
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")   # CANONICAL 5-class (klass5): CBS excl SNP-variants; strain-private = mm0-clean
kl={f"{s}|{t}|{q}":k for s,t,q,k in zip(d.strain,d.timepoint,d.sequence,d.klass5)}
def clean_space(X):
    lnc=defaultdict(list); pcg=defaultdict(list)
    for ln in open(f"{ANN}/{X}_v3.5.gff3"):
        if ln[0]=="#": continue
        f=ln.split("\t")
        if len(f)<9 or f[2] not in ("gene","ncRNA_gene"): continue
        c=stripc(f[0])
        if c not in MAIN: continue
        if "biotype=lncRNA" in f[8]: lnc[c].append((int(f[3]),int(f[4])))
        elif "biotype=protein_coding" in f[8]: pcg[c].append((int(f[3]),int(f[4])))
    def mrg(dd):
        o={}
        for c,iv in dd.items():
            iv.sort(); m=[]
            for s,e in iv:
                if m and s<=m[-1][1]: m[-1]=[m[-1][0],max(m[-1][1],e)]
                else: m.append([s,e])
            o[c]=m
        return o
    lm=mrg(lnc); pm=mrg(pcg); pss={c:[s for s,e in pm[c]] for c in pm}; pee={c:[e for s,e in pm[c]] for c in pm}
    def in_pc(c,s,e):
        if c not in pss: return False
        i=bisect.bisect_right(pss[c],e)-1; return i>=0 and pee[c][i]>s
    cb=0; ss={}; ee={}
    for c,m in lm.items():
        ss[c]=[s for s,e in m]; ee[c]=[e for s,e in m]
        for s,e in m:
            ov=0
            if c in pss:
                for k in range(max(0,bisect.bisect_left(pss[c],s)-1),bisect.bisect_left(pss[c],e)+1):
                    if k<len(pm[c]): a,b=pm[c][k]; ov+=max(0,min(b,e)-max(a,s))
            cb+=(e-s)-ov
    def inlnc(c,s,e):
        if c not in ss: return False
        i=bisect.bisect_right(ss[c],e)-1; return i>=0 and ee[c][i]>s and not in_pc(c,s,e)
    return cb,inlnc
rows=[]
for X in CANON:
    gsize=sum(int(l) for c,l in (ln.split() for ln in open(f"{ROOT}/results/indexs/{X}/chrNameLength.txt")) if stripc(c) in MAIN)
    cb,inlnc=clean_space(X); lncfrac=cb/gsize
    loci=defaultdict(list)
    for ln in open(f"{U}/unique16/loci/{X}.cand_loci.ens.bed"):
        p=ln.split("\t")
        if p[0] in MAIN: loci[p[3]].append((p[0],int(p[1]),int(p[2])))
    for cls,lab in zip(CLS,LAB):
        ids=[n for n,L in loci.items() if len(L)==1 and kl.get(n)==cls]
        n_in=sum(inlnc(*loci[n][0]) for n in ids); fold=(n_in/len(ids))/lncfrac if ids and lncfrac else np.nan
        rows.append(dict(strain=X,cls=lab,fold=round(fold,3) if ids else np.nan,n_uniq=len(ids),n_in_lncRNA=n_in,clean_lncRNA_pct=round(lncfrac*100,2)))
    print(f"{X}: clean-lncRNA={lncfrac*100:.1f}% | "+" | ".join(f"{lab}:{r['fold']}x(n={r['n_uniq']:,})" for lab,r in zip(LAB,rows[-2:])))
df=pd.DataFrame(rows); df.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/data/source_data/SourceData_ncrna_driven_test16.csv",index=False)
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(13,6),dpi=300); x=np.arange(len(CANON)); bw=0.38
for i,(cls,lab) in enumerate(zip(CLS,LAB)):
    sub=df[df.cls==lab].set_index("strain").reindex(CANON)
    ax.bar(x+(i-0.5)*bw,sub.fold,bw,color=COL[i],label=f"unique: {lab}")
    for xi,v,nn in zip(x+(i-0.5)*bw,sub.fold,sub.n_uniq):
        if pd.notna(v): ax.text(xi,(v if v>0 else 0)+0.06,f"n={int(nn)}",ha="center",va="bottom",fontsize=4.2,rotation=90,color=COL[i])
ax.axhline(1,color="#555",ls="--",lw=1); ax.text(len(CANON)-1,1.04,"random expectation (1×)",ha="right",fontsize=7,color="#555")
ax.set_xticks(x); ax.set_xticklabels([])   # strain labels carried by the classical/wild companion below
ax.set_ylabel("fold-enrichment at CLEAN lncRNA genes\n(observed ÷ genomic expectation)",fontsize=9); ax.legend(fontsize=8.5,frameon=False)
ax.set_title("16-strain ncRNA-driven test — unique piRNAs vs clean lncRNA genes (protein-coding-excluded, confounding-free)",fontsize=10.5,fontweight="bold")
ax.spines[['top','right']].set_visible(False)
fig.text(0.5,0.005,"Per strain (canonical order), fold-enrichment of uniquely-mapping NOVEL-class piRNA loci inside clean lncRNA genes (lncRNA minus protein_coding). Loci = cand_loci16 (own genome, Ensembl); class via join to final_classified. "
  "n above each bar = uniquely-mapping loci: CLASSICAL strains have very low n (low power — noisy folds); the robust signal is the high-n wild-derived strains (WSB/CAST/PWK/SPRET), modestly enriched (~1.4–2×), consistent with the pilot.",ha="center",fontsize=6.0,color="#555")
fig.tight_layout(rect=[0,0.02,1,1])
# classical(blue)/wild(orange) companion: uniquely-mapping unique-piRNA loci per strain (subspecies colour scheme)
fig.subplots_adjust(bottom=0.34)
_tot=df.groupby("strain").n_uniq.sum().reindex(CANON).fillna(0).values
_cax=add_classical_wild_companion(fig,ax,CANON,_tot,gap=0.13,height_frac=0.20,ylabel="uniq-map\nloci (log)")
_cax.set_xticks(np.arange(len(CANON))); _cax.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=6.5)
for lab,s in zip(_cax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — uniquely-mapping unique-piRNA loci per strain",fontsize=7.5,fontweight="bold",loc="left")
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_ncrna_driven_test16.{e}",bbox_inches="tight")
print("wrote Fig_ncrna_driven_test16.{png,pdf,svg} + source data")
