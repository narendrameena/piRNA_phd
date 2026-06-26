#!/usr/bin/env python3
"""THEME 22 step 12 — MULTIMAPPING confounder. Non-ref clusters are 93% TE-rich; TE piRNAs multimap, so all-primary
FPM could be inflated. clusters_fpm.bed carries BOTH allFPM (col4, multimap-incl.) and uniqFPM (col5, unique reads).
Test: (a) are non-ref clusters MORE multimapping than reference (multimap fraction = 1-uniq/all)? (b) does the
'non-ref well-expressed' finding survive on UNIQUE reads? (c) does the divergence correlation survive on unique reads?
NB this also CORRECTS step 8/figure, which mistakenly summed allFPM+uniqFPM."""
import pandas as pd, numpy as np, json
from scipy.stats import mannwhitneyu, spearmanr
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"; D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}; div=json.load(open("/tmp/divergence.json"))
rows=[]; mm_nr=[]; mm_ref=[]; uq_nr=[]; uq_ref=[]; al_nr=[]; al_ref=[]
for X in S:
    f=pd.read_csv(f"{CP}/{X}.clusters_fpm.bed",sep="\t",header=None,names=["c","s","e","allF","uniqF","strand","tp"],dtype={"c":str})
    nr=pd.read_csv(f"{D}/nonref/{X}.nonref.bed",sep="\t",header=None,names=["c","s","e","id"],dtype={"c":str}); nrset=set(zip(nr.c.astype(str),nr.s,nr.e))
    f["nr"]=[(c,s,e) in nrset for c,s,e in zip(f.c.astype(str),f.s,f.e)]
    cl=f.groupby(["c","s","e","nr"],as_index=False).agg(allF=("allF","sum"),uniqF=("uniqF","sum"))
    cl["mm"]=1-cl.uniqF/cl.allF.replace(0,np.nan)
    a,b=cl[cl.nr],cl[~cl.nr]
    mm_nr+=list(a.mm.dropna()); mm_ref+=list(b.mm.dropna()); uq_nr+=list(a.uniqF); uq_ref+=list(b.uniqF); al_nr+=list(a.allF); al_ref+=list(b.allF)
    rows.append(dict(strain=X,nr_pct_all=100*a.allF.sum()/cl.allF.sum(),nr_pct_uniq=100*a.uniqF.sum()/cl.uniqF.sum(),
                     mm_nr=a.mm.median(),mm_ref=b.mm.median()))
r=pd.DataFrame(rows); r["divergence"]=r.strain.map(div); r["wild"]=r.strain.isin(WILD)
print("=== (a) are non-reference clusters MORE multimapping than reference? ===")
print(f"  multimap fraction (1-uniq/all): non-ref median={np.median(mm_nr):.3f}  vs reference={np.median(mm_ref):.3f}")
u,p=mannwhitneyu(mm_nr,mm_ref,alternative='greater'); print(f"  non-ref MORE multimapping (MW greater): p={p:.1e}  {'-> YES, TE-driven multimapping' if p<0.05 else ''}")
print("=== (b) does 'non-ref well-expressed' SURVIVE on UNIQUE reads? ===")
print(f"  allFPM : non-ref median={np.median(al_nr):.1f} vs ref={np.median(al_ref):.1f}")
print(f"  uniqFPM: non-ref median={np.median(uq_nr):.1f} vs ref={np.median(uq_ref):.1f}")
u,p=mannwhitneyu(uq_nr,uq_ref,alternative='less'); print(f"  uniqFPM non-ref < ref (MW less): p={p:.3f}  ({'non-ref STILL not lower -> genuine' if p>0.05 else 'non-ref now LOWER on unique reads -> was multimap-inflated'})")
print("=== (c) does the DIVERGENCE correlation survive on unique reads? ===")
for col in ["nr_pct_all","nr_pct_uniq"]:
    rho,pp=spearmanr(r.divergence,r[col]); print(f"  divergence vs {col:11s}: rho={rho:+.2f} p={pp:.3f} {'*' if pp<0.05 else ''}")
for col in ["nr_pct_all","nr_pct_uniq"]:
    u,pp=mannwhitneyu(r[r.wild][col],r[~r.wild][col],alternative='greater'); print(f"  wild>classical {col:11s}: wild={r[r.wild][col].median():.3f} vs {r[~r.wild][col].median():.3f} p={pp:.3f} {'*' if pp<0.05 else ''}")
r.to_csv(f"{D}/multimapping_test.csv",index=False); print("\nwrote multimapping_test.csv  (nr_pct_uniq = unique-read non-ref piRNA share)")
