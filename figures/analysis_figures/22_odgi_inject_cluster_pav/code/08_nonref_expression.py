#!/usr/bin/env python3
"""THEME 22 step 8 — do the NON-REFERENCE piRNA clusters matter for OUTPUT? For each strain: per-cluster expression =
sum(FPM+ + FPM-) over timepoints (clusters_fpm.bed); flag non-reference (coord match to nonref.bed); compute (a) the
fraction of total cluster-piRNA from non-reference clusters, (b) whether any non-ref cluster falls in the top-90%-
cumulative set (the dominant producers), (c) Mann-Whitney test that non-ref expression < reference. Creative angle:
rank of the single biggest non-reference cluster, and the cumulative coverage already reached by then."""
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"; T22=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
res=[]; allnr=[]; allref=[]
for X in S:
    fpm=pd.read_csv(f"{CP}/{X}.clusters_fpm.bed",sep="\t",header=None,names=["chrom","start","end","fp","fm","strand","tp"],dtype={"chrom":str})
    fpm["expr"]=fpm.fp+fpm.fm
    nr=pd.read_csv(f"{T22}/nonref/{X}.nonref.bed",sep="\t",header=None,names=["chrom","start","end","id"],dtype={"chrom":str})
    nrset=set(zip(nr.chrom.astype(str),nr.start,nr.end))
    fpm["nonref"]=[(c,s,e) in nrset for c,s,e in zip(fpm.chrom.astype(str),fpm.start,fpm.end)]
    cl=fpm.groupby(["chrom","start","end","nonref"],as_index=False).expr.sum()
    tot=cl.expr.sum(); nrexpr=cl[cl.nonref].expr.sum()
    cl=cl.sort_values("expr",ascending=False).reset_index(drop=True); cl["cum"]=cl.expr.cumsum()/tot
    nr_in_top90=int(cl.loc[cl.cum<=0.9,"nonref"].sum())
    best_rank=int(cl[cl.nonref].index.min())+1 if cl.nonref.any() else None
    best_cum=100*float(cl[cl.nonref].cum.min()) if cl.nonref.any() else None
    res.append(dict(strain=X,n_clusters=len(cl),n_nonref=int(cl.nonref.sum()),nonref_expr_pct=round(100*nrexpr/tot,3),
                    n_top90=int((cl.cum<=0.9).sum()),nr_in_top90=nr_in_top90,best_nonref_rank=best_rank,cum_at_best_nonref=round(best_cum,1) if best_cum else None))
    allnr+=list(cl[cl.nonref].expr); allref+=list(cl[~cl.nonref].expr)
r=pd.DataFrame(res); print(r.to_string(index=False))
print(f"\n=== OVERALL ===")
print(f"non-reference clusters carry {100*sum(allnr)/(sum(allnr)+sum(allref)):.3f}% of total cluster piRNA  (per-strain {r.nonref_expr_pct.min():.2f}-{r.nonref_expr_pct.max():.2f}%)")
print(f"non-reference clusters inside the top-90%-cumulative (dominant) set: {r.nr_in_top90.sum()} of {r.n_nonref.sum()} total non-ref")
print(f"the single biggest non-ref cluster per strain ranks {r.best_nonref_rank.min()}-{r.best_nonref_rank.max()} of ~{int(r.n_clusters.mean())}; by then {r.cum_at_best_nonref.min():.0f}-{r.cum_at_best_nonref.max():.0f}% of piRNA is already covered by other (reference) clusters")
u,p=mannwhitneyu(allnr,allref,alternative='less')
print(f"TEST non-ref expression < reference (Mann-Whitney U): median non-ref={np.median(allnr):.2f} vs reference={np.median(allref):.2f} FPM, p={p:.2e}")
r.to_csv(f"{T22}/nonref_expression_summary.csv",index=False)
