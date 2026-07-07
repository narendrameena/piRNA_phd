#!/usr/bin/env python3
"""THEME 22 step 8 — do the NON-REFERENCE piRNA clusters matter for OUTPUT? For each strain: per-MERGED-cluster
expression = all-primary FPM summed over timepoints+strands (via _nonref_util.merged_cluster_expr, which aggregates
the unmerged clusters_fpm.bed to the merged clusters.bed intervals — see that module for the exact-match bug this
fixes); flag non-reference (nonref.bed, a merged subset). Compute (a) the fraction of total cluster-piRNA from
non-reference clusters, (b) whether any non-ref cluster falls in the top-90%-cumulative (dominant) set, (c) a
two-sided Mann-Whitney of non-ref vs reference per-cluster expression — the non-ref clusters are FEW (~0.7% of total
piRNA output) but individually WELL-expressed (median > reference). Creative angle: rank of the single biggest
non-reference cluster, and the cumulative coverage already reached by then."""
import pandas as pd, numpy as np, sys, os as _os
from scipy.stats import mannwhitneyu
sys.path.insert(0,_os.path.dirname(_os.path.abspath(__file__))); from _nonref_util import merged_cluster_expr
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"; T22=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
res=[]; allnr=[]; allref=[]
for X in S:
    cl=merged_cluster_expr(CP,T22,X)   # per MERGED cluster (fixes the exact (chrom,start,end) match that dropped ~12% of non-ref clusters: nonref.bed is MERGED, clusters_fpm.bed UNMERGED)
    cl["expr"]=cl.allF   # all-primary FPM (col4); NOT allF+uniqF (that was the double-count)
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
u,p=mannwhitneyu(allnr,allref,alternative='two-sided')
_dir="higher" if np.median(allnr)>np.median(allref) else "lower"
print(f"TEST non-ref vs reference per-cluster expression (Mann-Whitney U, two-sided): median non-ref={np.median(allnr):.2f} vs reference={np.median(allref):.2f} FPM (non-ref individually {_dir}), p={p:.2e}")
r.to_csv(f"{T22}/nonref_expression_summary.csv",index=False)
