#!/usr/bin/env python3
"""THEME 22 step 11 — CONFOUNDING checks for the 'non-ref piRNA tracks divergence' claim.
(1) total-output scaling: partial Spearman controlling for total cluster count + total FPM.
(2) assembly fragmentation: does per-strain graph-fragment count drive it?
(3) mappability/lift-artifact: does the non-ref RATE (non-ref/native) just track divergence (= possible lift failure)?
(4) phylogenetic non-independence: does the correlation hold WITHIN the classical strains only (drops the wild block)?"""
import pandas as pd, numpy as np, json, collections
from scipy.stats import spearmanr, rankdata, pearsonr
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"; D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
ev=pd.read_csv(f"{D}/evolution_test.csv")
# per-strain total clusters + total FPM
tot=[]
for X in S:
    fpm=pd.read_csv(f"{CP}/{X}.clusters_fpm.bed",sep="\t",header=None,names=["c","s","e","fp","fm","st","tp"]); tot.append((X,len(fpm),(fpm.fp+fpm.fm).sum()))
m=ev.merge(pd.DataFrame(tot,columns=["strain","n_clusters_total","total_fpm"]),on="strain")
# fragment count per strain (graph fragmentation)
frag=collections.Counter()
for l in open(f"{D}/graph_path_lengths.txt"):
    p=l.split("\t")[0].split("#")
    if len(p)>=4: frag[p[0]]+=1
m["n_fragments"]=m.strain.map(frag)
# non-ref RATE (mappability proxy): non-ref / native
def native(X): return sum(1 for _ in open(f"{CP}/{X}.clusters.bed"))
m["n_native"]=m.strain.map(native); m["nonref_rate"]=100*m.n_nonref/m.n_native
def pspear(df,x,y,z):
    rx,ry,rz=rankdata(df[x]),rankdata(df[y]),rankdata(df[z])
    res=lambda a:a-np.polyval(np.polyfit(rz,a,1),rz)
    return pearsonr(res(rx),res(ry))
print("=== (0) headline (uncontrolled) ===")
for y in ["n_nonref","nonref_expr_pct"]:
    r,p=spearmanr(m.divergence,m[y]); print(f"  divergence vs {y}: rho={r:+.2f} p={p:.3f}")
print("=== (1) are non-ref metrics confounded by TOTAL output/clusters? ===")
for c in ["n_clusters_total","total_fpm","n_fragments"]:
    rd,pd_=spearmanr(m.divergence,m[c]); re,pe=spearmanr(m.nonref_expr_pct,m[c])
    print(f"  {c:18s}: vs divergence rho={rd:+.2f}(p={pd_:.2f}) | vs nonref_expr% rho={re:+.2f}(p={pe:.2f})")
print("  PARTIAL divergence vs nonref_expr% controlling for:")
for z in ["n_clusters_total","total_fpm","n_fragments"]:
    r,p=pspear(m,"divergence","nonref_expr_pct",z); print(f"     | {z:18s}: partial r={r:+.2f} p={p:.3f}")
print("=== (3) mappability/lift-artifact: non-ref RATE vs divergence (strong => suspicious) ===")
r,p=spearmanr(m.divergence,m.nonref_rate); print(f"  nonref_RATE(%) vs divergence: rho={r:+.2f} p={p:.3f} | rate range {m.nonref_rate.min():.2f}-{m.nonref_rate.max():.2f}% (all tiny)")
print("=== (4) phylogenetic non-independence: WITHIN classical strains only (n=12) ===")
cl=m[~m.strain.isin(WILD)]
for y in ["n_nonref","nonref_expr_pct"]:
    r,p=spearmanr(cl.divergence,cl[y]); print(f"  classical-only divergence vs {y}: rho={r:+.2f} p={p:.3f} (n={len(cl)})")
m.to_csv(f"{D}/confounding.csv",index=False)
print("\nwrote confounding.csv")
