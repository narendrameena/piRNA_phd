#!/usr/bin/env python3
"""THEME 22 step 9 — RIGOROUS test that piRNA-cluster novelty TRACKS genome/TE evolution across the 16 strains.
Per strain: divergence (1 - theme-21 graph-PAV mean = genome evolution, TE/SV-dominated), non-reference TE-insertion
burden (VCF, when available) vs the non-reference piRNA-cluster count and their piRNA-expression share. Spearman
correlation + wild-vs-classical Mann-Whitney."""
import pandas as pd, numpy as np, json, os
from scipy.stats import spearmanr, mannwhitneyu
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
div=json.load(open("/tmp/divergence.json"))
te=pd.read_csv(f"{D}/nonref_te_summary.csv")          # strain, n_nonref, n_TE_overlap
ex=pd.read_csv(f"{D}/nonref_expression_summary.csv")  # strain, nonref_expr_pct, ...
m=te.merge(ex[["strain","nonref_expr_pct"]],on="strain")
m["divergence"]=m.strain.map(div); m["wild"]=m.strain.isin(WILD)
if os.path.exists(f"{D}/te_insertion_burden.json"):
    m["te_burden"]=m.strain.map(json.load(open(f"{D}/te_insertion_burden.json")))
print("per-strain table:\n",m[["strain","divergence","n_nonref","nonref_expr_pct","wild"]+(["te_burden"] if "te_burden" in m else [])].sort_values("divergence",ascending=False).to_string(index=False))
print("\n=== SPEARMAN: does genome/TE evolution predict non-reference piRNA clusters? ===")
for x in ["divergence"]+(["te_burden"] if "te_burden" in m else []):
    for y in ["n_nonref","nonref_expr_pct"]:
        r,p=spearmanr(m[x],m[y]); print(f"  {x:11s} vs {y:16s}: rho={r:+.2f}  p={p:.3f} {'*' if p<0.05 else ''}")
if "te_burden" in m:
    r,p=spearmanr(m.divergence,m.te_burden); print(f"  divergence  vs te_burden       : rho={r:+.2f}  p={p:.3f}  (sanity: TE burden tracks divergence)")
print("\n=== WILD vs CLASSICAL (Mann-Whitney, wild > classical) ===")
for y in ["divergence","n_nonref","nonref_expr_pct"]+(["te_burden"] if "te_burden" in m else []):
    u,p=mannwhitneyu(m[m.wild][y],m[~m.wild][y],alternative="greater")
    print(f"  {y:16s}: wild median={m[m.wild][y].median():.3f} vs classical={m[~m.wild][y].median():.3f}  p={p:.3f} {'*' if p<0.05 else ''}")
m.to_csv(f"{D}/evolution_test.csv",index=False); print("\nwrote evolution_test.csv")
