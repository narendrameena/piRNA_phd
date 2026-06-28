#!/usr/bin/env python3
"""THEME 22 step 15 — aggregate the all-vs-all co-location. For each non-reference cluster (strain X), how many OTHER
strains it lifts to (share the insertion SEQUENCE) and how many also have a piCB CLUSTER there (share the CLUSTER).
Classify strain-private vs shared; sharing-degree distribution; private fraction wild vs classical; and whether sharing
is within-(wild/classical)-group (= phylogenetic, ancestral-insertion signal)."""
import collections, os, subprocess, pandas as pd
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"; D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
STR=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
share_seq=collections.defaultdict(set); share_clu=collections.defaultdict(set)
for X in STR:
    for Y in STR:
        if X==Y: continue
        f=f"{D}/colo/{X}__{Y}.lifted.bed"
        if not os.path.exists(f) or os.path.getsize(f)==0: continue
        for l in open(f): share_seq[l.rstrip("\n").split('\t')[3]].add(Y)
        cmd=f"bedtools intersect -a <(sort -k1,1 -k2,2n {f}) -b <(sort -k1,1 -k2,2n {CP}/{Y}.clusters.bed) -u 2>/dev/null | cut -f4 | sort -u"
        for cid in subprocess.run(cmd,shell=True,capture_output=True,text=True,executable='/bin/bash').stdout.split(): share_clu[cid].add(Y)
rows=[]
for X in STR:
    for l in open(f"{D}/nonref/{X}.nonref.bed"):
        cid=l.rstrip("\n").split('\t')[3]; sq=share_seq.get(cid,set()); cl=share_clu.get(cid,set()); ws=X in WILD
        rows.append(dict(cid=cid,strain=X,wild=ws,n_other_seq=len(sq),n_other_clu=len(cl),
                         n_same_group=sum(1 for y in sq if (y in WILD)==ws)))
df=pd.DataFrame(rows)
print("=== CO-LOCATION of the non-reference clusters ===")
print(f"total non-reference clusters: {len(df)}")
print(f"strain-PRIVATE  (insertion in NO other strain): {(df.n_other_seq==0).sum():4d} ({100*(df.n_other_seq==0).mean():.0f}%)")
print(f"shared insertion (>=1 other strain carries it): {(df.n_other_seq>=1).sum():4d} ({100*(df.n_other_seq>=1).mean():.0f}%)")
print(f"shared CLUSTER  (>=1 other strain clusters too): {(df.n_other_clu>=1).sum():4d} ({100*(df.n_other_clu>=1).mean():.0f}%)")
print("\nsharing degree (#other strains carrying the insertion):")
vc=df.n_other_seq.value_counts().sort_index()
for k,v in vc.items(): print(f"  +{k:2d} other strains: {v:4d}")
print("\nprivate fraction, wild vs classical:")
for w in [True,False]:
    s=df[df.wild==w]; print(f"  {'wild    ' if w else 'classical'}: {100*(s.n_other_seq==0).mean():.0f}% private  (n={len(s)})")
sh=df[df.n_other_seq>=1]
if len(sh): print(f"\nphylogenetic signal: of shared-insertion edges, {sh.n_same_group.sum()}/{sh.n_other_seq.sum()} ({100*sh.n_same_group.sum()/sh.n_other_seq.sum():.0f}%) are within the SAME (wild/classical) group")
df.to_csv(f"{D}/colocation.csv",index=False); print("\nwrote colocation.csv")
