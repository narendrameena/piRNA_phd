#!/usr/bin/env python3
"""THEME 22 step 16 — TE-AGE test. Hypothesis: CONSERVED/reference clusters sit on OLD TEs, NON-REFERENCE clusters on
YOUNG TEs (recent insertions = the arms race). TE age proxy = RepeatMasker % divergence from consensus (col 2; LOWER =
younger). For each strain, intersect non-ref vs reference (=lifted) clusters with the RM annotation; compare the
divergence of the overlapping TEs (Mann-Whitney, non-ref < reference). Also WITHIN the major families (L1, ERVK) to
show it is insertion-age, not just family composition."""
import pandas as pd, glob, subprocess, numpy as np, collections, os, json
from scipy.stats import mannwhitneyu
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"; D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
nr=collections.defaultdict(list); rf=collections.defaultdict(list)
for X in S:
    rmf=glob.glob(f"{B}/resources/repeatMasker/{X}_*.out")
    if not rmf: continue
    rmb=f"/tmp/rmage_{X}.bed"
    if not os.path.exists(rmb) or os.path.getsize(rmb)==0:
        subprocess.run(f"awk 'NR>3{{print $5\"\\t\"$6\"\\t\"$7\"\\t\"$2\"\\t\"$11}}' {rmf[0]} | sort -k1,1 -k2,2n > {rmb}",shell=True,executable='/bin/bash')
    nrset=set((c,s,e) for c,s,e in (l.rstrip().split('\t')[:3] for l in open(f"{D}/nonref/{X}.nonref.bed")))
    allc=pd.read_csv(f"{CP}/{X}.clusters.bed",sep="\t",header=None,usecols=[0,1,2],names=["c","s","e"],dtype={"c":str})
    allc["isnr"]=[(c,str(s),str(e)) in nrset for c,s,e in zip(allc.c.astype(str),allc.s,allc.e)]
    for grp,sub,store in [("nr",allc[allc.isnr],nr),("rf",allc[~allc.isnr],rf)]:
        b=f"/tmp/{grp}age_{X}.bed"
        sub.assign(c=X+"#1#chr"+sub.c.astype(str))[["c","s","e"]].sort_values(["c","s"]).to_csv(b,sep="\t",header=False,index=False)
        out=subprocess.run(f"bedtools intersect -a {b} -b {rmb} -wa -wb 2>/dev/null | cut -f7,8",shell=True,capture_output=True,text=True,executable='/bin/bash').stdout
        for ln in out.splitlines():
            d,fam=ln.split('\t'); store["ALL"].append(float(d)); store[fam].append(float(d))
    print(f"  {X}: non-ref TEs {len(nr['ALL'])}  reference TEs {len(rf['ALL'])}",flush=True)
print("\n=== TE AGE (RepeatMasker % divergence from consensus; LOWER = YOUNGER) ===")
print(f"ALL TEs: non-reference median {np.median(nr['ALL']):.1f}%  vs  reference/conserved median {np.median(rf['ALL']):.1f}%")
u,p=mannwhitneyu(nr['ALL'],rf['ALL'],alternative='less'); print(f"  non-reference TEs YOUNGER than reference (MW less): p={p:.2e}  {'-> CONFIRMED' if p<0.05 else ''}")
print("WITHIN major families (controls for family composition):")
famrows=[]
for fam in ["LINE/L1","LTR/ERVK","LTR/ERVL-MaLR","LTR/ERV1","SINE/B2","SINE/Alu","SINE/B4","LTR/ERVL"]:
    if len(nr.get(fam,[]))>=20 and len(rf.get(fam,[]))>=20:
        u,p=mannwhitneyu(nr[fam],rf[fam],alternative='less')
        print(f"  {fam:16s}: non-ref {np.median(nr[fam]):5.1f}% vs ref {np.median(rf[fam]):5.1f}%  p={p:.1e} {'*' if p<0.05 else ''}  (n {len(nr[fam])}/{len(rf[fam])})")
        famrows.append((fam,round(np.median(nr[fam]),2),round(np.median(rf[fam]),2),len(nr[fam]),len(rf[fam]),p))
pd.DataFrame(famrows,columns=["family","nonref_div","ref_div","n_nonref","n_ref","p"]).to_csv(f"{D}/te_age_byfamily.csv",index=False)
pd.DataFrame({"set":["non_reference","reference"],"median_div_pct":[np.median(nr['ALL']),np.median(rf['ALL'])],"n_TEs":[len(nr['ALL']),len(rf['ALL'])]}).to_csv(f"{D}/te_age_test.csv",index=False)
import json; json.dump({"nonref_div":nr['ALL'][:200000],"ref_div":rf['ALL'][:200000]},open(f"{D}/te_age_dists.json","w"))
print("\nwrote te_age_test.csv + te_age_dists.json")
