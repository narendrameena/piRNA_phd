#!/usr/bin/env python3
"""THEME 22 step 17 — characterise the 'present-in-most-but-not-GRCm39' subset (non-ref clusters shared with >=14 other
strains = present in 15-16/16 strains yet absent from the C57BL/6J reference). (A) collapse the strain-entries to
DISTINCT loci via the co-location graph; (B) what they are — dominant TE family, TE age (RM %div), unique-read
expression; vs the rest of the non-reference set. The 'why absent from GRCm39' (deletion vs gap vs technical) is step 18."""
import pandas as pd, numpy as np, collections, os, subprocess
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"; D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
STR=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
TH=14
co=pd.read_csv(f"{D}/colocation.csv")
coords={}
for X in STR:
    for l in open(f"{D}/nonref/{X}.nonref.bed"):
        f=l.rstrip("\n").split('\t'); coords[f[3]]=(f[0],int(f[1]),int(f[2]))
nu=co[co.n_other_seq>=TH].copy()
nu["chrom"]=nu.cid.map(lambda c:coords[c][0]); nu["start"]=nu.cid.map(lambda c:coords[c][1]); nu["end"]=nu.cid.map(lambda c:coords[c][2])
nu_ids=set(nu.cid)
print(f"=== near-universal subset (shared with >={TH} other strains = present in {TH+1}-16 strains, absent GRCm39): {len(nu)} strain-entries ===")
# --- collapse via co-location graph ---
parent={c:c for c in nu_ids}
def find(x):
    r=x
    while parent[r]!=r: r=parent[r]
    while parent[x]!=r: parent[x],x=r,parent[x]
    return r
def union(a,b): parent[find(a)]=find(b)
nu_by=collections.defaultdict(list)
for _,r in nu.iterrows(): nu_by[r.strain].append((r.chrom,r.start,r.end,r.cid))
for X in STR:
    for Y in STR:
        if X==Y: continue
        f=f"{D}/colo/{X}__{Y}.lifted.bed"
        if not os.path.exists(f): continue
        for l in open(f):
            c=l.rstrip("\n").split('\t')
            if c[3] not in nu_ids: continue
            lc,ls,le=c[0],int(c[1]),int(c[2])
            for (yc,ys,ye,ycid) in nu_by[Y]:
                if yc==lc and not(le<ys or ls>ye): union(c[3],ycid); break
comp=collections.defaultdict(list)
for c in nu_ids: comp[find(c)].append(c)
print(f"collapse -> {len(comp)} DISTINCT loci")
print(f"distinct-strains-per-locus: {sorted([len(set(c.split('|')[0] for c in v)) for v in comp.values()],reverse=True)}")
# --- TE + expression per entry (batched per strain) ---
te_fam={}; te_div={}; uq={}
for X in STR:
    sub=nu[nu.strain==X]
    if not len(sub): continue
    bed="/tmp/nu_%s.bed"%X
    sub.assign(c=X+"#1#chr"+sub.chrom.astype(str))[["c","start","end","cid"]].sort_values(["c","start"]).to_csv(bed,sep="\t",header=False,index=False)
    rmb=f"/tmp/rmage_{X}.bed"
    if os.path.exists(rmb):
        out=subprocess.run(f"bedtools intersect -a {bed} -b {rmb} -wa -wb 2>/dev/null",shell=True,capture_output=True,text=True,executable='/bin/bash').stdout
        byid=collections.defaultdict(list)
        for ln in out.splitlines():
            f=ln.split('\t'); byid[f[3]].append((float(f[7]),f[8],int(f[6])-int(f[5])))  # div, fam, len
        for cid,hits in byid.items():
            hits.sort(key=lambda h:-h[2]); te_fam[cid]=hits[0][1]; te_div[cid]=min(h[0] for h in hits)
    fpm=pd.read_csv(f"{CP}/{X}.clusters_fpm.bed",sep="\t",header=None,names=["c","s","e","allF","uniqF","st","tp"],dtype={"c":str})
    g=fpm.groupby(["c","s","e"],as_index=False).uniqF.sum()
    gmap={(str(c),s,e):u for c,s,e,u in zip(g.c,g.s,g.e,g.uniqF)}
    for _,r in sub.iterrows(): uq[r.cid]=gmap.get((str(r.chrom),r.start,r.end),0.0)
nu["te_family"]=nu.cid.map(te_fam); nu["te_div"]=nu.cid.map(te_div); nu["uniqFPM"]=nu.cid.map(uq)
print(f"\nTE-overlapping: {nu.te_family.notna().sum()}/{len(nu)} ({100*nu.te_family.notna().mean():.0f}%)")
print("dominant TE families:");
for fam,n in nu.te_family.value_counts().head(8).items(): print(f"    {n:3d}  {fam}")
print(f"\nTE age (RM %div): near-universal median {nu.te_div.median():.1f}%  (cf. all-non-ref 14.3%, reference 16.4%)")
print(f"expression (uniqFPM): near-universal median {nu.uniqFPM.median():.1f}  max {nu.uniqFPM.max():.0f}")
print(f"  high-expression (uniqFPM>50): {(nu.uniqFPM>50).sum()} entries")
# representative per locus (highest-expression member)
reps=[]
for k,members in comp.items():
    m=nu[nu.cid.isin(members)].sort_values("uniqFPM",ascending=False).iloc[0]
    reps.append(dict(locus=k,n_strains=len(set(c.split('|')[0] for c in members)),rep_cid=m.cid,strain=m.strain,chrom=m.chrom,start=int(m.start),end=int(m.end),
                     kb=round((m.end-m.start)/1000,1),te_family=m.te_family,te_div=m.te_div,uniqFPM=round(m.uniqFPM,1)))
rep=pd.DataFrame(reps).sort_values("uniqFPM",ascending=False)
print(f"\n=== {len(rep)} distinct loci (top by expression) ===")
print(rep[["n_strains","strain","chrom","start","kb","te_family","te_div","uniqFPM"]].head(20).to_string(index=False))
rep.to_csv(f"{D}/shared_subset_loci.csv",index=False); nu.to_csv(f"{D}/shared_subset_entries.csv",index=False)
print("\nwrote shared_subset_loci.csv + shared_subset_entries.csv")
