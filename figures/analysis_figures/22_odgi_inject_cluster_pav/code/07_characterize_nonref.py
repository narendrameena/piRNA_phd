#!/usr/bin/env python3
"""THEME 22 step 7 — characterise the 1,393 NON-REFERENCE piRNA clusters (sequence absent from GRCm39). Are they
TE-insertion-driven? Intersect each strain's non-reference cluster BED (native coords) with its RepeatMasker .out
(own-genome TE annotation). Expectation (mirroring theme 21): young L1/IAP/ERV insertions = strain-specific new
piRNA-cluster loci absent from the reference."""
import pandas as pd, glob, subprocess, os, collections
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
fam=collections.Counter(); tot=hit=0; per=[]
for X in S:
    f=f"{D}/nonref/{X}.nonref.bed"
    if not os.path.exists(f) or os.path.getsize(f)==0: continue
    nr=pd.read_csv(f,sep="\t",header=None,names=["chrom","start","end","id"],dtype={"chrom":str})
    rmf=glob.glob(f"{B}/resources/repeatMasker/{X}_*.out")
    if not rmf: per.append((X,len(nr),None)); tot+=len(nr); continue
    nrf=f"/tmp/nr_{X}.bed"
    nr.assign(chrom=X+"#1#chr"+nr.chrom.astype(str)).to_csv(nrf,sep="\t",header=False,index=False)
    cmd=f"awk 'NR>3{{print $5\"\\t\"$6\"\\t\"$7\"\\t\"$11}}' {rmf[0]} | sort -k1,1 -k2,2n | bedtools intersect -a <(sort -k1,1 -k2,2n {nrf}) -b stdin -wa -wb 2>/dev/null"
    out=subprocess.run(cmd,shell=True,capture_output=True,text=True,executable="/bin/bash").stdout
    ids=set();
    for ln in out.splitlines():
        c=ln.split("\t"); ids.add(c[3]); fam[c[7]]+=1
    tot+=len(nr); hit+=len(ids); per.append((X,len(nr),len(ids))); os.remove(nrf)
    print(f"  {X}: {len(nr)} non-ref | TE-overlapping {len(ids)} ({100*len(ids)/max(len(nr),1):.0f}%)")
print(f"\n=== NON-REFERENCE clusters: {tot} total | TE-overlapping {hit} ({100*hit/max(tot,1):.0f}%) ===")
print("top TE families at non-reference piRNA clusters:")
for f2,n in fam.most_common(12): print(f"    {n:5d}  {f2}")
pd.DataFrame(per,columns=["strain","n_nonref","n_TE_overlap"]).to_csv(f"{D}/nonref_te_summary.csv",index=False)
pd.DataFrame(fam.most_common(),columns=["TE_family","n"]).to_csv(f"{D}/nonref_te_families.csv",index=False)
