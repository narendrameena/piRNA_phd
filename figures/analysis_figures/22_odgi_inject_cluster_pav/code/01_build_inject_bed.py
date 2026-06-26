#!/usr/bin/env python3
"""THEME 22 step 1 — map each strain's NATIVE piCB cluster coords onto the fragmented graph strain paths
(STRAIN#0#chr#offset) and emit an odgi-inject BED (path-space: path_name, local_start, local_end, cluster_name).
Fragment chr-span = [offset, offset+length); offset from the path name, length from vg paths -E. A native cluster
(chr,cs,ce) -> the fragment with offset<=cs and ce<=offset+length; local = (cs-offset, ce-offset). Clusters spanning a
fragment boundary / assembly gap / non-graph chrom are dropped (rate reported)."""
import pandas as pd, bisect, os
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
STRAINS=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
# fragment map: strain -> chrom -> sorted [(offset, end, path_name)]
frag={}
for ln in open(f"{D}/graph_path_lengths.txt"):
    parts=ln.rstrip("\n").split("\t")
    if len(parts)!=2: continue
    name,length=parts[0],int(parts[1]); p=name.split("#")
    if len(p)<4: continue   # GRCm39 single-path-per-chrom (no offset) — not an inject target
    frag.setdefault(p[0],{}).setdefault(p[2],[]).append((int(p[3]),int(p[3])+length,name))
for s in frag:
    for c in frag[s]: frag[s][c].sort()
rows=[]; mapped=dropped_nochrom=dropped_boundary=0
for X in STRAINS:
    cl=pd.read_csv(f"{CP}/{X}.clusters.bed",sep="\t",header=None,usecols=[0,1,2],names=["chrom","start","end"],dtype={"chrom":str})
    fg=frag.get(X,{})
    for i,r in enumerate(cl.itertuples(index=False)):
        c=str(r.chrom); cs=int(r.start); ce=int(r.end); arr=fg.get(c)
        if not arr: dropped_nochrom+=1; continue
        j=bisect.bisect_right([a[0] for a in arr],cs)-1; found=None
        for k in (j,j-1,j+1):
            if 0<=k<len(arr) and arr[k][0]<=cs and ce<=arr[k][1]: found=arr[k]; break
        if found is None: dropped_boundary+=1; continue
        off,_,name=found; rows.append((name,cs-off,ce-off,f"CLU|{X}|{i}")); mapped+=1
tot=mapped+dropped_nochrom+dropped_boundary
print(f"clusters: mapped={mapped:,} | dropped no-graph-chrom={dropped_nochrom:,} | dropped fragment-boundary={dropped_boundary:,} | drop-rate={100*(dropped_nochrom+dropped_boundary)/tot:.1f}%")
df=pd.DataFrame(rows,columns=["path","start","end","name"]).sort_values(["path","start"])
df.to_csv(f"{D}/inject.bed",sep="\t",header=False,index=False)
df[["name"]].to_csv(f"{D}/cluster_paths.txt",header=False,index=False)
print(f"wrote inject.bed ({len(df):,} records) + cluster_paths.txt | per-strain mapped:",{X:int((df.name.str.split('|').str[1]==X).sum()) for X in STRAINS[:4]},"...")
