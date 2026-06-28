#!/usr/bin/env python3
"""THEME 22 step 19 — set up the graph-native cross-check (step 20) of the shared-subset deletion/present calls.
Map each shared-subset representative cluster to its injected CLU path by RECONSTRUCTED native coords from inject.bed
(native = fragment_offset[4th PanSN field] + local); this also confirms the inject's CLU|strain|i == clusters.bed row i
(i = re-lift NR-1). Write the cluster-path loci bed (CLU|strain|id over [0,len]) and the path-groups TSV (GRCm39 + the
16 strains, from graph_path_lengths.txt) for `odgi pav`."""
import pandas as pd, collections
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
inj=collections.defaultdict(list); clu_len=collections.Counter()
for l in open(f"{D}/inject.bed"):
    f=l.rstrip("\n").split('\t'); p=f[0].split('#'); off=int(p[3]); ns=off+int(f[1]); ne=off+int(f[2])
    inj[(p[0],p[2])].append((ns,ne,f[3])); clu_len[f[3]]+=int(f[2])-int(f[1])
rep=pd.read_csv(f"{D}/shared_subset_loci.csv"); bed=[]
for _,r in rep.iterrows():
    cs=[(ns,ne,c) for ns,ne,c in inj.get((r.strain,str(r.chrom)),[]) if not(ne<r.start or ns>r.end)]
    if cs:
        b=max(cs,key=lambda c:min(c[1],r.end)-max(c[0],r.start)); bed.append((b[2],clu_len[b[2]],r.rep_cid))
with open(f"{D}/graph_check_loci.bed","w") as o:
    for clu,L,cid in bed: o.write(f"{clu}\t0\t{L}\t{cid}\n")
with open(f"{D}/graph_check_groups.tsv","w") as o:
    for l in open(f"{D}/graph_path_lengths.txt"):
        pth=l.split('\t')[0]; o.write(f"{pth}\t{'GRCm39' if pth.startswith('GRCm39') else pth.split('#')[0]}\n")
print(f"mapped {len(bed)}/{len(rep)} reps -> graph_check_loci.bed; groups from graph_path_lengths.txt")
