#!/usr/bin/env python3
"""THEME 22 step 21 — POSITIVE CONTROLS for the graph check. Step-20's first odgi pav gave GRCm39=0 for ALL 19
shared-subset loci (alarming, given C57BL_6NJ~GRCm39). To decide 'broken GRCm39 group' vs 'real absence', add two
coordinate-grounded control tiers and re-run odgi pav (graph_check_loci3.bed -> graph_check_pav3.tsv via the step-20
command):
  (1) GRCm39-FRAME loci  (regions in GRCm39's own coords) -> GRCm39 MUST cover ~1.0,
  (2) C57BL_6NJ clusters that DID lift to GRCm39 (= present) -> GRCm39 should cover ~1.0.
RESULT: frame 1.0, lifted 1.0, shared subset 0.0 -> the GRCm39 group is correct and the shared subset is GENUINELY
absent from GRCm39's path (the actual pangenome MSA). This corrects step-18's minimap2 'present' calls, which were the
TE-sequence-matches-elsewhere confound."""
import collections
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"; CP=f"{B}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
inj=collections.defaultdict(list); clu_len=collections.Counter()
for l in open(f"{D}/inject.bed"):
    f=l.rstrip("\n").split('\t'); p=f[0].split('#'); off=int(p[3])
    inj[(p[0],p[2])].append((off+int(f[1]),off+int(f[2]),f[3])); clu_len[f[3]]+=int(f[2])-int(f[1])
nonref=set((l.split('\t')[0],int(l.split('\t')[1]),int(l.split('\t')[2])) for l in open(f"{D}/nonref/C57BL_6NJ.nonref.bed"))
ctrl=[]
for l in open(f"{CP}/C57BL_6NJ.clusters.bed"):
    c=l.rstrip("\n").split('\t'); chrom,s,e=c[0],int(c[1]),int(c[2])
    if (chrom,s,e) in nonref or not (3000<e-s<12000): continue
    cs=[(ns,ne,cl) for ns,ne,cl in inj.get(('C57BL_6NJ',chrom),[]) if not(ne<s or ns>e)]
    if cs: b=max(cs,key=lambda x:min(x[1],e)-max(x[0],s)); ctrl.append((b[2],clu_len[b[2]],f"CTRL_lifted_{chrom}:{s}"))
    if len(ctrl)>=5: break
with open(f"{D}/graph_check_loci3.bed","w") as o:
    o.write("GRCm39#0#1\t3000000\t3005000\tGRCM39FRAME_chr1_3Mb\n")
    o.write("GRCm39#0#11\t60000000\t60008000\tGRCM39FRAME_chr11_60Mb\n")
    for clu,L,nm in ctrl: o.write(f"{clu}\t0\t{L}\t{nm}\n")
    o.write(open(f"{D}/graph_check_loci.bed").read())
print(f"graph_check_loci3.bed: 2 GRCm39-frame + {len(ctrl)} C57BL_6NJ-lifted controls + 19 subset")
print("re-run step-20 odgi pav with graph_check_loci3.bed -> graph_check_pav3.tsv")
