#!/usr/bin/env python3
"""THEME 22 step 18 — WHY are the 'present-in-most' loci absent from GRCm39? For each distinct locus: align the unique
FLANKS (+/-5kb) and the CLUSTER body to GRCm39 (minimap2). Flanks anchor the homologous locus (unique seq, robust);
the cluster is TE-rich so we require it to align AT the flank locus (not TE copies elsewhere). N-check the bracket.
  flanks don't align            -> DIVERGENT_absent (locus not in GRCm39)
  flanks align, cluster at locus-> PRESENT_in_GRCm39 (technical lift failure / not truly novel)
  flanks align, cluster absent, bracket has Ns  -> ASSEMBLY_GAP
  flanks align, cluster absent, bracket clean   -> C57BL6J_DELETION (reference lineage lost a cluster the panel kept)"""
import pandas as pd, subprocess, os
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"; SIF=f"{B}/cactus_v2.9.3.sif"
GR=f"{B}/results/pangenome/prepared/GRCm39.fa"; MMI="/tmp/GRCm39.mmi"; PREP=lambda X:f"{B}/results/pangenome/prepared/{X}.fa"
def sx(c): return subprocess.run(f"singularity exec --bind /mnt,/tmp {SIF} bash -c '{c}'",shell=True,capture_output=True,text=True).stdout
if not os.path.exists(MMI): sx(f"minimap2 -x asm10 -d {MMI} {GR}")
def aln(qfa):  # -> best query cov, list of (tchrom,tstart,tend)
    cov=0.0; tp=[]
    for ln in sx(f"minimap2 -x asm10 {MMI} {qfa} 2>/dev/null").splitlines():
        f=ln.split('\t'); cov=max(cov,(int(f[3])-int(f[2]))/int(f[1])); tp.append((f[5],int(f[7]),int(f[8])))
    return cov,tp
def fracN(ch,a,b):
    a=max(1,a); s="".join(sx(f"samtools faidx {GR} {ch}:{a}-{b}").splitlines()[1:]).upper(); return (s.count('N')/len(s)) if s else 1.0
rep=pd.read_csv(f"{D}/shared_subset_loci.csv"); rows=[]
for _,r in rep.iterrows():
    X,ch,st,en=r.strain,str(r.chrom),int(r.start),int(r.end); sz=en-st
    sx(f"samtools faidx {PREP(X)} {ch}:{st}-{en} > /tmp/clu.fa")
    sx(f"samtools faidx {PREP(X)} {ch}:{max(1,st-5000)}-{st} {ch}:{en}-{en+5000} > /tmp/flk.fa")
    cc,ct=aln("/tmp/clu.fa"); fc,ft=aln("/tmp/flk.fa")
    cls="DIVERGENT_absent"; brk=""; fN=""
    if fc>0.4 and ft:
        chrom=max(set(t[0] for t in ft),key=lambda c:sum(1 for t in ft if t[0]==c))
        ps=[p for t in ft if t[0]==chrom for p in (t[1],t[2])]; a,b=min(ps),max(ps); brk=f"{chrom}:{a}-{b}"
        at_locus=any(t[0]==chrom and t[2]>a-sz-5000 and t[1]<b+sz+5000 for t in ct) and cc>0.4
        if at_locus: cls="PRESENT_in_GRCm39"
        else:
            fN=round(fracN(chrom,a,b),2); cls="ASSEMBLY_GAP" if fN>0.3 else "C57BL6J_DELETION"
    rows.append(dict(strain=X,n_strains=r.n_strains,chrom=ch,start=st,kb=r.kb,te_family=r.te_family,uniqFPM=r.uniqFPM,
                     qcov_flank=round(fc,2),qcov_cluster=round(cc,2),bracket=brk,fracN=fN,classification=cls))
res=pd.DataFrame(rows).sort_values("uniqFPM",ascending=False)
print(res.to_string(index=False)); print("\n",res.classification.value_counts().to_string())
res.to_csv(f"{D}/shared_subset_grcm39_status.csv",index=False); print("\nwrote shared_subset_grcm39_status.csv")
