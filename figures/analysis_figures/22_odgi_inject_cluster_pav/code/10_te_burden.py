#!/usr/bin/env python3
"""Per-strain NON-REFERENCE TE-insertion burden: count, from the deconstructed pangenome VCF, the TE-sized
(>=300 bp longer than REF) INSERTION alleles each strain carries relative to GRCm39 = the strain's genome-wide
non-reference TE-insertion load (a direct 'TE evolution' measure, independent of the piRNA clusters)."""
import gzip, json
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
VCF=f"{B}/results/pangenome/output/mouse_17strain_pangenome.vcf.gz"
D=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
cnt={s:0 for s in S}; sidx={}; n=0
with gzip.open(VCF,"rt") as f:
    for ln in f:
        if ln[0]=="#":
            if ln.startswith("#CHROM"):
                cols=ln.rstrip("\n").split("\t"); sidx={s:cols.index(s) for s in S if s in cols}
            continue
        n+=1; F=ln.rstrip("\n").split("\t"); ref=F[3]; alts=F[4].split(",")
        te=set(str(i+1) for i,a in enumerate(alts) if len(a)-len(ref)>=300)
        if not te: continue
        for s,j in sidx.items():
            for al in F[j].split(":")[0].replace("|","/").split("/"):
                if al in te: cnt[s]+=1; break
json.dump(cnt,open(f"{D}/te_insertion_burden.json","w"))
print(f"variants scanned: {n:,}")
for s in sorted(cnt,key=lambda x:-cnt[x]): print(f"  {s:13s} {cnt[s]:6d} non-reference TE-sized insertions")
