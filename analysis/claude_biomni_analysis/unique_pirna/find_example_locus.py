#!/usr/bin/env python3
"""Find a clean, REAL example locus for the schematic: a strain-private TE insertion that produces
several strain-private piRNAs (a TE-driven piRNA source). For strain X: strain-private-locus candidates'
own-genome loci (cand_self) that fall inside an X-private insertion (ins_loci) AND a TE (stranded RM .out)
-> group by insertion -> pick an insertion carrying many strain-private piRNAs, report TE family/strand
and the piRNA strands (antisense-to-TE = silencing). Usage: <strain>"""
import os,subprocess,sys,tempfile
from collections import defaultdict
import pandas as pd, pysam
S4="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
SA="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/sense_antisense"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
X=sys.argv[1] if len(sys.argv)>1 else "SPRET_EiJ"
d=pd.read_csv(f"{S4}/{X}.step4_classified.csv.gz"); sp=set(d.loc[d.klass=="unique: strain-private locus","id"])
# strain-private piRNA loci (all alignments) -> BED
bam=pysam.AlignmentFile(f"{S4}/{X}.cand_self.Aligned.sortedByCoord.out.bam","rb")
pb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG)
for a in bam.fetch(until_eof=True):
    if a.is_unmapped or a.query_name not in sp: continue
    pb.write(f"{a.reference_name}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\t0\t{'-' if a.is_reverse else '+'}\n")
pb.close(); bam.close()
# intersect piRNA loci with private insertions (which insertion each piRNA sits in)
ins=subprocess.run([BT,"intersect","-a",pb.name,"-b",f"{PG}/{X}.ins_loci.bed","-wa","-wb"],capture_output=True,text=True).stdout
# group: insertion (chrom,start,end of ins) -> list of (piRNA pos, strand)
byins=defaultdict(list)
for ln in ins.splitlines():
    f=ln.split("\t"); pirna=(f[0],int(f[1]),int(f[2]),f[5],f[3]); ins_key=(f[6],int(f[7]),int(f[8]))
    byins[ins_key].append(pirna)
# rank insertions by # distinct piRNAs
ranked=sorted(byins.items(),key=lambda kv:-len({p[4] for p in kv[1]}))
# annotate top insertions with the overlapping TE (stranded RM) and pick a clean one
te=pd.read_csv(f"{SA}/{X}.TE_stranded.bed",sep="\t",header=None,names=["c","s","e","fam","sc","st"]) if os.path.exists(f"{SA}/{X}.TE_stranded.bed") else None
print(f"[{X}] strain-private piRNAs in private insertions; top insertions by #piRNAs:")
chosen=None
for (c,s,e),plist in ranked[:12]:
    npi=len({p[4] for p in plist}); strands=[p[3] for p in plist]
    tef="NA"; test="NA"
    if te is not None:
        sub=te[(te.c==c)&(te.s<e)&(te.e>s)]
        if len(sub): row=sub.iloc[(sub.e.clip(upper=e)-sub.s.clip(lower=s)).values.argmax()]; tef=row.fam; test=row.st
    plus=strands.count("+"); minus=strands.count("-")
    print(f"  {c}:{s}-{e} ({e-s}bp) | piRNAs={npi} | strands +{plus}/-{minus} | TE={tef} ({test})")
    if chosen is None and npi>=4 and tef!="NA": chosen=((c,s,e),plist,tef,test)
if chosen:
    (c,s,e),plist,tef,test=chosen
    import json
    pos=sorted({(p[1],p[2],p[3]) for p in plist})
    ex=dict(strain=X,chrom=c,ins_start=s,ins_end=e,ins_bp=e-s,n_pirna=len({p[4] for p in plist}),
            te_family=tef,te_strand=test,pirna_positions=[list(p) for p in pos])
    json.dump(ex,open(f"{PG}/example_locus.json","w"),indent=1)
    print(f"\nCHOSEN example -> example_locus.json: {c}:{s}-{e}, TE={tef}({test}), {ex['n_pirna']} strain-private piRNAs")
