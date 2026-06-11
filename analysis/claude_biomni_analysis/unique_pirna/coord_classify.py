#!/usr/bin/env python3
"""Coordinate-based TE-driven test (per strain X), removing the sequence-match confound: does a
candidate's PRODUCTION LOCUS in X's own genome fall WITHIN a strain-private insertion's X-location?
Insertion X-loci = minimap2 of X-private-insertion sequences back to X's genome (primary, mapq>=20).
Candidate loci = cand_self BAM (X coords, PanSN — same frame as the genome FASTA, no stripping). Reports
locus-in-private-insertion by Step-4 class (other classes = null) and the TE-driven count (strain-private
locus + locus in private insertion + TE-annotated)."""
import os,subprocess,sys,tempfile
from collections import Counter
import pandas as pd, pysam
S4="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
X=sys.argv[1]
d=pd.read_csv(f"{S4}/{X}.step4_classified.csv.gz")
te=pd.read_csv(f"{S4}/{X}.private_TE_percandidate.csv.gz"); id2te=dict(zip(te.id,te.classfam))
# candidate loci BED from cand_self (all alignments)
bam=pysam.AlignmentFile(f"{S4}/{X}.cand_self.Aligned.sortedByCoord.out.bam","rb")
cb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG)
for a in bam.fetch(until_eof=True):
    if not a.is_unmapped:
        cb.write(f"{a.reference_name}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
cb.close(); bam.close()
out=subprocess.run([BT,"intersect","-a",cb.name,"-b",f"{PG}/{X}.ins_loci.bed","-wa","-u"],
                   capture_output=True,text=True).stdout
os.unlink(cb.name)
in_ins=set(l.split("\t")[3] for l in out.splitlines() if l)
rows=[]
for klass in sorted(d.klass.unique()):
    ids=set(d.loc[d.klass==klass,"id"]); hit=ids & in_ins
    rows.append(dict(strain=X,klass=klass,n=len(ids),locus_in_priv_ins=len(hit),
                     pct=round(100*len(hit)/max(len(ids),1),3)))
sp=set(d.loc[d.klass=="unique: strain-private locus","id"])
tedrv=[i for i in (sp & in_ins) if i in id2te]; fam=Counter(id2te[i] for i in tedrv)
pd.DataFrame(rows).to_csv(f"{PG}/{X}.coord_byclass.csv",index=False)
print(f"[{X}] locus-in-private-insertion by class:")
for r in rows: print("   ",r["klass"],"->",r["locus_in_priv_ins"],"/",r["n"],f"({r['pct']}%)")
print(f"[{X}] TE-DRIVEN (strain-private locus + locus in private insertion + TE): {len(tedrv)}")
print(f"       families: {dict(fam.most_common(8))}")
