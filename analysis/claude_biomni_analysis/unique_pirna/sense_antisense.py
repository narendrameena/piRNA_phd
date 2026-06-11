#!/usr/bin/env python3
"""Sense/antisense of unique piRNAs relative to the overlapping TE (per strain X). piRNA locus+strand
(cand_self) vs stranded TE annotation (RM .out): SAME strand = sense (from the TE's own transcript),
OPPOSITE = antisense (silencing-competent — can base-pair the TE mRNA). Reports antisense fraction by
Step-4 class and by TE family. Coords are PanSN in both (no stripping). Usage: <strain>"""
import os,subprocess,tempfile
from collections import defaultdict
import pandas as pd, pysam
S4="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
SA="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/sense_antisense"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
X=os.sys.argv[1]
d=pd.read_csv(f"{S4}/{X}.step4_classified.csv.gz"); id2class=dict(zip(d.id,d.klass))
# piRNA loci BED (strand-aware) from cand_self
bam=pysam.AlignmentFile(f"{S4}/{X}.cand_self.Aligned.sortedByCoord.out.bam","rb")
pb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=SA)
for a in bam.fetch(until_eof=True):
    if a.is_unmapped: continue
    pb.write(f"{a.reference_name}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\t0\t{'-' if a.is_reverse else '+'}\n")
pb.close(); bam.close()
# intersect with stranded TE BED (-wo: a=6 cols, b=6 cols, overlap=last)
out=subprocess.run([BT,"intersect","-a",pb.name,"-b",f"{SA}/{X}.TE_stranded.bed","-wo"],capture_output=True,text=True).stdout
os.unlink(pb.name)
# per candidate: primary TE by max overlap; sense/antisense vs that TE
best=defaultdict(lambda:(-1,None,None))   # id -> (overlap, family, orientation)
for ln in out.splitlines():
    f=ln.split("\t")
    cid=f[3]; p_str=f[5]; fam=f[9].split("|")[-1]; t_str=f[11]; ov=int(f[-1])
    orient="sense" if p_str==t_str else "antisense"
    if ov>best[cid][0]: best[cid]=(ov,fam,orient)
rec=[]
for cid,(ov,fam,orient) in best.items():
    rec.append(dict(id=cid,klass=id2class.get(cid,"NA"),family=fam,orientation=orient))
r=pd.DataFrame(rec)
r.to_csv(f"{SA}/{X}.sense_antisense_percandidate.csv.gz",index=False,compression="gzip")
# antisense fraction by class
print(f"[{X}] TE-overlapping candidates: {len(r):,}")
by=r.groupby("klass").orientation.apply(lambda s:(s=="antisense").mean()*100).round(1)
n =r.groupby("klass").size()
print("  antisense % by Step-4 class (n):")
for k in by.index: print(f"    {k}: {by[k]}%  (n={n[k]})")
# antisense fraction for genuinely-unique, by top family
gu=r[r.klass.str.startswith("unique")]
fam=gu.groupby("family").orientation.apply(lambda s:(s=="antisense").mean()*100).round(1)
fn =gu.groupby("family").size(); top=fn.sort_values(ascending=False).head(8).index
print("  genuinely-unique antisense % by top TE family:")
for k in top: print(f"    {k}: {fam[k]}%  (n={fn[k]})")
pd.DataFrame({"klass":by.index,"antisense_pct":by.values,"n":n.reindex(by.index).values}).to_csv(f"{SA}/{X}.antisense_byclass.csv",index=False)
