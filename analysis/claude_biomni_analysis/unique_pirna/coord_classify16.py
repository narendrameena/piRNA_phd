#!/usr/bin/env python3
"""16-strain coordinate TE-driven classify (per strain X): does each candidate's production locus in X's own
genome fall within an X-private insertion locus? cand_self16 BAM (PanSN) intersect ins16/{X}.ins_loci.bed
(minimap2 of X-private insertions back to X, PanSN). Per-class locus-in-private-insertion + the random-locus
null EXP = merged private-insertion bp / genome. Class from final_classified (id = X|tp|seq)."""
import os,subprocess,sys,tempfile
import pandas as pd, pysam
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; X=sys.argv[1]
d=pd.read_csv(f"{U}/unique16/final_classified.csv.gz"); d=d[d.strain==X].copy()
d["id"]=X+"|"+d.timepoint.astype(str)+"|"+d.sequence
bam=pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam","rb")
cb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG)
for a in bam.fetch(until_eof=True):
    if not a.is_unmapped: cb.write(f"{a.reference_name}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
cb.close(); bam.close()
out=subprocess.run([BT,"intersect","-a",cb.name,"-b",f"{PG}/ins16/{X}.ins_loci.bed","-wa","-u"],capture_output=True,text=True).stdout
os.unlink(cb.name)
in_ins=set(l.split("\t")[3] for l in out.splitlines() if l)
rows=[]
for klass in sorted(d.klass.unique()):
    ids=set(d.loc[d.klass==klass,"id"]); hit=ids & in_ins
    rows.append(dict(strain=X,klass=klass,n=len(ids),locus_in_priv_ins=len(hit),pct=round(100*len(hit)/max(len(ids),1),3)))
mrg=subprocess.run(f"sort -k1,1 -k2,2n {PG}/ins16/{X}.ins_loci.bed | {BT} merge",shell=True,capture_output=True,text=True).stdout
insbp=sum(int(p[2])-int(p[1]) for p in (l.split('\t') for l in mrg.splitlines()) if len(p)>=3)
gsize=sum(int(l) for c,l in (ln.split() for ln in open(f"{ROOT}/results/indexs/{X}/chrNameLength.txt")))
exp=round(100*insbp/gsize,4)
pd.DataFrame(rows).assign(exp_pct=exp).to_csv(f"{PG}/{X}.coord_byclass16.csv",index=False)
print(f"[{X}] EXP(null)={exp}% (priv-ins {insbp/1e6:.1f} Mb / genome {gsize/1e9:.2f} Gb)")
for r in rows: print(f"   {r['klass']} -> {r['locus_in_priv_ins']}/{r['n']} ({r['pct']}%) fold={round(r['pct']/exp,2) if exp else 'NA'}")
