#!/usr/bin/env python3
"""16-strain coordinate TE-driven classify (per strain X): does each candidate's production locus in X's own
genome fall within an X-private insertion locus? cand_self16 BAM (PanSN) intersect ins16/{X}.ins_loci.bed
(minimap2 of X-private insertions back to X, PanSN). Per-class locus-in-private-insertion vs TWO nulls: the
single-locus null EXP = merged private-insertion bp / genome, AND a MULTIPLICITY-MATCHED null 1-(1-p)^NH. HEADLINE =
the NH==1 uniquely-mapping subset (pct_uniqmap/n_uniqmap; fold_uniqmap = pct_uniqmap/EXP), the clean, un-inflatable
signal whose single production locus is unambiguously localizable. NH>1 multimappers are un-localizable (no single
production locus) and are kept ONLY in the matched-null and naive columns (fold_matched on the 1-(1-p)^NH null,
fold on the naive single-locus null over ALL candidates). Class from final_classified_clean_2read klass5 (id = X|tp|seq)."""
import os,subprocess,sys,tempfile
import pandas as pd, pysam, numpy as np
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; X=sys.argv[1]
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz"); d=d[d.strain==X].copy()   # ADOPTED ≥2-read 5-class
d["id"]=X+"|"+d.timepoint.astype(str)+"|"+d.sequence
bam=pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam","rb")
cb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG)
nh={}   # id -> mapping multiplicity (NH). KEEP all candidates: TE-derived piRNAs are multi-copy (young/active
        # insertions have near-identical copies -> high NH, and are the MOST insertion-associated), so restricting to
        # NH==1 discards the very TE drivers. Instead correct the null for each candidate's k chances (exp_matched).
for a in bam.fetch(until_eof=True):
    if a.is_unmapped: continue
    nh[a.query_name]=max(nh.get(a.query_name,1), a.get_tag("NH") if a.has_tag("NH") else 1)
    cb.write(f"{a.reference_name}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
cb.close(); bam.close()
out=subprocess.run([BT,"intersect","-a",cb.name,"-b",f"{PG}/ins16/{X}.ins_loci.bed","-wa","-u"],capture_output=True,text=True).stdout
os.unlink(cb.name)
in_ins=set(l.split("\t")[3] for l in out.splitlines() if l)
mrg=subprocess.run(f"sort -k1,1 -k2,2n {PG}/ins16/{X}.ins_loci.bed | {BT} merge",shell=True,capture_output=True,text=True).stdout
insbp=sum(int(pp[2])-int(pp[1]) for pp in (l.split('\t') for l in mrg.splitlines()) if len(pp)>=3)
gsize=sum(int(l) for c,l in (ln.split() for ln in open(f"{ROOT}/results/indexs/{X}/chrNameLength.txt")))
p=insbp/gsize; exp=round(100*p,4)
rows=[]
for klass in sorted(d.klass5.unique()):
    ids=set(d.loc[d.klass5==klass,"id"]); hit=ids & in_ins
    kk=np.array([nh.get(i,1) for i in ids]) if ids else np.array([1])   # per-candidate mapping multiplicity (NH)
    exp_m=100*float((1-(1-p)**kk).mean())                               # MULTIPLICITY-MATCHED null: P(>=1 of k random loci in insertion)
    pct=round(100*len(hit)/max(len(ids),1),3)
    idu=[i for i in ids if nh.get(i,1)==1]; hitu=set(idu)&in_ins        # NH==1 uniquely-mapping subset = the CLEAN, un-inflatable signal (headline)
    pct_u=round(100*len(hitu)/max(len(idu),1),3)
    rows.append(dict(strain=X,klass=klass,n=len(ids),n_uniqmap=len(idu),locus_in_priv_ins=len(hit),
                     pct=pct,pct_uniqmap=pct_u,
                     exp_pct=exp,exp_matched_pct=round(exp_m,4),
                     fold=round(pct/exp,2) if exp else None,               # naive single-locus null on ALL candidates (INFLATED by multi-mapping)
                     fold_matched=round(pct/exp_m,2) if exp_m else None,   # multiplicity-matched null on ALL (inclusive but blends un-localizable NH>1)
                     fold_uniqmap=round(pct_u/exp,2) if exp else None))    # HEADLINE: NH==1 uniquely-mapping fold (clean, unambiguous)
pd.DataFrame(rows).to_csv(f"{PG}/{X}.coord_byclass16.csv",index=False)
print(f"[{X}] EXP(single)={exp}%  priv-ins {insbp/1e6:.1f} Mb / genome {gsize/1e9:.2f} Gb")
for r in rows: print(f"   {r['klass']} -> NH==1 {r['pct_uniqmap']}% fold={r['fold_uniqmap']} (n={r['n_uniqmap']})  |  ALL {r['pct']}% fold_naive={r['fold']} fold_matched={r['fold_matched']} (n={r['n']})")
