#!/usr/bin/env python3
"""Which TE families do the GENUINELY-UNIQUE strain-PRIVATE-LOCUS piRNAs originate from?
For each strain X: take Step-4 'unique: strain-private locus' candidates -> their loci in X's own
genome (cand_self BAM, mm0) -> strip PanSN prefix ({X}#1#chrN -> chrN) to match X's RepeatMasker BED
(col1=chrN, col4=TE_name|class/family) -> bedtools intersect -> per candidate assign the TE family with
the largest overlap. Reports TE-derived fraction, class breakdown, and top families. (Sense/antisense
needs the TE GFF3 with strand and is a Step-6 refinement; here = family identity only.)
"""
import os,subprocess,sys,tempfile
import pandas as pd, pysam
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
S4=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/step4"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
STR=["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]
allfam=[]; allcls=[]; summ=[]
for X in STR:
    rm=f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
    d=pd.read_csv(f"{S4}/{X}.step4_classified.csv.gz")
    priv=set(d.loc[d.klass=="unique: strain-private locus","id"])
    pref=f"{X}#1#"
    # write candidate-locus BED (chrN), one line per alignment of a private candidate
    bam=pysam.AlignmentFile(f"{S4}/{X}.cand_self.Aligned.sortedByCoord.out.bam","rb")
    cbed=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=S4)
    n_aln=0; mapped_ids=set()
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.query_name not in priv: continue
        chrom=a.reference_name
        chrom=chrom[len(pref):] if chrom.startswith(pref) else chrom
        st="-" if a.is_reverse else "+"
        cbed.write(f"{chrom}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\t0\t{st}\n")
        n_aln+=1; mapped_ids.add(a.query_name)
    cbed.close(); bam.close()
    # intersect with RM (report candidate id + TE col4 + overlap bp)
    out=subprocess.run([BT,"intersect","-a",cbed.name,"-b",rm,"-wo"],capture_output=True,text=True).stdout
    os.unlink(cbed.name)
    # parse: a=6 cols, b=7 cols, -wo appends overlap bp as last col; b col4 (TE_name|class/family)=index 9
    rows=[]
    for ln in out.splitlines():
        f=ln.split("\t")
        cid=f[3]; te_field=f[9]; ov=int(f[-1])
        cf=te_field.split("|")[-1]                  # class/family
        rows.append((cid,cf,ov))
    ov=pd.DataFrame(rows,columns=["id","classfam","ov"])
    te_ids=set(ov.id)
    # per candidate -> family with max total overlap
    prim=(ov.groupby(["id","classfam"]).ov.sum().reset_index()
            .sort_values("ov",ascending=False).drop_duplicates("id"))
    prim[["id","classfam"]].to_csv(f"{S4}/{X}.private_TE_percandidate.csv.gz",index=False,compression="gzip")
    famc=prim.classfam.value_counts()
    clsc=prim.classfam.str.split("/").str[0].value_counts()
    npriv=len(priv)
    summ.append(dict(strain=X,private=npriv,mapped=len(mapped_ids),
                     TE_derived=len(te_ids),TE_frac=round(100*len(te_ids)/npriv,1)))
    for k,v in famc.items(): allfam.append(dict(strain=X,classfam=k,n=int(v)))
    for k,v in clsc.items(): allcls.append(dict(strain=X,te_class=k,n=int(v)))
    print(f"\n[{X}] private={npriv:,} mapped-to-own-genome={len(mapped_ids):,} TE-derived={len(te_ids):,} ({100*len(te_ids)/npriv:.1f}%)")
    print("  top TE classes:",clsc.head(6).to_dict())
    print("  top TE families:",famc.head(10).to_dict())

pd.DataFrame(summ).to_csv(f"{S4}/TE_private_summary.csv",index=False)
pd.DataFrame(allfam).to_csv(f"{S4}/TE_private_families.csv",index=False)
pd.DataFrame(allcls).to_csv(f"{S4}/TE_private_classes.csv",index=False)
print("\nwrote TE_private_summary.csv, TE_private_families.csv, TE_private_classes.csv")
