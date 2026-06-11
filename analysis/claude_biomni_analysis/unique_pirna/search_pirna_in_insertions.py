#!/usr/bin/env python3
"""TE-DRIVEN test (pangenome): is each Step-4 candidate's sequence CONTAINED in an X-private pangenome
insertion (+/- revcomp)? The four Step-4 classes serve as the built-in null: strain-private-locus
piRNAs should be embedded in strain-private insertions far more than expressed-elsewhere / SNP-variant /
conserved-but-silent ones. Hits are crossed with the per-candidate TE family (Step 6a) -> TE-driven.
Substring search via grep -F -o (Aho-Corasick); counts are a slight lower bound (non-overlapping matches).
"""
import os,subprocess,tempfile
import pandas as pd
S4="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
comp=str.maketrans("ACGTN","TGCAN")
def rc(s): return s.translate(comp)[::-1]
rate=[]; tedrv=[]
for X in ["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]:
    ins=f"{PG}/{X}.private_insertions.fasta"
    d=pd.read_csv(f"{S4}/{X}.step4_classified.csv.gz")
    seqs=d.sequence.unique().tolist()
    pats=set(seqs)|{rc(s) for s in seqs}
    with tempfile.NamedTemporaryFile("w",suffix=".txt",delete=False,dir=PG) as pf:
        pf.write("\n".join(pats)+"\n"); pfn=pf.name
    res=subprocess.run(["grep","-F","-o","-f",pfn,ins],capture_output=True,text=True)
    os.unlink(pfn)
    found=set(res.stdout.split())
    d["in_priv_ins"]=d.sequence.map(lambda s: s in found or rc(s) in found)
    for klass,g in d.groupby("klass"):
        rate.append(dict(strain=X,klass=klass,n=len(g),in_priv_ins=int(g.in_priv_ins.sum()),
                         pct=round(100*g.in_priv_ins.mean(),3)))
    # TE-driven = strain-private-locus piRNA, embedded in an X-private insertion, AND TE-annotated (6a)
    tef=f"{S4}/{X}.private_TE_percandidate.csv.gz"
    sp=d[d.klass=="unique: strain-private locus"].copy()
    in_ins=int(sp.in_priv_ins.sum())
    te_drv=te_total=0; fam={}
    if os.path.exists(tef):
        te=pd.read_csv(tef)
        sp=sp.merge(te,on="id",how="left")
        hit=sp[sp.in_priv_ins & sp.classfam.notna()]
        te_drv=len(hit); te_total=int(sp.classfam.notna().sum())
        fam=hit.classfam.value_counts().head(8).to_dict()
    tedrv.append(dict(strain=X,strain_private=len(sp),in_private_insertion=in_ins,
                      TE_annotated=te_total,TE_driven=te_drv,top_families=fam))
    print(f"[{X}] strain-private={len(sp):,} | in X-private insertion={in_ins:,} | TE-driven (in priv-ins & TE)={te_drv:,}")
    print(f"        top TE-driven families: {fam}")

rr=pd.DataFrame(rate); rr.to_csv(f"{PG}/pirna_in_private_insertion_byclass.csv",index=False)
pd.DataFrame(tedrv).to_csv(f"{PG}/TE_driven_summary.csv",index=False)
print("\n=== % of candidates embedded in an X-PRIVATE insertion, by Step-4 class (enrichment/null) ===")
print(rr.pivot(index="klass",columns="strain",values="pct").to_string())
print("\nwrote pirna_in_private_insertion_byclass.csv + TE_driven_summary.csv")
