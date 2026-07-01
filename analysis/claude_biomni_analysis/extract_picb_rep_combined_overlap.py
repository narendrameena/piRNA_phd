#!/usr/bin/env python3
"""
Replicate-vs-combined PICB cluster OVERLAP (does a single replicate find the SAME
clusters as the combined run?). For each strain x timepoint, genomic-interval overlap
(same chromosome + same strand, ANY overlap >=1 bp) between the 3 per-replicate cluster
sets and the combined set:
  - combined clusters supported by 3 / 2 / 1 / 0 of the single replicates
  - fraction of single-replicate clusters retained in the combined run (recovery)
  - replicate-to-replicate reproducibility
ANY-overlap is used (not reciprocal) because the combined run can consolidate adjacent
clusters; the question is "is there a cluster at this locus in the replicate?".
Writes per-group metrics CSV.
"""
import glob, os, warnings
import numpy as np, pandas as pd
warnings.simplefilter("ignore")

ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
ISV=["129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","CAST_EiJ","CBA_J","DBA_2J",
     "FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ","C57BL_6NJ"]
TPS=["16.5dpc","12.5dpp","20.5dpp"]

def clusters(path):
    if not os.path.exists(path): return None
    d=pd.read_excel(path, sheet_name="clusters")
    return d[["seqnames","start","end","strand"]].copy()

def overlap_flags(q, r):
    """boolean per q-interval: overlaps >=1 r-interval on same (seqnames,strand)."""
    out=np.zeros(len(q),bool)
    if q is None or r is None or len(q)==0 or len(r)==0: return out
    q=q.reset_index(drop=True)
    for (chrom,strand),rg in r.groupby(["seqnames","strand"]):
        qm=((q.seqnames==chrom)&(q.strand==strand)).values
        if not qm.any(): continue
        st=rg.start.values; en=rg.end.values
        o=np.argsort(st); st=st[o]; pmax=np.maximum.accumulate(en[o])
        qs=q.loc[qm,"start"].values; qe=q.loc[qm,"end"].values
        k=np.searchsorted(st,qe,side="right")
        res=np.zeros(len(qs),bool); ok=k>0
        res[ok]=pmax[k[ok]-1]>=qs[ok]
        out[np.where(qm)[0]]=res
    return out

rows=[]
for strain in ISV:
    for tp in TPS:
        c=clusters(f"{ROOT}/results/picb_result_combined/{strain}/{strain}-{tp}.combined.xlsx")
        reps=[clusters(f"{ROOT}/results/picb_result/{strain}/{strain}-{tp}.{i}.xlsx") for i in (1,2,3)]
        reps=[r for r in reps if r is not None]
        if c is None or len(reps)<2:
            print(f"  skip {strain}-{tp} (missing files: combined={c is not None}, reps={len(reps)})"); continue
        nrep=len(reps)
        support=np.zeros(len(c),int)
        for r in reps: support+=overlap_flags(c,r).astype(int)
        rec=[overlap_flags(r,c).mean() for r in reps]                 # rep clusters found in combined
        rr=[overlap_flags(reps[i],reps[j]).mean() for i in range(nrep) for j in range(nrep) if i!=j]
        rows.append(dict(strain=strain,timepoint=tp,n_reps=nrep,
            n_combined=len(c), n_rep_mean=float(np.mean([len(r) for r in reps])),
            comb_support3=int((support==3).sum()),comb_support2=int((support==2).sum()),
            comb_support1=int((support==1).sum()),comb_support0=int((support==0).sum()),
            frac_comb_in_anyrep=float((support>=1).mean()),
            frac_comb_in_allrep=float((support==nrep).mean()),
            frac_comb_unique=float((support==0).mean()),
            frac_rep_in_combined=float(np.mean(rec)),
            rep_rep_reproducibility=float(np.mean(rr))))
        print(f"  {strain}-{tp}: nC={len(c)} | comb in any rep={100*(support>=1).mean():.1f}% "
              f"in all reps={100*(support==nrep).mean():.1f}% unique={100*(support==0).mean():.1f}% | "
              f"rep->comb={100*np.mean(rec):.1f}% rep-rep={100*np.mean(rr):.1f}%", flush=True)

df=pd.DataFrame(rows)
out=f"{ROOT}/analysis/claude_biomni_analysis/source_data/SourceData_PICB_rep_combined_overlap.csv"
df.to_csv(out,index=False)
print(f"\nWROTE {out} ({len(df)} groups)")
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
df["tp"]=df.timepoint.map(TPMAP)
print("\n=== mean by timepoint (%) ===")
print((df.groupby("tp")[["frac_comb_in_anyrep","frac_comb_in_allrep","frac_comb_unique",
       "frac_rep_in_combined","rep_rep_reproducibility"]].mean()*100).reindex(["E16.5","P12.5","P20.5"]).round(1).to_string())
