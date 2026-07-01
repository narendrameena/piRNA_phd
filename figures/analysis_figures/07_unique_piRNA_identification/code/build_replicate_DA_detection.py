#!/usr/bin/env python3
"""Per-REPLICATE detection of the strain-specific DA piRNA set (for Fig_strain_specific_DA16 error bars).

The strain-specific count is computed ONCE per strain x timepoint by edgeR on all replicate libraries
together (edgeR cannot run on a single replicate). To show biological-replicate variability we recompute,
for each replicate library, how many of that strain's strain-specific piRNAs it individually detects
(>=1 read) — the per-replicate "recovery" of the set. Mean ~ set size (each member is present in >=2/3 reps);
SD = replicate consistency. Same philosophy as theme-03's replicate error bars.

Streams the big count matrix (6.7M seqs x 48 libs) but only parses the ~150k strain-specific rows
(row indices found via the line-aligned seqs.txt.gz). Output: edger16/replicate_detected_strain_specific.csv
"""
import gzip, os
import pandas as pd, numpy as np

DIR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}

rows=[]
for f,tp in TPMAP.items():
    samp=pd.read_csv(f"{DIR}/{f}.samples.tsv",sep="\t")              # row order == count-matrix column order
    # verify header (library names) matches samples order exactly
    with gzip.open(f"{DIR}/{f}.counts.tsv.gz","rt") as cf:
        header=cf.readline().rstrip("\n").split("\t")
    assert header==samp["sample"].tolist(), f"{f}: counts header != samples order"
    strain_cols={}                                                  # strain -> [(col_idx, rep), ...]
    for ci,(st,rp) in enumerate(zip(samp["strain"],samp["rep"])):
        strain_cols.setdefault(st,[]).append((ci,int(rp)))
    da=pd.read_csv(f"{DIR}/{f}.strain_specific_DA_2read.csv.gz")   # ADOPTED ≥2-read absence set
    seq2strain=dict(zip(da.sequence,da.strain))
    set_size=da.strain.value_counts().to_dict()
    # target row indices (line-aligned with counts data rows)
    target={}
    with gzip.open(f"{DIR}/{f}.seqs.txt.gz","rt") as sf:
        for i,line in enumerate(sf):
            st=seq2strain.get(line.rstrip("\n"))
            if st is not None: target[i]=st
    det={st:{rp:0 for (_,rp) in cols} for st,cols in strain_cols.items()}
    with gzip.open(f"{DIR}/{f}.counts.tsv.gz","rt") as cf:
        cf.readline()                                               # skip header
        for i,line in enumerate(cf):                                # i == seqs row index
            st=target.get(i)
            if st is None: continue
            vals=line.split("\t")
            for ci,rp in strain_cols[st]:
                if int(vals[ci])>=1: det[st][rp]+=1
    for st,reps in det.items():
        for rp,d in reps.items():
            rows.append(dict(strain=st,timepoint=tp,rep=rp,detected=d,set_size=int(set_size.get(st,0))))
    print(f"[{f}] done: {len(da)} strain-specific seqs over {da.strain.nunique()} strains")

out=pd.DataFrame(rows)
out.to_csv(f"{DIR}/replicate_detected_strain_specific_2read.csv",index=False)
print(f"\nwrote {DIR}/replicate_detected_strain_specific.csv ({len(out)} rows)")

# summary: set size vs per-rep mean/SD (decide bar/error representation)
print("\nstrain x tp : set_size  rep-mean  rep-SD  (per-rep detected of the strain-specific set)")
for tp in TPMAP.values():
    for st in out.strain.unique():
        sub=out[(out.strain==st)&(out.timepoint==tp)]
        if not len(sub): continue
        ss=sub.set_size.iloc[0]
        print(f"  {st:13s} {tp}: {ss:7d}  {sub.detected.mean():8.0f}  {sub.detected.std(ddof=1):7.0f}  "
              f"(mean/size={sub.detected.mean()/ss:.2f})" if ss else "")
