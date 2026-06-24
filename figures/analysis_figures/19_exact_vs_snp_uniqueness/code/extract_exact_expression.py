#!/usr/bin/env python3
"""THEME 19 — extract expression of the EXACT-sequence genuinely-unique piRNAs (15,118 entries incl. the 4,394
SNP-alleles) for the heatmap + PCA. Streams the edger16 count matrices ONCE; caches per-rep CPM (seq × 144 =
3 tp × 48 libraries, cols 'tp|sample', linear CPM by libsize_window). Heatmap derives mean-per-strain×tp; PCA
uses per-rep per tp. Output: data/exact_cpm_perrep.csv.gz."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import gzip, numpy as np, pandas as pd
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; ED=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/edger16"
T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
TPS=["16.5dpc","12.5dpp","20.5dpp"]
GU=["unique: conserved-but-silent","unique: strain-private locus","unique: stage-shifted (heterochronic)"]
d=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz"); targets=set(d[d.klass_exact.isin(GU)].sequence)
print(f"exact-unique distinct sequences: {len(targets):,}",flush=True)
frames=[]
for tp in TPS:
    with gzip.open(f"{ED}/{tp}.counts.tsv.gz","rt") as cf: header=cf.readline().rstrip("\n").split("\t")
    smp=pd.read_csv(f"{ED}/{tp}.samples.tsv",sep="\t").set_index("sample")
    lib=np.array([smp.loc[h,"libsize_window"] for h in header],float)
    rows={}
    with gzip.open(f"{ED}/{tp}.seqs.txt.gz","rt") as sf, gzip.open(f"{ED}/{tp}.counts.tsv.gz","rt") as cf:
        cf.readline()
        for sl,cl in zip(sf,cf):
            if sl[:-1] in targets: rows[sl[:-1]]=np.array(cl.rstrip("\n").split("\t"),float)
    M=np.vstack(list(rows.values())); cpm=M/lib*1e6
    frames.append(pd.DataFrame(cpm,index=list(rows.keys()),columns=[f"{tp}|{h}" for h in header]))
    print(f"{tp}: {len(rows):,} seqs",flush=True)
E=pd.concat(frames,axis=1).reindex(sorted(targets)).fillna(0.0)
E.to_csv(f"{T}/data/exact_cpm_perrep.csv.gz")
print(f"wrote exact_cpm_perrep.csv.gz {E.shape}")
