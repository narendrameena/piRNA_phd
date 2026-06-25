import gzip, pandas as pd, numpy as np
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; ED=f"{U}/edger16"
clean=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","timepoint"]); clean["L"]=clean.sequence.str.len()
WIN={"16.5dpc":[27],"12.5dpp":[27,30],"20.5dpp":[30]}; res=[]
for tp in ["16.5dpc","12.5dpp","20.5dpp"]:
    sub=clean[clean.timepoint==tp]; SETS={L:set(sub[sub.L==L].sequence) for L in range(23,35)}; cleanset=set(sub.sequence); abund={}
    with gzip.open(f"{ED}/{tp}.seqs.txt.gz","rt") as sf, gzip.open(f"{ED}/{tp}.counts.tsv.gz","rt") as cf:
        cf.readline()
        for sl,cl in zip(sf,cf):
            s=sl[:-1]
            if s in cleanset: abund[s]=sum(int(x) for x in cl.rstrip("\n").split("\t"))
    for L in WIN[tp]:
        sset=SETS[L]; pre={k:{y[:L] for y in SETS[L+k]} for k in (1,2,3)}
        def israg(x): return any(x[:L-k] in SETS[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3))
        rag=[x for x in sset if israg(x)]; tot=sum(abund.get(x,0) for x in sset); ragab=sum(abund.get(x,0) for x in rag)
        res.append(dict(tp=tp,L=L,n=len(sset),rag_pct_seq=round(100*len(rag)/len(sset),1),rag_pct_reads=round(100*ragab/tot,1) if tot else 0,
                        med_ab_rag=float(np.median([abund.get(x,0) for x in rag])) if rag else 0,
                        med_ab_norag=float(np.median([abund.get(x,0) for x in sset if not israg(x)])) or 0))
        print(f"{tp} {L}nt: seq-ragged={100*len(rag)/len(sset):.0f}%  READ-ragged={100*ragab/tot:.0f}%  median-abund ragged {np.median([abund.get(x,0) for x in rag]):.0f} vs non-ragged {np.median([abund.get(x,0) for x in sset if not israg(x)]):.0f}",flush=True)
pd.DataFrame(res).to_csv("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/20_pirna_3prime_length_heterogeneity/data/SourceData_expression_ragged.csv",index=False); print("done")
