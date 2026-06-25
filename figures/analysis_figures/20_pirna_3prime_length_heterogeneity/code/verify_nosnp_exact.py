"""Does the ragged-3' signal depend on SNP-variant piRNAs? Re-run window raggedness on (a) ALL clean_2read,
(b) EXCLUDING klass=='SNP-variant', (c) EXCLUDING SNP-variant + low-quality (i.e. only exact/conserved + strain-private).
Exact prefix matching throughout (set membership = zero mismatch)."""
import pandas as pd, numpy as np
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","timepoint","klass"]); d["L"]=d.sequence.str.len()
print("klass counts:",dict(d.klass.value_counts()))
WIN={"16.5dpc":[27],"12.5dpp":[27,30],"20.5dpp":[30]}; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
def rag_pct(sub,tp,L):
    S={k:set(sub[(sub.timepoint==tp)&(sub.L==k)].sequence) for k in range(23,35)}; sset=S[L]
    if not sset: return 0,0
    pre={k:{y[:L] for y in S[L+k]} for k in (1,2,3)}
    r=sum(1 for x in sset if any(x[:L-k] in S[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3)))
    return len(sset),100*r/len(sset)
subsets={"ALL":d,"no SNP-variant":d[d.klass!="SNP-variant"],"no SNP-var + no low-qual":d[~d.klass.isin(["SNP-variant","low-quality"])]}
print(f"\n{'tp/win':<12}"+"".join(f"{name:>26}" for name in subsets))
for tp in ["16.5dpc","12.5dpp","20.5dpp"]:
    for L in WIN[tp]:
        cells=[]
        for name,sub in subsets.items():
            n,p=rag_pct(sub,tp,L); cells.append(f"{p:.0f}% (n={n:,})")
        print(f"{TPN[tp]+' '+str(L)+'nt':<12}"+"".join(f"{c:>26}" for c in cells))
