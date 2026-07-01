import pandas as pd, numpy as np
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
clean=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint"]); clean["L"]=clean.sequence.str.len()
WIN={"16.5dpc":[27],"12.5dpp":[27,30],"20.5dpp":[30]}; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}; rows=[]
for tp in ["16.5dpc","12.5dpp","20.5dpp"]:
    for st in clean.strain.unique():
        sub=clean[(clean.timepoint==tp)&(clean.strain==st)]; SETS={L:set(sub[sub.L==L].sequence) for L in range(23,35)}
        for L in WIN[tp]:
            sset=SETS[L]
            if len(sset)<30: continue
            pre={k:{y[:L] for y in SETS[L+k]} for k in (1,2,3)}
            rag=sum(1 for x in sset if any(x[:L-k] in SETS[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3)))
            rows.append(dict(strain=st,tp=TPN[tp],L=L,n=len(sset),rag_pct=100*rag/len(sset),wild=st in WILD))
df=pd.DataFrame(rows); df.to_csv("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/20_pirna_3prime_length_heterogeneity/data/SourceData_perstrain_ragged.csv",index=False)
print("=== per-strain ragged-3' %, by tp x window (does the finding hold across strains?) ===")
print(df.groupby(["tp","L"]).rag_pct.agg(n_strains="size",mean="mean",std="std",min="min",max="max").round(1))
print("\n=== developmental increase per strain (E16.5-27 -> P20.5-30), strains with both ===")
piv=df.pivot_table(index="strain",columns=["tp","L"],values="rag_pct")
print("strains where P20.5-30nt > E16.5-27nt raggedness:",int((piv[("P20.5",30)]>piv[("E16.5",27)]).sum()),"of",piv[[("E16.5",27),("P20.5",30)]].dropna().shape[0])
print("\nwild vs classical mean ragged% (P20.5 30nt):")
print(df[(df.tp=="P20.5")&(df.L==30)].groupby("wild").rag_pct.agg(['mean','size']).round(1))
