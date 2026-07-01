import pandas as pd, numpy as np
from scipy.stats import fisher_exact
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
clean=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint"])
clean["L"]=clean.sequence.str.len(); clean["cand_id"]=clean.strain+"|"+clean.timepoint+"|"+clean.sequence
oc=pd.read_csv(f"{U}/sense_antisense/SourceData_sense_antisense16_percand.csv.gz")
TE=dict(zip(oc.id,oc.family.notna())); ORI=dict(zip(oc.id,oc.orientation))
WIN={"16.5dpc":[27],"12.5dpp":[27,30],"20.5dpp":[30]}; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
print("=== does TE-derivation / antisense-orientation affect 3' raggedness? ===")
for tp in ["16.5dpc","12.5dpp","20.5dpp"]:
    sub=clean[clean.timepoint==tp]; SETS={L:set(sub[sub.L==L].sequence) for L in range(23,35)}
    for L in WIN[tp]:
        sset=SETS[L]; pre={k:{y[:L] for y in SETS[L+k]} for k in (1,2,3)}
        israg=lambda x: any(x[:L-k] in SETS[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3))
        subL=sub[sub.L==L]; seqte={}; seqori={}
        for r in subL.itertuples():
            seqte[r.sequence]=seqte.get(r.sequence,False) or TE.get(r.cand_id,False)
            o=ORI.get(r.cand_id)
            if isinstance(o,str): seqori.setdefault(r.sequence,[]).append(o)
        df=pd.DataFrame([(s,israg(s),seqte.get(s,False),(max(set(seqori[s]),key=seqori[s].count) if s in seqori else None)) for s in sset],columns=["seq","rag","te","ori"])
        rt=df.groupby("te").rag.mean()*100; p_te=fisher_exact(pd.crosstab(df.te,df.rag).values)[1]
        dto=df[df.ori.notna()];
        ro=dto.groupby("ori").rag.mean()*100 if len(dto)>10 else pd.Series(dtype=float)
        p_ori=fisher_exact(pd.crosstab(dto.ori,dto.rag).values)[1] if dto.ori.nunique()==2 else float("nan")
        print(f"{TPN[tp]} {L}nt: ragged by TE-derived {dict(rt.round(0))} (Fisher p={p_te:.1e}) | TE-overlapping fraction {100*df.te.mean():.0f}% | ragged by orientation {dict(ro.round(0))} (p={p_ori:.1e})")
