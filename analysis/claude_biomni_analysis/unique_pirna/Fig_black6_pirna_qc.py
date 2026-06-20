#!/usr/bin/env python3
"""Definitive piRNA-enrichment QC for the EXTERNAL black6 (C57BL/6) sRNA vs our C57BL_6NJ. For each sample,
subsample 3M reads, adapter-trim (black6 = TruSeq small-RNA -a TCGTATGCCGTCTTCTGCTTGT -g GTTCAGAGTTCTACAGTCCGACGATC,
its real prep; 6NJ = project -a TGGAATTCTCGGGTGCCAAGG), keep 18-45 nt inserts, then the LENGTH distribution
(piRNA 24-32 nt peak vs miRNA ~22 nt) and 1U bias — the signatures of oxidised, piRNA-enriched libraries.
black6 timepoints come from different SRA studies. E16.5 = ALL 6 controls (3 TEX15/GSE150350 SRR11774396-398 + 3 SPOCD1/GSE131377
SRR9077640-642, both O'Carroll lab; Illumina/NEBnext adapter, prior-pipeline-trimmed); P12.5/P20.5 = SRR7720xx (TruSeq, trimmed here)."""
import warnings; warnings.filterwarnings("ignore")
import subprocess, os, tempfile, gzip
from collections import Counter
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; SR=f"{ROOT}/resources/mice16_data/srna"; PG=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
B6B=f"{ROOT}/resources/black6"   # E16.5 black6 = ALL 6 controls from 2 O'Carroll-lab studies, prior-pipeline-trimmed (correct Illumina/NEBnext adapter; the project TruSeq adapter does NOT match these libraries)
B6_E165=[f"{B6B}/E16_5/sRNA/trimmed/SRR1177439{n}_trimmed.fastq.gz" for n in (6,7,8)]+[f"{B6B}/GSE131377/trimmed/SRR907764{n}_trimmed.fastq.gz" for n in (0,1,2)]
CUT="/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/cutadapt"
TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
B6A=["-a","TCGTATGCCGTCTTCTGCTTGT","-g","GTTCAGAGTTCTACAGTCCGACGATC"]; NJA=["-a","TGGAATTCTCGGGTGCCAAGG"]
UNIQ={"16.5dpc":["1","2","3"],"12.5dpp":["1","2"],"20.5dpp":["1","2"]}   # black6 unique reps (dropping the byte-identical duplicates)
N=3_000_000; LR=list(range(18,46))
def qc(fq,adapt):
    if not os.path.exists(fq): return None
    sub=tempfile.NamedTemporaryFile(suffix=".fq",delete=False,dir=PG); trim=tempfile.NamedTemporaryFile(suffix=".fq",delete=False,dir=PG)
    subprocess.run(f"zcat {fq} | head -{N*4} > {sub.name}",shell=True)
    subprocess.run([CUT,*adapt,"--minimum-length","18","--maximum-length","45","--discard-untrimmed","-o",trim.name,sub.name],capture_output=True,text=True)
    lens=Counter(); first=Counter(); tot=0
    with open(trim.name) as fh:
        for i,l in enumerate(fh):
            if i%4==1:
                s=l.strip();
                if s: lens[len(s)]+=1; first[s[0]]+=1; tot+=1
    os.unlink(sub.name); os.unlink(trim.name)
    return lens,first,tot
def lens_1u(fq,n=N):   # insert-length + 5'nt from an ALREADY adapter-trimmed fastq (E16.5 black6 controls were trimmed by the prior pipeline with the correct Illumina/NEBnext adapter)
    if not os.path.exists(fq): return None
    lens=Counter(); first=Counter(); tot=0
    with gzip.open(fq,"rt") as fh:
        for i,l in enumerate(fh):
            if i>=n*4: break
            if i%4==1:
                s=l.strip()
                if s and 18<=len(s)<=45: lens[len(s)]+=1; first[s[0]]+=1; tot+=1
    return lens,first,tot
res={}  # (strain,tp) -> (lenpct over LR, u1, n, passfrac)
for lab,tp in TPS:
    for strain,adapt,reps in [("black6",B6A,UNIQ[tp]),("C57BL_6NJ",NJA,["1"])]:
        ls=Counter(); fs=Counter(); tot=0
        if strain=="black6" and tp=="16.5dpc":   # ALL 6 E16.5 controls (3 TEX15/GSE150350 + 3 SPOCD1/GSE131377), already adapter-trimmed by the prior pipeline
            for fq in B6_E165:
                o=lens_1u(fq)
                if o: ls+=o[0]; fs+=o[1]; tot+=o[2]
        else:
            for r in reps:
                fq=f"{SR}/{('C57BL_6' if strain=='black6' else 'C57BL_6NJ')}-{tp}.{r}.1s.fastq.gz"
                o=qc(fq,adapt)
                if o: ls+=o[0]; fs+=o[1]; tot+=o[2]
        if tot==0: continue
        lp=np.array([100*ls.get(L,0)/tot for L in LR]); u1=100*fs.get("T",0)/sum(fs.values())
        res[(strain,lab)]=(lp,u1,tot)
        pk=LR[int(np.argmax(lp))]; pir=sum(100*ls.get(L,0)/tot for L in range(24,33))
        print(f"{strain} {lab}: modal len {pk} nt | 24-32nt = {pir:.0f}% | 1U = {u1:.0f}% | n(trimmed,sampled)={tot:,}")
# ---- figure ----
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,ax=plt.subplots(1,3,figsize=(14,4.6),dpi=300); xs=np.array(LR); SC={"black6":"#D55E00","C57BL_6NJ":"#0072B2"}
for j,(lab,tp) in enumerate(TPS):
    a=ax[j]; a.axvspan(23.5,32.5,color="#0072B2",alpha=0.06,zorder=0); a.axvline(22,ls=":",lw=0.9,color="#888")
    a.text(22,a.get_ylim()[1] if False else 1,"miRNA\n~22",fontsize=6,color="#888",ha="center",va="bottom")
    for strain in ("black6","C57BL_6NJ"):
        if (strain,lab) in res:
            lp,u1,tot=res[(strain,lab)]; a.plot(xs,lp,color=SC[strain],lw=1.9,marker="o",ms=2.6,label=f"{strain} (1U {u1:.0f}%)")
    a.set_title(f"{lab}",fontsize=9.5,fontweight="bold"); a.set_xlabel("insert length (nt)",fontsize=8.5)
    if j==0: a.set_ylabel("% of trimmed reads",fontsize=8.5)
    a.set_xticks(xs[::3]); a.tick_params(labelsize=6.5); a.legend(fontsize=6.6,frameon=False); a.spines[['top','right']].set_visible(False)
    a.text(28,a.get_ylim()[1]*0.96,"piRNA 24-32",fontsize=6,color="#0072B2",ha="center",va="top")
fig.suptitle("Definitive piRNA-enrichment QC — external black6 (C57BL/6) vs C57BL_6NJ (project): insert-length + 1U",fontsize=10.2,fontweight="bold",y=1.02)
fig.text(0.5,-0.04,"piRNA-enriched libraries peak at 24-32 nt with 5'-U (1U). 3M reads/sample. E16.5 black6 = ALL 6 controls pooled (3 TEX15/GSE150350 + 3 SPOCD1/GSE131377), prior-pipeline-trimmed with the Illumina/NEBnext adapter (the project's TruSeq adapter does NOT match these libraries); P12.5/P20.5 black6 = unique reps (TruSeq).",ha="center",fontsize=5.8,color="#555")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_black6_pirna_qc.{e}",bbox_inches="tight")
import pandas as pd
pd.DataFrame({"length":LR,**{f"{s}_{l}_pct":res[(s,l)][0] for (s,l) in res}}).to_csv(f"{PG}/SourceData_black6_pirna_qc.csv",index=False)
print("wrote Fig_black6_pirna_qc.{png,pdf,svg} + source data")
