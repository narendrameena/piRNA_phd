#!/usr/bin/env python3
"""Per-SAMPLE black6 piRNA-enrichment QC: each individual sample = its own line (insert-length + 1U), labelled by name.
E16.5 black6 = ALL 6 controls (3 TEX15/GSE150350 SRR11774396-398 + 3 SPOCD1/GSE131377 SRR9077640-642; Illumina/NEBnext,
prior-pipeline-trimmed). P12.5/P20.5 black6 = unique reps (TruSeq, trimmed here). C57BL_6NJ = project reps 1-3 (TruSeq)."""
import warnings; warnings.filterwarnings("ignore")
import subprocess, os, tempfile, gzip
from collections import Counter
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; SR=f"{ROOT}/resources/mice16_data/srna"; PG=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"; B6B=f"{ROOT}/resources/black6"
CUT="/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/cutadapt"
B6A=["-a","TCGTATGCCGTCTTCTGCTTGT","-g","GTTCAGAGTTCTACAGTCCGACGATC"]; NJA=["-a","TGGAATTCTCGGGTGCCAAGG"]
B6_E165=[(f"TEX15·{s}", f"{B6B}/E16_5/sRNA/trimmed/SRR1177439{s}_trimmed.fastq.gz") for s in (6,7,8)]+[(f"SPOCD1·{s}", f"{B6B}/GSE131377/trimmed/SRR907764{s}_trimmed.fastq.gz") for s in (0,1,2)]
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
                s=l.strip()
                if s: lens[len(s)]+=1; first[s[0]]+=1; tot+=1
    os.unlink(sub.name); os.unlink(trim.name)
    return lens,first,tot
def lens_1u(fq,n=N):
    if not os.path.exists(fq): return None
    lens=Counter(); first=Counter(); tot=0
    with gzip.open(fq,"rt") as fh:
        for i,l in enumerate(fh):
            if i>=n*4: break
            if i%4==1:
                s=l.strip()
                if s and 18<=len(s)<=45: lens[len(s)]+=1; first[s[0]]+=1; tot+=1
    return lens,first,tot
def curve(o):
    lens,first,tot=o; lp=np.array([100*lens.get(L,0)/tot for L in LR]); u1=100*first.get("T",0)/(sum(first.values()) or 1); return lp,u1
# per-tp sample list: (group, name, loader)
def b6raw(tp,r): return ("black6", f"black6·r{r}", (lambda: qc(f"{SR}/C57BL_6-{tp}.{r}.1s.fastq.gz", B6A)))
def njraw(tp,r): return ("C57BL_6NJ", f"C57BL_6NJ·r{r}", (lambda: qc(f"{SR}/C57BL_6NJ-{tp}.{r}.1s.fastq.gz", NJA)))
SAMPLES={
  "16.5dpc": [("black6",f"black6 {nm}", (lambda f=f: lens_1u(f))) for nm,f in B6_E165] + [njraw("16.5dpc",r) for r in (1,2,3)],
  "12.5dpp": [b6raw("12.5dpp",r) for r in (1,2)] + [njraw("12.5dpp",r) for r in (1,2,3)],
  "20.5dpp": [b6raw("20.5dpp",r) for r in (1,2)] + [njraw("20.5dpp",r) for r in (1,2,3)],
}
TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,ax=plt.subplots(1,3,figsize=(15.5,5.2),dpi=300); xs=np.array(LR)
_src=[]   # source data: per-sample insert-length distribution + 1U%
for j,(lab,tp) in enumerate(TPS):
    a=ax[j]; a.axvspan(24.5,32.5,color="#0072B2",alpha=0.06,zorder=0); a.axvline(22,ls=":",lw=0.9,color="#888")
    samps=SAMPLES[tp]
    b6=[s for s in samps if s[0]=="black6"]; nj=[s for s in samps if s[0]=="C57BL_6NJ"]
    cols={"black6":cm.Oranges(np.linspace(0.45,0.92,max(len(b6),1))), "C57BL_6NJ":cm.Blues(np.linspace(0.45,0.92,max(len(nj),1)))}
    gi={"black6":0,"C57BL_6NJ":0}
    for grp,name,load in samps:
        o=load()
        if not o or o[2]==0: gi[grp]+=1; continue
        lp,u1=curve(o); c=cols[grp][gi[grp]]; gi[grp]+=1
        _src.append({"timepoint":lab,"group":grp,"sample":name,"oneU_pct":round(float(u1),2),**{f"len{L}nt_pct":round(float(lp[L-18]),4) for L in LR}})
        a.plot(xs,lp,color=c,lw=1.4,marker="o",ms=2.0,label=f"{name} (1U {u1:.0f}%)")
        print(f"{lab} {name}: 1U {u1:.0f}% | 25-32nt {sum(lp[L-18] for L in range(25,33)):.0f}%")
    a.set_title(lab,fontsize=10,fontweight="bold"); a.set_xlabel("insert length (nt)",fontsize=9)
    if j==0: a.set_ylabel("% of trimmed reads",fontsize=9)
    a.set_xticks(xs[::3]); a.tick_params(labelsize=7); a.legend(fontsize=5.6,frameon=False,loc="upper right"); a.spines[['top','right']].set_visible(False)
    a.text(28,a.get_ylim()[1]*0.97,"piRNA 25-32",fontsize=6,color="#0072B2",ha="center",va="top")
fig.suptitle("Per-sample piRNA-enrichment QC — black6 (C57BL/6) vs C57BL_6NJ: each sample, insert-length + 1U",fontsize=11,fontweight="bold",y=1.02)
fig.text(0.5,-0.03,"Each line = ONE sample. E16.5 black6 = all 6 controls (3 TEX15/GSE150350 + 3 SPOCD1/GSE131377; NEBnext, prior-pipeline-trimmed). P12.5/P20.5 black6 = unique reps; C57BL_6NJ = project reps (TruSeq, trimmed here). 3M reads/sample.",ha="center",fontsize=6.2,color="#555")
fig.tight_layout()
import pandas as _pd, os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/06_zamore_coverage/data/source_data"; _os.makedirs(_SD,exist_ok=True)
_pd.DataFrame(_src).to_csv(f"{_SD}/SourceData_Fig_black6_pirna_qc_persample.csv",index=False)   # per-sample insert-length distribution (% per nt, 18-45) + 1U%
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_black6_pirna_qc_persample.{e}",bbox_inches="tight")
print("wrote Fig_black6_pirna_qc_persample.{png,pdf,svg}")
