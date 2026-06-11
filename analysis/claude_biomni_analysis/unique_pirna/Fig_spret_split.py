#!/usr/bin/env python3
"""SPRET SNP-variant-vs-locus split, per timepoint. Classifies each SPRET presence/absence
candidate by its min mismatch to a piRNA EXPRESSED in C57BL_6NJ/CAST:
  0mm   = exact sequence expressed-low in others (below RPM threshold) -> NOT truly unique
  1-3mm = SNP-variant of a conserved expressed piRNA                   -> NOT truly unique
  no hit= genuinely SPRET-novel sequence (no <=3mm expressed match)    -> REAL
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
# sN -> sequence
sN2seq={}
with open(f"{U}/SPRET_candidates.fasta") as fh:
    sid=None
    for ln in fh:
        if ln[0]==">": sid=ln[1:].strip()
        else: sN2seq[sid]=ln.strip()
mm=pd.read_csv(f"{U}/SPRET_candidate_minmm_to_others_expressed.tsv",sep="\t",header=None,names=["sN","mm"])
seq2mm={sN2seq[s]:m for s,m in zip(mm.sN,mm.mm)}
cand=pd.read_csv(f"{U}/strain_specific_presenceAbsence_candidates.csv.gz")
sp=cand[cand.strain=="SPRET_EiJ"].copy()
def cls(seq):
    m=seq2mm.get(seq)
    if m is None: return "genuinely novel (no <=3mm)"
    return "expressed-low in others (0mm)" if m==0 else "SNP-variant (1-3mm)"
sp["klass"]=sp.sequence.map(cls)
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
sp["tp"]=sp.timepoint.map(TPMAP)
ORDER=["expressed-low in others (0mm)","SNP-variant (1-3mm)","genuinely novel (no <=3mm)"]
COLS={"expressed-low in others (0mm)":"#bdbdbd","SNP-variant (1-3mm)":"#E69F00","genuinely novel (no <=3mm)":"#C0392B"}
tab=sp.groupby(["tp","klass"]).size().unstack(fill_value=0).reindex(TPO)[ORDER]
fig,ax=plt.subplots(figsize=(6.2,4.2),dpi=300)
plt.rcParams.update({"font.family":"Liberation Sans"})
x=np.arange(3); bottom=np.zeros(3)
for k in ORDER:
    ax.bar(x,tab[k].values,bottom=bottom,width=0.62,color=COLS[k],edgecolor="white",linewidth=0.5,label=k,zorder=3)
    bottom+=tab[k].values
for xi,t in zip(x,TPO):
    nov=tab.loc[t,"genuinely novel (no <=3mm)"]; tot=tab.loc[t].sum()
    ax.text(xi,tot+500,f"novel: {nov:,}\n({100*nov/tot:.1f}%)",ha="center",va="bottom",fontsize=6.5,color="#C0392B",fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(TPO,fontsize=9); ax.set_ylabel("SPRET presence/absence candidates",fontsize=9)
ax.legend(fontsize=6.8,frameon=False,loc="upper left",bbox_to_anchor=(0.0,-0.10),ncol=1)
ax.set_title("SPRET_EiJ: naive 'unique' candidates split by expression in other strains",fontsize=9.5,fontweight="bold")
fig.text(0.5,-0.20,"~99% are expressed (exactly or as a SNP-variant) in C57BL_6NJ/CAST -> not truly unique. "
  "Only the red sliver is genuinely SPRET-novel (no <=3mm expressed match). Mismatch cutoff <=3 (sensitivity: <=1mm gives more novel).",
  ha="center",fontsize=5.3,color="#666")
fig.tight_layout()
out=f"{U}/Fig_spret_split"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
tab.to_csv(f"{U}/SourceData_spret_split.csv")
print(tab.to_string()); print("\nTOTAL by class:\n", sp.klass.value_counts().to_string())
print("\nSensitivity (genuinely-novel = no expressed match within cutoff):")
for c in (1,2,3):
    nov=sum(1 for seq in sp.sequence.unique() if seq2mm.get(seq,99)> c-1 and seq2mm.get(seq,99)>=c)  # placeholder
print("  <=1mm cutoff novel (no 0/1mm match):", (~sp.sequence.map(lambda s: seq2mm.get(s,99)<=1)).sum(),
      "| <=3mm novel:", (~sp.sequence.map(lambda s: seq2mm.get(s,99)<=3)).sum())
