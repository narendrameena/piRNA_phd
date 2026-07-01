#!/usr/bin/env python3
"""TEST FIGURE — justifies the SNP-variant mismatch cutoff (≤3 mm = allelic variant of a piRNA expressed
at the homologous locus elsewhere -> NOT novel). (A) overall mismatch distribution of the 217,559 SNP-variants:
monotonic 1>2>3 = genuine allelic SNPs (a chance 29-mer match would INCREASE with mm). (B) per-strain mismatch
composition (canonical order): the most divergent (wild-derived) strains carry the largest 2-3 mm fraction
(SPRET 5.0%/CAST 4.0% at 3 mm vs classical <2%) -> a 2-mm cutoff would misclassify their allelic variants as
novel. Source: unique16/snp_variant_refinement.csv (mm = mismatches to the conserved homolog)."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from strain_order import STRAIN_ORDER, WILD
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
d=pd.read_csv(f"{U}/unique16/snp_variant_refinement.csv")
MM=[1,2,3]; MMC={1:"#bdbdbd",2:"#fdae61",3:"#d73027"}   # 3-mm slice highlighted (the one a 2-mm cutoff would drop)
order=[s for s in STRAIN_ORDER if s in set(d.home)]
fig,(axA,axB)=plt.subplots(1,2,figsize=(11,4.4),dpi=300,gridspec_kw={"width_ratios":[0.8,1.5]})

# A — overall mm distribution
oc=d.mm.value_counts().reindex(MM).fillna(0).astype(int); tot=oc.sum()
axA.bar([str(m) for m in MM],oc.values,width=0.66,color=[MMC[m] for m in MM],edgecolor="white",linewidth=0.5,zorder=3)
for i,m in enumerate(MM): axA.text(i,oc[m]+tot*0.012,f"{oc[m]:,}\n{100*oc[m]/tot:.1f}%",ha="center",va="bottom",fontsize=8,fontweight="bold")
axA.set_xlabel("mismatches to the conserved homolog (mm)",fontsize=8.5); axA.set_ylim(0,oc.max()*1.2)
axA.set_ylabel("SNP-variant piRNAs (count)",fontsize=9); axA.set_title("A  mismatch distribution",fontsize=9,fontweight="bold",loc="left")
axA.spines[["top","right"]].set_visible(False)

# B — per-strain mm composition (proportion), canonical order
comp=d.groupby(["home","mm"]).size().unstack().reindex(order)[MM].fillna(0)
prop=comp.div(comp.sum(1),axis=0); x=np.arange(len(order)); bottom=np.zeros(len(order))
for m in MM:
    axB.bar(x,prop[m].values,0.8,bottom=bottom,color=MMC[m],edgecolor="white",linewidth=0.3,label=f"{m} mm",zorder=3); bottom+=prop[m].values
axB.set_xticks(x); axB.set_xticklabels([s.replace("_","/") for s in order],rotation=45,ha="right",fontsize=7.5)
for lab,s in zip(axB.get_xticklabels(),order): lab.set_color("#C0392B" if s in WILD else "#222")
axB.set_ylim(0,1); axB.set_ylabel("mismatch composition (proportion)",fontsize=9)
axB.set_title("B  per-strain mismatch composition (wild-derived labels red)",fontsize=9,fontweight="bold",loc="left")
axB.legend(title="cutoff ≤3 mm",fontsize=8,title_fontsize=8,frameon=False,loc="upper left",bbox_to_anchor=(1.005,1.0),ncol=1)
axB.spines[["top","right"]].set_visible(False)

fig.suptitle("Choosing the SNP-variant mismatch cutoff (≤3 mm = allelic variant, NOT a novel strain-private piRNA)",fontsize=10.5,fontweight="bold",y=1.0)
fig.text(0.5,-0.04,"A 29-nt piRNA matching a conserved piRNA (expressed at the homologous locus in another strain) at 1–3 mm is an allelic SNP variant, not novel. "
  "The monotonic 1>2>3 distribution confirms genuine allelic variation (chance matches would rise with mm). The most divergent strains (SPRET 5.0%, CAST 4.0% at 3 mm) need the ≤3 cutoff — "
  "a 2-mm cutoff would mis-call their allelic variants as strain-private.",ha="center",fontsize=6.6,color="#555")
fig.tight_layout(rect=[0,0.02,1,0.99])
base=f"{U}/Fig_snp_variant_cutoff"
for e in ("pdf","svg","png"): fig.savefig(f"{base}.{e}",bbox_inches="tight")
sd=comp.assign(total=comp.sum(1)); sd.columns=[f"mm{c}" if isinstance(c,int) else c for c in sd.columns]
sd.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data/SourceData_Fig_snp_variant_cutoff.csv")
print(oc.to_string()); print("wrote",base)
