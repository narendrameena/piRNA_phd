#!/usr/bin/env python3
"""THEME 18 — edgeR vs DESeq2 benchmark for strain-specific piRNA differential abundance (data-based, BioMNI
3/3-verified design). One-vs-rest per strain on a 300k-feature representative sample, 3 tp x 5 permutations.
Decisive result: edgeR is ~2x anti-conservative and produces ~10-20x more NULL false-positives; DESeq2 is
well-calibrated (null p<0.05 ~= 0.05). -> DESeq2 adopted. Panels: A real-data calls, B null false-positives
(log), C null p<0.05 calibration vs target, D null p-value distribution (uniform expected).
Data: unique_pirna/bench_da/{tp}.{concordance,permnull,nullpv_hist}.csv"""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
B=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/bench_da"
T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
CE="#D55E00"; CD="#0072B2"   # edgeR vermillion, DESeq2 blue (colourblind-safe)
conc={tp:pd.read_csv(f"{B}/{tp}.concordance.csv") for tp in TPS}
perm={tp:pd.read_csv(f"{B}/{tp}.permnull.csv") for tp in TPS}
hist={tp:pd.read_csv(f"{B}/{tp}.nullpv_hist.csv") for tp in TPS}
x=np.arange(len(TPS)); bw=0.38
fig,axes=plt.subplots(2,2,figsize=(13.5,10),dpi=300); (axA,axB),(axC,axD)=axes
# A: real-data calls (overall = sum of per-strain calls; both methods, grouped)
eA=[conc[tp].edgeR.sum() for tp in TPS]; dA=[conc[tp].DESeq2.sum() for tp in TPS]
axA.bar(x-bw/2,eA,bw,color=CE,label="edgeR"); axA.bar(x+bw/2,dA,bw,color=CD,label="DESeq2")
for xi,e,d in zip(x,eA,dA):
    axA.text(xi-bw/2,e,f"{e/1000:.0f}k",ha="center",va="bottom",fontsize=8,color=CE,fontweight="bold")
    axA.text(xi+bw/2,d,f"{d/1000:.0f}k",ha="center",va="bottom",fontsize=8,color=CD,fontweight="bold")
axA.set_xticks(x); axA.set_xticklabels([TPN[t] for t in TPS]); axA.set_ylabel("strain-specific calls (sum over 16 strains)",fontsize=9.5)
axA.set_ylim(top=max(eA)*1.18); axA.legend(frameon=False,fontsize=9); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  Real-data calls — edgeR calls ~20–25% more\n(those extra calls are mostly false positives → panel B)",fontsize=10.5,fontweight="bold",loc="left")
# B: NULL false-positives (mean±SD), LOG
eB=[perm[tp].edgeR_FP.mean() for tp in TPS]; eBs=[perm[tp].edgeR_FP.std() for tp in TPS]
dB=[perm[tp].DESeq2_FP.mean() for tp in TPS]; dBs=[perm[tp].DESeq2_FP.std() for tp in TPS]
axB.bar(x-bw/2,eB,bw,yerr=eBs,color=CE,capsize=3,error_kw=dict(elinewidth=0.8),label="edgeR")
axB.bar(x+bw/2,dB,bw,yerr=dBs,color=CD,capsize=3,error_kw=dict(elinewidth=0.8),label="DESeq2")
axB.set_yscale("log"); axB.set_xticks(x); axB.set_xticklabels([TPN[t] for t in TPS])
axB.set_ylabel("NULL false-positive calls (log; mean±SD, 5 perms)",fontsize=9.5)
for xi,e,d in zip(x,eB,dB): axB.text(xi,max(e,d)*1.5,f"{e/d:.0f}×",ha="center",fontsize=9,fontweight="bold",color="#222")
axB.legend(frameon=False,fontsize=9); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Permutation NULL false-positives — edgeR ~10–21× more\n(strain labels shuffled → all calls are false; fewer = better)",fontsize=10.5,fontweight="bold",loc="left")
# C: null p<0.05 calibration vs target
eC=[perm[tp].edgeR_p05.mean() for tp in TPS]; eCs=[perm[tp].edgeR_p05.std() for tp in TPS]
dC=[perm[tp].DESeq2_p05.mean() for tp in TPS]; dCs=[perm[tp].DESeq2_p05.std() for tp in TPS]
axC.bar(x-bw/2,eC,bw,yerr=eCs,color=CE,capsize=3,error_kw=dict(elinewidth=0.8),label="edgeR")
axC.bar(x+bw/2,dC,bw,yerr=dCs,color=CD,capsize=3,error_kw=dict(elinewidth=0.8),label="DESeq2")
axC.axhline(0.05,ls="--",color="#444",lw=1.2,label="target = 0.05 (calibrated)")
axC.set_xticks(x); axC.set_xticklabels([TPN[t] for t in TPS]); axC.set_ylabel("fraction of NULL tests with p<0.05",fontsize=9.5)
axC.legend(frameon=False,fontsize=8.5); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  Calibration — edgeR ~2× anti-conservative; DESeq2 ≈ 0.05\n(well-calibrated under the null)",fontsize=10.5,fontweight="bold",loc="left")
# D: null p-value distribution (16.5dpc representative): should be uniform (flat at 0.05)
h=hist["16.5dpc"]; xb=h.bin_lo+0.025; w=0.022
axD.bar(xb-w/2,h.edgeR_frac,w,color=CE,label="edgeR"); axD.bar(xb+w/2,h.DESeq2_frac,w,color=CD,label="DESeq2")
axD.axhline(0.05,ls="--",color="#444",lw=1.2,label="uniform (well-calibrated)")
axD.set_xlabel("null p-value",fontsize=9.5); axD.set_ylabel("fraction of tests (E16.5)",fontsize=9.5)
axD.legend(frameon=False,fontsize=8.5); axD.spines[["top","right"]].set_visible(False)
axD.set_title("D  Null p-value distribution (E16.5) — edgeR spikes near 0\n(excess false signal); DESeq2 ≈ flat/uniform",fontsize=10.5,fontweight="bold",loc="left")
fig.suptitle("edgeR vs DESeq2 for strain-specific piRNA differential abundance — DESeq2 is correctly calibrated (adopted)",fontsize=13,fontweight="bold",y=0.995)
fig.text(0.5,0.005,"One-vs-rest per strain, 16 strains × 3 reps · 300k-feature representative sample · 5 label permutations · benchmark design BioMNI 3/3-verified "
  "(Soneson&Delorenzi 2013, Schurch 2016, Love 2014, Robinson 2010). edgeR QL F-test is anti-conservative for this many-group one-vs-rest design.",ha="center",fontsize=7,color="#666")
fig.tight_layout(rect=[0,0.02,1,0.975])
out=f"{T}/figures/Fig_benchmark_da"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
# source data
sd=pd.DataFrame({"tp":[TPN[t] for t in TPS],"real_edgeR":eA,"real_DESeq2":dA,
  "null_FP_edgeR":eB,"null_FP_DESeq2":dB,"null_p05_edgeR":eC,"null_p05_DESeq2":dC,
  "null_KS_edgeR":[perm[t].edgeR_KS.mean() for t in TPS],"null_KS_DESeq2":[perm[t].DESeq2_KS.mean() for t in TPS]})
sd.to_csv(f"{T}/data/source_data/SourceData_Fig_benchmark_da.csv",index=False)
print("wrote",out); print(sd.to_string(index=False))
