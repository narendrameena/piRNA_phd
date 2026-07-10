#!/usr/bin/env python3
"""Figure: SNP-variant producer method test (129S1->C57, delivered n=330). See SNP_VARIANT_METHOD_TEST.md."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/figures"
# measured on the 129S1->C57 task (see SNP_VARIANT_METHOD_TEST.md); numbers frozen here so the figure is self-contained
g={"genome_recall_pct":84.5455,"genome_rec":279,"n_deliv":330,"miss":51,"no_align":43,"genome_ne":8}
bt={"clean":9931,"substr":41905,"total":51836}
STAR={"same":9549,"substr":41072,"gapped":34086,"rev":23527}; STAR_tot=sum(STAR.values())
BT={"same":bt["clean"],"substr":bt["substr"],"gapped":0,"rev":0}; BT_tot=bt["total"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9})
fig=plt.figure(figsize=(12.6,4.8),dpi=300)
gs=fig.add_gridspec(1,3,wspace=0.42,width_ratios=[1.15,0.85,1.05],bottom=0.22,top=0.82)

# ---- A: recall of the delivered SNP-variant set ----
axA=fig.add_subplot(gs[0,0])
labels=["genomic STAR\n(classify_step416)","bowtie -v3\n(direct pool)","STAR\n(direct pool)"]
recall=[g["genome_recall_pct"],100.0,100.0]; colors=["#D55E00","#009E73","#0072B2"]
x=np.arange(3); axA.bar(x,recall,color=colors,edgecolor="white",width=0.72)
for xi,r in zip(x,recall): axA.text(xi,r+1.6,f"{r:.1f}%",ha="center",fontweight="bold",fontsize=10)
axA.axhline(100,ls=":",c="#999",lw=0.8)
axA.set_xticks(x); axA.set_xticklabels(labels,fontsize=7.6)
for t,c in zip(axA.get_xticklabels(),colors): t.set_color(c); t.set_fontweight("bold")
axA.set_ylabel("recall of delivered SNP-variants\n(129S1 -> C57, n = 330)",fontsize=8.6)
axA.set_ylim(0,116)
axA.set_title("A   Reference determines correctness\n(orange = genome 84.5%;  green/blue = expressed-pool 100%)",fontsize=8.4,fontweight="bold",loc="left")
axA.spines[["top","right"]].set_visible(False)

# ---- B: why the genome proxy misses ----
axB=fig.add_subplot(gs[0,1])
sizes=[g["no_align"],g["genome_ne"]]
axB.pie(sizes,labels=[f"no alignment\nto genome\n(n={g['no_align']})",f"aligns, but\ngenome ≠ expressed\n(n={g['genome_ne']})"],
        colors=["#E69F00","#CC79A7"],autopct=lambda p:f"{p:.0f}%",startangle=90,
        textprops={"fontsize":7.6},wedgeprops={"edgecolor":"white","linewidth":1.2})
axB.set_title(f"B   Why the genome proxy misses\n({g['miss']} missed; genomic presence ≠ expression)",fontsize=9.0,fontweight="bold",loc="left")

# ---- C: aligner output composition ----
axC=fig.add_subplot(gs[0,2])
cats=["same-length\n(usable)","substring\n(length-isoform)","gapped-indel\n(wrong for SNP)","reverse\n(wrong strand)"]
cc=["#009E73","#bbbbbb","#D55E00","#7a3b9a"]
tools=["bowtie -v3","STAR"]; x2=np.arange(2)
btv=np.array([100*BT[k]/BT_tot for k in ("same","substr","gapped","rev")])
stv=np.array([100*STAR[k]/STAR_tot for k in ("same","substr","gapped","rev")])
for i,(label,col) in enumerate(zip(cats,cc)):
    axC.bar(0,btv[i],bottom=btv[:i].sum(),color=col,width=0.62,label=label,edgecolor="white",linewidth=0.4)
    axC.bar(1,stv[i],bottom=stv[:i].sum(),color=col,width=0.62,edgecolor="white",linewidth=0.4)
axC.text(0,btv[0]/2,f"{btv[0]:.0f}%",ha="center",color="white",fontweight="bold",fontsize=8.5)
axC.text(1,stv[0]/2,f"{stv[0]:.0f}%",ha="center",color="white",fontweight="bold",fontsize=8.5)
axC.set_xticks(x2); axC.set_xticklabels(tools,fontweight="bold"); axC.set_ylabel("% of output records",fontsize=8.6); axC.set_ylim(0,100)
axC.set_title("C   bowtie: exhaustive (635 found), subst-only;\nSTAR: non-exhaustive (631, missed 4) + gapped/rev",fontsize=8.2,fontweight="bold",loc="left",pad=8)
axC.legend(fontsize=5.9,loc="upper center",bbox_to_anchor=(0.5,-0.14),ncol=2,frameon=False,handlelength=1.2,columnspacing=1.0)
axC.spines[["top","right"]].set_visible(False)

fig.suptitle("SNP-variant producer method test — the reference (expressed pool vs genome), not the tool, determines correctness",
             fontsize=10.2,fontweight="bold",y=0.96)
fig.text(0.5,0.015,"Task: which 129S1 candidate piRNAs are ≤3-substitution (Hamming) variants of a C57BL/6NJ-EXPRESSED piRNA (delivered n=330; C57 pool 34.9M seqs). "
  "Genomic-proxy (classify_step416) = 84.5%; direct pool-search = 100% by BOTH bowtie and STAR — so the delivered numbers are correct and the genome reference was the defect. "
  "Global genomic-proxy reproduction of the full delivered set: 50% (forward-only) -> 86% (minus-strand-fixed). bowtie -v3: 34.9M reads in 3 min.",
  ha="center",fontsize=6.3,color="#555",wrap=True)
for ext in ("pdf","png","svg"): fig.savefig(f"{OUT}/snp_method_comparison.{ext}",bbox_inches="tight")
print(f"wrote {OUT}/snp_method_comparison.{{pdf,png,svg}}")
print(f"A recall: genome={g['genome_recall_pct']:.1f}% bowtie=100% STAR=100%")
print(f"B miss: no_align={g['no_align']} genome_ne={g['genome_ne']}")
print(f"C bowtie usable={btv[0]:.0f}% STAR usable={stv[0]:.0f}% (STAR gapped={stv[2]:.0f}% reverse={stv[3]:.0f}%)")
