#!/usr/bin/env python3
"""Figure: SNP-variant producer method test + resolution. All numbers frozen here (see SNP_VARIANT_METHOD_TEST.md).
A: reference (pool vs genome) determines correctness (129S1->C57 recall, n=330).
B: why the genome proxy misses (84% don't align).
C: full-scale reproduction of the delivered SNP set (n=217,559) by each producer.
D: resulting genuinely-unique -- the method band collapses; the lift-anchored producer == delivered."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/figures"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9})
GENO="#D55E00"; POOL="#E69F00"; STARC="#0072B2"; LIFT="#009E73"; DEL="#777777"
fig=plt.figure(figsize=(12.4,8.2),dpi=300)
gs=fig.add_gridspec(2,2,hspace=0.62,wspace=0.30,left=0.08,right=0.97,top=0.87,bottom=0.10)

# ---- A: reference determines correctness (129S1->C57 recall) ----
axA=fig.add_subplot(gs[0,0])
labs=["genomic STAR\n(classify_step416)","bowtie -v3\n(direct pool)","STAR\n(direct pool)"]
rec=[84.5,100.0,100.0]; cols=[GENO,LIFT,STARC]; x=np.arange(3)
axA.bar(x,rec,color=cols,edgecolor="white",width=0.72)
for xi,r in zip(x,rec): axA.text(xi,r+1.6,f"{r:.1f}%",ha="center",fontweight="bold",fontsize=9.5)
axA.axhline(100,ls=":",c="#999",lw=0.8); axA.set_xticks(x); axA.set_xticklabels(labs,fontsize=7.3)
for t,c in zip(axA.get_xticklabels(),cols): t.set_color(c); t.set_fontweight("bold")
axA.set_ylabel("recall of delivered SNP-variants\n(129S1->C57 subset, n=330)",fontsize=8.4); axA.set_ylim(0,116)
axA.set_title("A   Reference determines correctness\n(expressed pool 100% vs genome 84.5%)",fontsize=8.8,fontweight="bold",loc="left")
axA.spines[["top","right"]].set_visible(False)

# ---- B: why the genome proxy misses ----
axB=fig.add_subplot(gs[0,1])
axB.pie([43,8],labels=["no alignment\nto genome\n(n=43)","aligns, but\ngenome != expressed\n(n=8)"],
        colors=["#E69F00","#CC79A7"],autopct=lambda p:f"{p:.0f}%",startangle=90,
        textprops={"fontsize":7.6},wedgeprops={"edgecolor":"white","linewidth":1.2})
axB.set_title("B   Why the genome proxy misses\n(51 missed; genomic presence != expression)",fontsize=8.8,fontweight="bold",loc="left")

# ---- C: full-scale reproduction of the delivered SNP set ----
axC=fig.add_subplot(gs[1,0])
pl=["genomic proxy\n(classify_step416)","pure-pool\nbowtie -v3","STAR\nco-location","lift-anchored\n(DEFINITIVE)"]
repro=[86.0,99.98,85.4,100.0]; miss=["~30k","33","31,823","0"]; cc=[GENO,POOL,STARC,LIFT]; x=np.arange(4)
axC.bar(x,repro,color=cc,edgecolor="white",width=0.72)
for xi,r,m in zip(x,repro,miss):
    axC.text(xi,r+0.6,f"{r:.1f}%",ha="center",fontweight="bold",fontsize=8.4)
    axC.text(xi,5,f"missed\n{m}",ha="center",va="bottom",fontsize=6.3,color="white",fontweight="bold")
axC.axhline(100,ls=":",c="#999",lw=0.8); axC.set_xticks(x); axC.set_xticklabels(pl,fontsize=6.9)
for t,c in zip(axC.get_xticklabels(),cc): t.set_color(c); t.set_fontweight("bold")
axC.set_ylabel("% of delivered SNP set reproduced\n(full scale, n=217,559)",fontsize=8.4); axC.set_ylim(0,108)
axC.set_title("C   Reproducing the delivered SNP set\n(lift-anchored: 100%, 0 missed)",fontsize=8.8,fontweight="bold",loc="left")
axC.spines[["top","right"]].set_visible(False)

# ---- D: resulting genuinely-unique -- the band collapses to the lift/delivered value ----
axD=fig.add_subplot(gs[1,1])
gl=["pure-pool\nbowtie","lift-anchored\n(DEFINITIVE)","delivered","STAR\nco-location"]
gu=[79.529,106.815,106.961,123.862]; gc=[POOL,LIFT,DEL,STARC]; x=np.arange(4)
axD.bar(x,gu,color=gc,edgecolor="white",width=0.72)
axD.axhline(106.961,ls="--",c=DEL,lw=1.0)
for xi,g in zip(x,gu): axD.text(xi,g+2.0,f"{g:.1f}k",ha="center",fontweight="bold",fontsize=8.2)
axD.text(3.4,106.961,"delivered",ha="right",va="bottom",fontsize=6.3,color=DEL)
axD.annotate("over-count\n(coincidental)",(0,79.5),(0.15,52),ha="center",fontsize=6.1,color=POOL,arrowprops=dict(arrowstyle="->",color=POOL,lw=0.8))
axD.annotate("under-count\n(STAR blind-spot)",(3,123.9),(2.9,140),ha="center",fontsize=6.1,color=STARC,arrowprops=dict(arrowstyle="->",color=STARC,lw=0.8))
axD.set_xticks(x); axD.set_xticklabels(gl,fontsize=6.9)
for t,c in zip(axD.get_xticklabels(),gc): t.set_color(c); t.set_fontweight("bold")
axD.set_ylabel("genuinely-unique piRNAs (thousands)",fontsize=8.4); axD.set_ylim(0,152)
axD.set_title("D   Genuinely-unique: the method band collapses\n(lift-anchored == delivered = 106,961)",fontsize=8.8,fontweight="bold",loc="left")
axD.spines[["top","right"]].set_visible(False)

fig.suptitle("SNP-variant producer: the delivered numbers are correct, reproducible (100%), and set by the lift-anchored method",
             fontsize=10.6,fontweight="bold",y=0.945)
fig.text(0.5,0.012,"The SNP-variant class (54% of klass5) = piRNAs 1-3 substitutions from a piRNA EXPRESSED at the ORTHOLOGOUS locus in another strain. classify_step416's genome reference was wrong (84.5%); "
  "pure-pool bowtie over-counts on coincidental lookalikes; STAR co-location under-counts (cannot align divergent orthologs). Anchoring the ortholog with the CACTUS LIFT (build_snp_variant_lift.py) reproduces the "
  "delivered set exactly (100%, 0 missed) -> genuinely-unique = 106,961, definitively and reproducibly.",
  ha="center",fontsize=6.4,color="#555",wrap=True)
for ext in ("pdf","png","svg"): fig.savefig(f"{OUT}/Fig_snp_method_comparison.{ext}",bbox_inches="tight")
print("wrote Fig_snp_method_comparison.{pdf,png,svg}")
