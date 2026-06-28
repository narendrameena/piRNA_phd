#!/usr/bin/env python3
"""THEME 22 figure 4 — the 'present-in-most-but-not-GRCm39' subset. (A) 117 strain-entries collapse to 22 distinct loci,
one present in all 16 strains; (B) distinct TE profile — SINE/B4-enriched and much YOUNGER than the bulk non-reference
set; (C) WHY absent from GRCm39 — deletion vs assembly-gap vs technical vs divergent; (D) the headline: a few are
high-expression piRNA clusters the reference entirely lacks (incl. C57BL/6J-lineage losses)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
T="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav"; D=f"{T}/data"
loci=pd.read_csv(f"{D}/shared_subset_loci.csv"); ent=pd.read_csv(f"{D}/shared_subset_entries.csv")
st=pd.read_csv(f"{D}/shared_subset_grcm39_status.csv")
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(12.8,9.4),dpi=300)
# A: distinct-strains-per-locus
vc=loci.n_strains.value_counts().sort_index()
axA.bar(vc.index,vc.values,color="#4393C3",edgecolor="white")
axA.set_xlabel("# strains sharing the locus (of 16)",fontsize=9); axA.set_ylabel("distinct loci",fontsize=9.2)
axA.spines[["top","right"]].set_visible(False)
axA.text(0.5,0.9,f"117 strain-entries → {len(loci)} distinct loci\n(flagged non-reference: did not halLiftover to GRCm39)",transform=axA.transAxes,ha="center",va="top",fontsize=7.6,color="#444",style="italic")
axA.set_title("A  Shared across the panel, flagged non-reference\n(117 entries → 22 loci)",fontsize=9.4,fontweight="bold",loc="left")
# B: TE family + age
fam=ent.te_family.value_counts().head(7)[::-1]
TC={"LTR":"#6a3d9a","LINE":"#E69F00","SINE":"#33a02c","DNA":"#b15928","Low":"#999","Sim":"#bbb"}
axB.barh(range(len(fam)),fam.values,color=[TC.get(str(f).split("/")[0][:3],"#888") for f in fam.index],edgecolor="white")
axB.set_yticks(range(len(fam))); axB.set_yticklabels(fam.index,fontsize=7)
axB.set_xlabel("strain-entries",fontsize=9); axB.spines[["top","right"]].set_visible(False)
axB.text(0.97,0.40,f"SINE/B4-enriched (vs ERVK/L1 in the bulk set)\nand YOUNGER: median {ent.te_div.median():.1f}% div\n(bulk non-ref 14.3%, reference 16.4%)",transform=axB.transAxes,ha="right",va="top",fontsize=7.2,color="#1B7837",fontweight="bold")
axB.set_title("B  A distinct, younger TE profile",fontsize=9.4,fontweight="bold",loc="left")
# C: GRCm39 status
order=["C57BL6J_DELETION","ASSEMBLY_GAP","PRESENT_in_GRCm39","DIVERGENT_absent"]
lab={"C57BL6J_DELETION":"C57BL/6J\ndeletion","ASSEMBLY_GAP":"assembly\ngap","PRESENT_in_GRCm39":"present\n(technical)","DIVERGENT_absent":"divergent /\nabsent"}
colmap={"C57BL6J_DELETION":"#C0392B","ASSEMBLY_GAP":"#E69F00","PRESENT_in_GRCm39":"#999","DIVERGENT_absent":"#4393C3"}
sc=st.classification.value_counts()
cats=[c for c in order if c in sc.index]+[c for c in sc.index if c not in order]
axC.bar(range(len(cats)),[sc[c] for c in cats],color=[colmap.get(c,"#777") for c in cats],edgecolor="white")
axC.set_xticks(range(len(cats))); axC.set_xticklabels([lab.get(c,c) for c in cats],fontsize=7.4)
axC.set_ylabel("distinct loci",fontsize=9.2); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  Why absent from GRCm39?",fontsize=9.4,fontweight="bold",loc="left")
# D: honest synthesis — mostly liftover artifact, a real minority are losses
axD.axis("off")
ndel=int((st.classification=="C57BL6J_DELETION").sum()); npres=int((st.classification=="PRESENT_in_GRCm39").sum())
dels=st[st.classification=="C57BL6J_DELETION"].sort_values("uniqFPM",ascending=False).head(4)
dl="\n".join(f"   chr{r.chrom}:{r.start/1e6:.0f}Mb  {str(r.te_family)[:14]}  {r.uniqFPM:.0f} FPM" for _,r in dels.iterrows())
axD.text(0.5,0.5,f"MOSTLY A LIFTOVER ARTIFACT — a real minority are losses\n\n"
  f"{npres}/22 are PRESENT in GRCm39 (cluster qcov≈1.0): the\nTE-rich body fails to halLiftover → falsely 'non-reference'\n(incl. the WSB chr4 L1, 247 FPM — it IS in the reference).\n\n"
  f"{ndel}/22 are genuine C57BL/6J-lineage DELETIONS — clusters\nthe panel kept but the reference lineage lost:\n{dl}\n\n"
  f"0 assembly gaps · 1 divergent.\n→ this subset is mostly a liftover limitation for repetitive\nclusters; a pangenome (odgi inject) sidesteps both.",
  transform=axD.transAxes,ha="center",va="center",fontsize=7.7,color="#222",bbox=dict(boxstyle="round,pad=0.6",fc="#f5f5f0",ec="#bbb"))
fig.suptitle("'Present-in-most-but-not-GRCm39' clusters: mostly a liftover artifact for TE-rich clusters, plus 6 genuine C57BL/6J-lineage losses",fontsize=10.0,fontweight="bold",y=0.99)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{T}/figures/Fig_shared_subset"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out)
