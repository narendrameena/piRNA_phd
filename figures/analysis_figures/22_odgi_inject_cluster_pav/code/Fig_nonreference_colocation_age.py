#!/usr/bin/env python3
"""THEME 22 figure 3 — non-reference clusters are YOUNG, strain-private insertions vs OLD, conserved ones.
(A) co-location: 66% strain-private, 34% shared (sharing follows phylogeny); (B) TE age — non-reference cluster TEs are
younger (lower RepeatMasker %div) than reference/conserved cluster TEs; (C) holds WITHIN the major driver families
(L1, ERVL-MaLR, ERV1), so it is insertion age not family composition; (D) the synthesis."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
T="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav"; D=f"{T}/data"
co=pd.read_csv(f"{D}/colocation.csv"); dist=json.load(open(f"{D}/te_age_dists.json")); famdf=pd.read_csv(f"{D}/te_age_byfamily.csv")
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(12.8,9.4),dpi=300)
# A: co-location sharing degree
vc=co.n_other_seq.value_counts().sort_index(); xs=list(range(16))
ys=[int(vc.get(k,0)) for k in xs]; cols=["#C0392B"]+["#4393C3"]*15
axA.bar(xs,ys,color=cols,edgecolor="white")
axA.set_xlabel("# OTHER strains carrying the same insertion",fontsize=9); axA.set_ylabel("non-reference clusters",fontsize=9.2)
axA.spines[["top","right"]].set_visible(False); axA.set_xticks([0,5,10,15])
axA.legend(handles=[Patch(fc="#C0392B",label="strain-private (66%)"),Patch(fc="#4393C3",label="shared (34%)")],fontsize=7.6,frameon=False,loc="upper center")
axA.text(0.97,0.55,"916 private (66%) · 477 shared (34%)\nsharing follows phylogeny:\n62% within wild/classical group",transform=axA.transAxes,ha="right",fontsize=7.3,color="#444",style="italic")
axA.set_title("A  Co-location: most are strain-private insertions\n(shared ones track phylogeny)",fontsize=9.4,fontweight="bold",loc="left")
# B: TE age overall
nr=np.array(dist["nonref_div"]); rf=np.array(dist["ref_div"])
bp=axB.boxplot([rf,nr],positions=[0,1],widths=0.55,patch_artist=True,showfliers=False,vert=True)
for p,c in zip(bp["boxes"],["#cccccc","#1B7837"]): p.set_facecolor(c)
for m in bp["medians"]: m.set_color("#111")
axB.set_xticks([0,1]); axB.set_xticklabels([f"reference /\nconserved\n(n={len(rf):,})",f"non-reference\n(n={len(nr):,})"],fontsize=8.4)
axB.set_ylabel("TE divergence from consensus (%)\n← YOUNGER          OLDER →",fontsize=8.8); axB.spines[["top","right"]].set_visible(False)
axB.text(0.5,0.93,f"non-ref median {np.median(nr):.1f}% vs reference {np.median(rf):.1f}%\nnon-reference TEs are YOUNGER (p=7×10⁻⁶⁸)",transform=axB.transAxes,ha="center",va="top",fontsize=7.5,color="#1B7837",fontweight="bold")
axB.set_title("B  Non-reference clusters sit on YOUNGER TEs\nthan conserved clusters",fontsize=9.4,fontweight="bold",loc="left")
# C: within-family TE age
fams=["LINE/L1","LTR/ERVK","LTR/ERVL-MaLR","LTR/ERV1"]; sub=famdf.set_index("family").loc[fams]
x=np.arange(len(fams)); w=0.38
axC.bar(x-w/2,sub.ref_div,w,color="#cccccc",edgecolor="white",label="reference")
axC.bar(x+w/2,sub.nonref_div,w,color="#1B7837",edgecolor="white",label="non-reference")
for i,(f,r) in enumerate(sub.iterrows()):
    if r.p<0.05: axC.text(i,max(r.ref_div,r.nonref_div)+0.4,"*",ha="center",fontsize=13,color="#B2182B")
axC.set_xticks(x); axC.set_xticklabels([f.replace("/","/\n") for f in fams],fontsize=7.6)
axC.set_ylabel("median TE divergence (%)",fontsize=9); axC.spines[["top","right"]].set_visible(False); axC.legend(fontsize=7.6,frameon=False,loc="upper left")
axC.text(0.5,0.04,"younger non-ref TEs WITHIN families (L1, ERVL-MaLR, ERV1: *)\n→ insertion age, not family composition. ERVK already youngest (n.s.)",transform=axC.transAxes,ha="center",fontsize=6.7,color="#444",style="italic")
axC.set_title("C  Holds within the driver families",fontsize=9.4,fontweight="bold",loc="left")
# D: synthesis
axD.axis("off")
axD.text(0.5,0.5,"SYNTHESIS — two ages of piRNA cluster\n\n"
  "CONSERVED clusters  →  OLD TEs, shared across strains\n  the ancient pachytene piRNA factories\n\n"
  "NON-REFERENCE clusters  →  YOUNG TEs, strain-private\n  recent LTR/ERVK + LINE-1 insertions, post-divergence,\n  caught by the host piRNA response (arms race)\n  = the leading edge of piRNA-cluster evolution\n\n"
  "Both predictions confirmed:\n  • 66% strain-private (co-location)\n  • younger TEs overall (p=7×10⁻⁶⁸) + within families\n  • novelty tracks the genome-wide TE-insertion burden (ρ=0.61)",
  transform=axD.transAxes,ha="center",va="center",fontsize=8.2,color="#222",
  bbox=dict(boxstyle="round,pad=0.7",fc="#f5f5f0",ec="#bbb"))
fig.suptitle("Non-reference piRNA clusters are YOUNG, strain-private TE insertions; conserved clusters are OLD and shared",fontsize=10.6,fontweight="bold",y=0.99)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{T}/figures/Fig_nonreference_colocation_age"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out)
