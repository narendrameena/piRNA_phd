#!/usr/bin/env python3
"""THEME 18 — effect of SNP-variants on the within-tp unique data (user request 2026-06-23).
A SNP-variant = a stage-peak candidate that aligns 1-3 mm to another strain's genome where that allele IS
expressed AT THE SAME stage -> it only LOOKS strain-specific; it is a strain ALLELE of a same-stage piRNA, NOT
genuinely unique. SNP-refinement removes them, PURIFYING the unique set. Panels:
  A  purifying effect: genuinely-unique WITHOUT SNP-refinement (SNP-variants counted as conserved-but-silent)
     vs WITH refinement, per stage -> how much the unique set shrinks.
  B  mismatch-count distribution (1 / 2 / 3 mm; per candidate = closest allele).
  C  SNP position along the piRNA (5'->3') -> where strain alleles differ.
  D  transition vs transversion of the strain SNPs.
Data: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz + unique16/snp_variant_refinement_withintp.csv"""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#B2182B"}
SNPK="SNP-variant (1-3mm, same stage)"; CBS="unique: conserved-but-silent"; SP="unique: strain-private locus"; SH="unique: stage-shifted (heterochronic)"
GU=[CBS,SP,SH]
d=pd.read_csv(f"{U}/deseq16_lenfilt/deseq_stagepeak_classified.csv.gz")
snpv_ids=set(d[d.klass==SNPK].cand_id)
snp=pd.read_csv(f"{U}/unique16/snp_variant_refinement_withintp.csv")
snp=snp[snp.cand_id.isin(snpv_ids)].copy()                       # restrict to DESeq2 stage-peak SNP-variants
# closest allele per candidate (min mm)
best=snp.sort_values("mm").drop_duplicates("cand_id")
def mismatches(a,b):
    return [(i+1,a[i],b[i]) for i in range(min(len(a),len(b))) if a[i]!=b[i]]
TI={frozenset("AG"),frozenset("CT")}
pos=[]; titv=[]; n_incons=0
for _,r in best.iterrows():
    ms=mismatches(str(r.home_seq),str(r.Y_allele))
    if len(ms)!=r.mm:            # data-integrity guard: skip the ~80 (1.8%) alleles whose stored home_seq/Y_allele
        n_incons+=1; continue    # disagree with the recorded mm (else the Ti/Tv & position panels gain spurious SNPs)
    for i,hb,yb in ms:
        pos.append(i); titv.append("transition" if frozenset((hb,yb)) in TI else "transversion")
print(f"[panels C/D] used {len(best)-n_incons} of {len(best)} candidates for SNP-level analysis; dropped {n_incons} (home_seq/Y_allele inconsistent with mm)")
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(13.5,9.5),dpi=300); xb=np.arange(len(TPS))
# A: purifying effect
with_ref=[d[(d.timepoint==t)&(d.klass.isin(GU))].shape[0] for t in TPS]
snpn=[d[(d.timepoint==t)&(d.klass==SNPK)].shape[0] for t in TPS]
without=[w+s for w,s in zip(with_ref,snpn)]
axA.bar(xb-0.2,without,0.4,color="#bbbbbb",label="WITHOUT SNP-refinement (naive unique)")
axA.bar(xb+0.2,with_ref,0.4,color="#7a3b9a",label="WITH SNP-refinement (true unique)")
for xi,wo,wr,s in zip(xb,without,with_ref,snpn):
    axA.text(xi-0.2,wo,f"{wo:,}",ha="center",va="bottom",fontsize=8,color="#666")
    axA.text(xi+0.2,wr,f"{wr:,}",ha="center",va="bottom",fontsize=8,color="#7a3b9a",fontweight="bold")
    axA.text(xi,max(wo,wr)*1.06,f"−{s:,}\n(−{100*s/wo:.0f}%)",ha="center",fontsize=7.5,color="#b00")
axA.set_xticks(xb); axA.set_xticklabels([TPN[t] for t in TPS]); axA.set_ylabel("genuinely-unique piRNAs",fontsize=9.5)
axA.set_ylim(top=max(without)*1.2); axA.legend(fontsize=7.8,frameon=False); axA.spines[["top","right"]].set_visible(False)
axA.set_title(f"A  SNP-refinement removes {sum(snpn):,} strain alleles ({100*sum(snpn)/sum(without):.0f}% of naive unique)\n→ purifies the unique set",fontsize=10.5,fontweight="bold",loc="left")
# B: mm distribution (per candidate)
mmc=best.mm.value_counts().reindex([1,2,3],fill_value=0)
axB.bar([1,2,3],mmc.values,color=["#E69F00","#D55E00","#A33"],width=0.6)
for m,v in zip([1,2,3],mmc.values): axB.text(m,v,f"{v:,}\n({100*v/mmc.sum():.0f}%)",ha="center",va="bottom",fontsize=8,fontweight="bold")
axB.set_xticks([1,2,3]); axB.set_xlabel("mismatches to nearest same-stage allele",fontsize=9.5); axB.set_ylabel("SNP-variant candidates",fontsize=9.5)
axB.set_ylim(top=mmc.max()*1.18); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Most strain alleles differ by a single SNP",fontsize=10.5,fontweight="bold",loc="left")
# C: SNP position along piRNA
ph=pd.Series(pos).value_counts().sort_index()
axC.bar(ph.index,ph.values,color="#117733",width=0.85)
axC.axvspan(1.5,7.5,color="#E69F00",alpha=0.12,zorder=0); axC.text(4.5,axC.get_ylim()[1] if False else max(ph.values)*0.95,"seed (2–7)",ha="center",fontsize=7.5,color="#7a5500")
axC.set_xlabel("SNP position along piRNA (5'→3')",fontsize=9.5); axC.set_ylabel("SNPs",fontsize=9.5); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  Where strain alleles differ (SNP position)",fontsize=10.5,fontweight="bold",loc="left")
# D: ti/tv
tc=pd.Series(titv).value_counts()
ti=tc.get("transition",0); tv=tc.get("transversion",0)
axD.bar(["transition\n(A↔G, C↔T)","transversion"],[ti,tv],color=["#0072B2","#D55E00"],width=0.6)
for i,v in enumerate([ti,tv]): axD.text(i,v,f"{v:,}\n({100*v/(ti+tv):.0f}%)",ha="center",va="bottom",fontsize=8,fontweight="bold")
axD.set_ylabel("SNPs",fontsize=9.5); axD.set_ylim(top=max(ti,tv)*1.18); axD.spines[["top","right"]].set_visible(False)
axD.set_title(f"D  Transition/transversion = {ti/max(tv,1):.2f}  (genome-wide ~2)",fontsize=10.5,fontweight="bold",loc="left")
fig.suptitle("Effect of SNP-variants on the within-tp unique piRNA set — strain alleles of same-stage piRNAs, removed by SNP-refinement",fontsize=12.5,fontweight="bold",y=0.995)
fig.tight_layout(rect=[0,0,1,0.97])
out=f"{T}/figures/Fig_snp_variant_effect"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pd.DataFrame({"tp":[TPN[t] for t in TPS],"unique_with_refinement":with_ref,"snp_variants_removed":snpn,"naive_unique_without":without}).to_csv(f"{T}/data/source_data/SourceData_Fig_snp_variant_effect.csv",index=False)
best[["cand_id","home","variant_strain","mm"]].to_csv(f"{T}/data/source_data/SourceData_Fig_snp_variant_effect_perSNP.csv.gz",index=False)
print(f"wrote {out} | SNP-variants={len(snpv_ids):,} removed; mm: {dict(mmc)}; ti/tv={ti}/{tv}")
