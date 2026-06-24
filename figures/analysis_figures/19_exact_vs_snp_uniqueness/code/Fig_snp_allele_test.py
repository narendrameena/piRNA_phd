#!/usr/bin/env python3
"""THEME 19 — TEST: are the piRNAs that the EXACT-sequence definition counts as 'unique' but the SNP-aware
definition rejects (the 4,394 SNP-alleles) genuine strain innovation, or standing genetic variation?
For each of the 4,394, count how many OTHER strains express a 1-3 mm allele (same stage) and the mismatch
distance. Result: they are shared (as a 1-3 mm allele) with ~11 of the other 15 strains on average, 88 % differ
by a single SNP → standing variation in widely-shared/conserved piRNAs, NOT strain-specific innovation. Supports
the SNP-aware definition (BioMNI 3/3). Panels: A #other strains sharing a 1-3mm allele; B mismatch distance;
C exact-unique composition (clean innovation vs SNP-allele standing variation).
Data: data/exact_stagepeak_classified.csv.gz + unique16/snp_variant_refinement_withintp.csv."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
CADD="#e7298a"; CSA="#1b9e77"
d=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz")
snpv=set(d[d.was_snp_variant].cand_id); n_clean=int(d.klass.isin(["unique: conserved-but-silent","unique: strain-private locus","unique: stage-shifted (heterochronic)"]).sum())
snp=pd.read_csv(f"{U}/unique16/snp_variant_refinement_withintp.csv"); snp=snp[snp.cand_id.isin(snpv)]
nshare=snp.groupby("cand_id").variant_strain.nunique()        # #other strains with a 1-3mm allele
mm=snp.sort_values("mm").drop_duplicates("cand_id").set_index("cand_id").mm
fig,(axA,axB,axC)=plt.subplots(1,3,figsize=(15,5),dpi=300,gridspec_kw=dict(width_ratios=[1.3,1,1],wspace=0.32))
# A: # other strains sharing a 1-3mm allele
vc=nshare.value_counts().sort_index();
axA.bar(vc.index,vc.values,color=CADD,width=0.85)
axA.axvline(nshare.mean(),ls="--",color="#333",lw=1.2); axA.text(nshare.mean()+0.3,vc.max()*0.9,f"mean {nshare.mean():.1f}\nof 15 other strains",fontsize=8.5,color="#333",fontweight="bold")
axA.set_xlabel("# OTHER strains expressing a 1–3 mm allele (same stage)",fontsize=9.5); axA.set_ylabel("SNP-allele candidates",fontsize=9.5)
axA.set_xticks(range(1,16)); axA.tick_params(labelsize=7.5); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  The exact-'unique' SNP-alleles are shared (±1–3 SNP)\nwith MANY other strains — not strain-specific",fontsize=10,fontweight="bold",loc="left")
# B: mismatch distance
mc=mm.value_counts().reindex([1,2,3],fill_value=0)
axB.bar([1,2,3],mc.values,color=["#fb6a4a","#cb181d","#67000d"],width=0.6)
for m,v in zip([1,2,3],mc.values): axB.text(m,v,f"{v:,}\n({100*v/mc.sum():.0f}%)",ha="center",va="bottom",fontsize=8.5,fontweight="bold")
axB.set_xticks([1,2,3]); axB.set_xlabel("mismatches to nearest expressed allele",fontsize=9.5); axB.set_ylabel("SNP-allele candidates",fontsize=9.5)
axB.set_ylim(top=mc.max()*1.18); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  88 % differ by a SINGLE SNP\n(standing variation, not new sequence)",fontsize=10,fontweight="bold",loc="left")
# C: exact-unique composition
axC.bar([0],[n_clean],0.6,color=CSA,label=f"clean innovation\n(no 1–3 SNP allele elsewhere): {n_clean:,}")
axC.bar([0],[len(snpv)],0.6,bottom=[n_clean],color=CADD,label=f"SNP-allele / standing variation\n(shared ±1–3 SNP w/ ~{nshare.mean():.0f} strains): {len(snpv):,}")
axC.text(0,n_clean+len(snpv),f"EXACT 'unique' = {n_clean+len(snpv):,}",ha="center",va="bottom",fontsize=9,fontweight="bold")
axC.text(0,n_clean/2,f"{100*n_clean/(n_clean+len(snpv)):.0f}%",ha="center",va="center",fontsize=10,color="white",fontweight="bold")
axC.text(0,n_clean+len(snpv)/2,f"{100*len(snpv)/(n_clean+len(snpv)):.0f}%",ha="center",va="center",fontsize=10,color="white",fontweight="bold")
axC.set_xticks([0]); axC.set_xticklabels(["exact-sequence\nunique"],fontsize=9); axC.set_ylabel("genuinely-unique piRNAs",fontsize=9.5)
axC.set_ylim(top=(n_clean+len(snpv))*1.15); axC.legend(fontsize=7.3,frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.12))
axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  29 % of exact-'unique' is standing\nvariation → SNP-aware preferred (BioMNI 3/3)",fontsize=10,fontweight="bold",loc="left")
fig.suptitle("TEST: exact-sequence 'unique' inflates by SNP-alleles of widely-shared piRNAs (standing variation), not by innovation",fontsize=12,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.95])
out=f"{T}/figures/Fig_snp_allele_test"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pd.DataFrame({"cand_id":nshare.index,"n_other_strains_with_allele":nshare.values}).to_csv(f"{T}/data/SourceData_Fig_snp_allele_test.csv.gz",index=False)
print(f"wrote {out} | SNP-alleles {len(snpv):,}, mean sharing {nshare.mean():.2f} strains, mm {dict(mc)}")
