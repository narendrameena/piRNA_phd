#!/usr/bin/env python3
"""THEME 18 — TE family of origin of the THREE within-tp unique mechanisms (DESeq2 stage-peak), feeding in the
NEW stage-shifted (heterochronic) class. Method = theme-08 Fig_TE_private_families16, extended to 3 classes on
the DESeq2 stage-peak within-tp set: per strain, intersect each class's candidate production loci (cand_self16
BAM, own genome, strip PanSN '#' prefix -> chrN) with the per-strain RepeatMasker BED (col4 = name|class/family);
assign each locus its largest-overlap TE family. CAVEAT: cand_self16 index = main chr + MT -> TE fraction is a
LOWER BOUND (same as theme 08). Compares strain-private (insertion) / conserved-but-silent (regulatory) /
stage-shifted (heterochronic): TE-derived fraction, family spectrum, classical-vs-wild.
Panels: A TE-derived % per mechanism (per-strain dots); B top TE-family spectrum per mechanism (heatmap, % of
class TE-loci); C stage-shifted TE-derived: classical vs wild.
Data: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz + cand_self16/*.bam + resources/repeatMasker/*.bed"""
import sys, subprocess, tempfile, os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np, pysam
from collections import Counter, defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mc
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
D=f"{U}/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"; PG=f"{U}/pangenome_te"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]
SP="unique: strain-private locus"; CBS="unique: conserved-but-silent"; SH="unique: stage-shifted (heterochronic)"
LAB={SP:"strain-private",CBS:"conserved-but-silent",SH:"stage-shifted"}; ORDER=["strain-private","conserved-but-silent","stage-shifted"]
KCOL={"strain-private":"#7a3b9a","conserved-but-silent":"#0072B2","stage-shifted":"#009E73"}
cache=f"{T}/data/SourceData_Fig_te_origin_stageshift.csv"
if os.path.exists(cache):
    Tb=pd.read_csv(cache); print("loaded cache")
else:
    d=pd.read_csv(f"{D}/deseq_stagepeak_classified.csv.gz")
    rows=[]
    for X in CANON:
        rm=f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
        if not os.path.exists(rm): print("no RM",X); continue
        g=d[(d.strain==X)&(d.klass.isin([SP,CBS,SH]))]; cls={r.cand_id:r.klass for r in g.itertuples()}
        if not cls: continue
        bam=pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam","rb")
        tb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG); seen=set()
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.query_name not in cls or a.query_name in seen: continue
            seen.add(a.query_name); tb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
        tb.close(); bam.close()
        out=subprocess.run(f"sort -k1,1 -k2,2n {tb.name} | {BT} intersect -a - -b {rm} -wa -wb",shell=True,capture_output=True,text=True).stdout
        os.unlink(tb.name)
        best={}
        for ln in out.splitlines():
            f=ln.split("\t")
            if len(f)<8: continue
            ov=min(int(f[2]),int(f[6]))-max(int(f[1]),int(f[5])); fam=f[7].split("|")[-1] if "|" in f[7] else f[7]
            if f[3] not in best or ov>best[f[3]][0]: best[f[3]]=(ov,fam)
        nseen=Counter(); nte=Counter(); cnt=defaultdict(Counter)
        for cid in seen: nseen[cls[cid]]+=1
        for cid,(ov,fam) in best.items(): cnt[cls[cid]][fam]+=1; nte[cls[cid]]+=1
        for k in [SP,CBS,SH]:
            for fam,n in cnt[k].items(): rows.append((LAB[k],fam,X,n))
            rows.append((LAB[k],"__nseen__",X,nseen[k])); rows.append((LAB[k],"__nte__",X,nte[k]))
        print(f"{X}: "+" ".join(f"{LAB[k]}={nseen[k]}(TE {nte[k]})" for k in [SP,CBS,SH]),flush=True)
    Tb=pd.DataFrame(rows,columns=["klass","family","strain","count"]); Tb.to_csv(cache,index=False)
# ---- metrics ----
def tefrac_perstrain(klab):
    a=Tb[(Tb.klass==klab)&(Tb.family=="__nte__")].set_index("strain")["count"]; b=Tb[(Tb.klass==klab)&(Tb.family=="__nseen__")].set_index("strain")["count"]
    return {X:(100*a.get(X,0)/b.get(X,np.nan) if b.get(X,0) else np.nan) for X in CANON}
def tefrac_agg(klab):
    a=Tb[(Tb.klass==klab)&(Tb.family=="__nte__")]["count"].sum(); b=Tb[(Tb.klass==klab)&(Tb.family=="__nseen__")]["count"].sum()
    return 100*a/b if b else np.nan, int(a), int(b)
fr={k:tefrac_perstrain(k) for k in ORDER}; agg={k:tefrac_agg(k) for k in ORDER}
# top family spectrum (% of each class's TE loci)
fams=Tb[~Tb.family.isin(["__nseen__","__nte__"])]
top=fams.groupby("family")["count"].sum().sort_values(ascending=False).head(13).index.tolist()
M=pd.DataFrame(index=top,columns=ORDER,dtype=float)
for k in ORDER:
    tot=fams[fams.klass==k]["count"].sum()
    for fm in top: M.loc[fm,k]=100*fams[(fams.klass==k)&(fams.family==fm)]["count"].sum()/tot if tot else 0
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(15,6),dpi=300); gs=fig.add_gridspec(1,3,width_ratios=[1.1,1.5,0.9],wspace=0.42)
# A: TE-derived % per mechanism
axA=fig.add_subplot(gs[0]); xb=np.arange(len(ORDER))
for i,k in enumerate(ORDER):
    axA.bar(i,agg[k][0],0.62,color=KCOL[k],zorder=2)
    vals=[fr[k][s] for s in CANON if not np.isnan(fr[k][s])]
    axA.scatter([i]*len(vals),vals,s=10,color="#222",alpha=0.5,zorder=3)
    axA.text(i,agg[k][0],f"{agg[k][0]:.0f}%\n({agg[k][1]} of {agg[k][2]})",ha="center",va="bottom",fontsize=7.5,fontweight="bold")
axA.set_xticks(xb); axA.set_xticklabels([k.replace("-","-\n",1) for k in ORDER],fontsize=8); axA.set_ylabel("% TE-derived (lower bound)",fontsize=9.5)
axA.set_ylim(0,max(agg[k][0] for k in ORDER)*1.3); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  TE-derived fraction per\nwithin-tp unique mechanism",fontsize=10.5,fontweight="bold",loc="left")
# B: family spectrum heatmap
axB=fig.add_subplot(gs[1]); CMAP=mc.LinearSegmentedColormap.from_list("v",["#f7f7f7","#9ecae8","#3a8fd4","#08306b"])
im=axB.imshow(M.values.astype(float),aspect="auto",cmap=CMAP,vmin=0,vmax=np.nanmax(M.values))
axB.set_xticks(range(len(ORDER))); axB.set_xticklabels(ORDER,rotation=20,ha="right",fontsize=8.5)
axB.set_yticks(range(len(top))); axB.set_yticklabels(top,fontsize=7.5)
for i in range(len(top)):
    for j in range(len(ORDER)):
        v=M.values[i,j]
        if v>0: axB.text(j,i,f"{v:.0f}",ha="center",va="center",fontsize=6.5,color="white" if v>np.nanmax(M.values)*0.55 else "#222")
axB.set_title("B  TE-family spectrum (% of each mechanism's TE-derived loci)",fontsize=10,fontweight="bold",loc="left")
fig.colorbar(im,ax=axB,fraction=0.04,pad=0.02).ax.tick_params(labelsize=6)
# C: stage-shifted classical vs wild TE-derived
axC=fig.add_subplot(gs[2])
cl=[fr["stage-shifted"][s] for s in CANON if s not in WILD and not np.isnan(fr["stage-shifted"][s])]
wl=[fr["stage-shifted"][s] for s in CANON if s in WILD and not np.isnan(fr["stage-shifted"][s])]
bp=axC.boxplot([cl,wl],labels=["classical","wild"],patch_artist=True,widths=0.55)
for patch,c in zip(bp["boxes"],["#0072B2","#D55E00"]): patch.set_facecolor(c); patch.set_alpha(0.6)
for i,(grp,c) in enumerate(zip([cl,wl],["#0072B2","#D55E00"]),1): axC.scatter(np.random.normal(i,0.05,len(grp)),grp,s=14,color=c,zorder=3,edgecolor="white",linewidth=0.3)
axC.set_ylabel("% TE-derived (stage-shifted)",fontsize=9.5); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  Stage-shifted TE-origin:\nclassical vs wild",fontsize=10,fontweight="bold",loc="left")
fig.suptitle("TE family of origin of the within-tp unique mechanisms — feeding in the new stage-shifted (heterochronic) class",fontsize=12.5,fontweight="bold",y=1.02)
fig.tight_layout(rect=[0,0,1,0.96])
out=f"{T}/figures/Fig_te_origin_stageshift"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("\n=== TE-derived fraction per mechanism (aggregated, lower bound) ===")
for k in ORDER: print(f"  {k}: {agg[k][0]:.1f}%  ({agg[k][1]}/{agg[k][2]} mappable loci)")
print("top families per mechanism (%):"); print(M.round(1).to_string())
print("wrote",out)
