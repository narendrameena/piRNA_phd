#!/usr/bin/env python3
"""16-strain TE families of origin for strain-private (genuinely-unique) piRNAs. For each strain, intersect
its strain-private-locus piRNA loci (cand_self16 BAM, own genome) with the per-strain RepeatMasker BED
(col4 = name|class/family), assign each locus its largest-overlap TE family, and count. (A) heatmap: top TE
families x 16 strains (log counts); (B) TE-derived fraction per strain (canonical order). CAVEAT: index =
main chr + MT only; private piRNAs on unplaced contigs are excluded -> TE fraction is a LOWER BOUND."""
import warnings; warnings.filterwarnings("ignore")
import sys,subprocess,tempfile,os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np, pysam
from collections import Counter, defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]
d=pd.read_csv(f"{U}/unique16/final_classified_clean.csv.gz")   # mm0-clean strain-private (klass5); cand_self16 already mm0-filters, so output is the clean set either way
fammat=defaultdict(dict); tefrac={}; nsp={}
for X in CANON:
    g=d[(d.strain==X)&(d.klass5=="unique: strain-private locus")].copy(); g["id"]=X+"|"+g.timepoint+"|"+g.sequence
    sp=set(g.id); nsp[X]=len(sp)
    rm=f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
    if not sp or not os.path.exists(rm): tefrac[X]=np.nan; continue
    bam=pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam","rb")
    tb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG); seen=set()
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.query_name not in sp or a.query_name in seen: continue
        seen.add(a.query_name); tb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
    tb.close(); bam.close()
    out=subprocess.run(f"sort -k1,1 -k2,2n {tb.name} | {BT} intersect -a - -b {rm} -wa -wb",shell=True,capture_output=True,text=True).stdout
    os.unlink(tb.name)
    best={}  # id -> (overlap,family)
    for ln in out.splitlines():
        f=ln.split("\t")
        if len(f)<8: continue
        ov=min(int(f[2]),int(f[6]))-max(int(f[1]),int(f[5])); fam=f[7].split("|")[-1] if "|" in f[7] else f[7]
        if f[3] not in best or ov>best[f[3]][0]: best[f[3]]=(ov,fam)
    fc=Counter(v[1] for v in best.values())
    for fam,n in fc.items(): fammat[fam][X]=n
    tefrac[X]=100*len(best)/len(seen) if seen else np.nan
    print(f"{X}: private={len(seen):,} TE-annotated={len(best):,} ({tefrac[X]:.0f}%) top={fc.most_common(3)}")
M=pd.DataFrame(fammat).T.reindex(columns=CANON).fillna(0)
M=M.loc[M.sum(1).sort_values(ascending=False).head(16).index]
M.to_csv(f"{PG}/SourceData_TE_private_families16.csv")
# ---- TE-family RNA-seq EXPRESSION at the SAME family granularity (Panel C) ----
# TE expression from SMALL RNA (results/TEtranscriptCount/all_featureCounts_TE.tab = sRNA reads on TE copies; NOT RNA-seq).
# Geneid = "<name> <class/family> <strain#1#chrom> <coords>" -> family = 2nd token; strain from col2; count col7.
famexpr=defaultdict(lambda: defaultdict(float))
F=f"{ROOT}/results/TEtranscriptCount/all_featureCounts_TE.tab"
if os.path.exists(F):
    for ln in open(F):
        if ln.startswith("#") or ln.startswith("Geneid"): continue
        c=ln.rstrip("\n").split("\t")
        if len(c)<7: continue
        parts=c[0].split()
        if len(parts)<2: continue
        X=c[1].split("#")[0]
        if X not in CANON: continue
        try: famexpr[parts[1]][X]+=float(c[6])
        except ValueError: pass
E=pd.DataFrame(famexpr).T.reindex(index=M.index,columns=CANON).fillna(0.0)   # same family rows/order as A
E.to_csv(f"{PG}/SourceData_TE_private_families16_expression.csv")
print("expression: strains with RNA-TE data =", sum(1 for X in CANON if (E[X]>0).any()), "/", len(CANON))
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(17,8.8),dpi=300)
gs=fig.add_gridspec(2,2,height_ratios=[3,1.05],hspace=0.45,wspace=0.16)
axH=fig.add_subplot(gs[0,0]); axE=fig.add_subplot(gs[0,1],sharey=axH); axB=fig.add_subplot(gs[1,:])
import matplotlib.colors as mc
_Mlog=np.log10(M.values+1); _cmap=plt.get_cmap("magma").copy(); _cmap.set_bad("white")   # empty (zero) cells -> white
im=axH.imshow(np.ma.masked_where(M.values==0,_Mlog),aspect="auto",cmap=_cmap,vmin=0,vmax=_Mlog.max())
axH.set_xticks(range(len(CANON))); axH.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.5)
axH.set_yticks(range(len(M))); axH.set_yticklabels(M.index,fontsize=7.5)
axH.set_xticks(np.arange(-.5,len(CANON),1),minor=True); axH.set_yticks(np.arange(-.5,len(M),1),minor=True)
axH.grid(which="minor",color="#e6e6e6",linewidth=0.5); axH.tick_params(which="minor",length=0)
for i in range(len(M)):
    for j in range(len(CANON)):
        v=int(M.values[i,j])
        if v>0: axH.text(j,i,f"{v:,}" if v<1000 else f"{v//1000}k",ha="center",va="center",fontsize=4.6,color="white" if np.log10(v+1)<np.log10(M.values.max()+1)*0.6 else "black")
axH.set_title("A  TE families of origin for strain-private piRNAs (top 16 families × 16 strains; cell = # private piRNAs, colour = log)",fontsize=9.6,fontweight="bold",loc="left")
cb=fig.colorbar(im,ax=axH,fraction=0.025,pad=0.01); cb.set_label("log10(# private piRNAs +1)",fontsize=7); cb.ax.tick_params(labelsize=6)
# Panel C — small-RNA (sRNA) TE-family expression (same families x strains)
_Elog=np.log10(E.values+1); _cmapE=__import__("seaborn").color_palette("mako_r",as_cmap=True).copy(); _cmapE.set_bad("#f0f0f0")   # light (low) -> dark (high) expression; empty cells = light grey
imE=axE.imshow(np.ma.masked_where(E.values<=0,_Elog),aspect="auto",cmap=_cmapE,vmin=0,vmax=(_Elog.max() if _Elog.max()>0 else 1))
axE.set_xticks(range(len(CANON))); axE.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.5)
axE.tick_params(labelleft=False)
axE.set_xticks(np.arange(-.5,len(CANON),1),minor=True); axE.set_yticks(np.arange(-.5,len(E),1),minor=True)
axE.grid(which="minor",color="#e6e6e6",linewidth=0.5); axE.tick_params(which="minor",length=0)
for i in range(len(E)):
    for j in range(len(CANON)):
        v=E.values[i,j]
        if v>0:
            lab=f"{v/1e6:.1f}M" if v>=1e6 else (f"{v/1e3:.0f}k" if v>=1e3 else f"{v:.0f}")
            axE.text(j,i,lab,ha="center",va="center",fontsize=4.0,color="white" if np.log10(v+1)<_Elog.max()*0.6 else "black")
axE.set_title("C  TE-family SMALL-RNA expression (sRNA reads on TEs; same families × strains; cell = summed sRNA counts, colour = log)",fontsize=9.6,fontweight="bold",loc="left")
cbE=fig.colorbar(imE,ax=axE,fraction=0.025,pad=0.01); cbE.set_label("log10(sRNA-on-TE +1)",fontsize=7); cbE.ax.tick_params(labelsize=6)
if not (E["NZO_HlLtJ"]>0).any(): axE.text(CANON.index("NZO_HlLtJ"),len(E)-0.2,"NZO RNA-seq\npending",ha="center",va="top",fontsize=4.2,color="#999",style="italic")
x=np.arange(len(CANON)); axB.bar(x,[tefrac[s] for s in CANON],color="#6a3d9a",edgecolor="white",linewidth=0.4)
for xi,s in zip(x,CANON):
    if pd.notna(tefrac[s]): axB.text(xi,tefrac[s]+0.5,f"{tefrac[s]:.0f}%",ha="center",va="bottom",fontsize=5.4,color="#6a3d9a")
axB.set_xticks(x); axB.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.5)
axB.set_ylabel("% TE-derived\n(lower bound)",fontsize=8); axB.set_title("B  TE-derived fraction of strain-private piRNAs per strain (mappable loci; lower bound — unplaced contigs excluded)",fontsize=9,fontweight="bold",loc="left")
axB.spines[['top','right']].set_visible(False)
fig.suptitle("TE families driving strain-private piRNAs across 16 strains — piRNA origin (A) vs small-RNA TE expression (C); LTR/ERVK(IAP) + LINE/L1 dominant, more so in wild-derived strains",fontsize=10.2,fontweight="bold",y=1.0)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_TE_private_families16.{e}",bbox_inches="tight")
print("wrote Fig_TE_private_families16.{png,pdf,svg} + source data")
