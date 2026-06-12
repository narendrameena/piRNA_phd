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
d=pd.read_csv(f"{U}/unique16/final_classified.csv.gz")
fammat=defaultdict(dict); tefrac={}; nsp={}
for X in CANON:
    g=d[(d.strain==X)&(d.klass=="unique: strain-private locus")].copy(); g["id"]=X+"|"+g.timepoint+"|"+g.sequence
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
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,(axH,axB)=plt.subplots(2,1,figsize=(12,8.4),dpi=300,gridspec_kw={"height_ratios":[3,1],"hspace":0.32})
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
cb=fig.colorbar(im,ax=axH,fraction=0.025,pad=0.01); cb.set_label("log10(count+1)",fontsize=7); cb.ax.tick_params(labelsize=6)
x=np.arange(len(CANON)); axB.bar(x,[tefrac[s] for s in CANON],color="#6a3d9a",edgecolor="white",linewidth=0.4)
for xi,s in zip(x,CANON):
    if pd.notna(tefrac[s]): axB.text(xi,tefrac[s]+0.5,f"{tefrac[s]:.0f}%",ha="center",va="bottom",fontsize=5.4,color="#6a3d9a")
axB.set_xticks(x); axB.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=7.5)
axB.set_ylabel("% TE-derived\n(lower bound)",fontsize=8); axB.set_title("B  TE-derived fraction of strain-private piRNAs per strain (mappable loci; lower bound — unplaced contigs excluded)",fontsize=9,fontweight="bold",loc="left")
axB.spines[['top','right']].set_visible(False)
fig.suptitle("TE families driving strain-private piRNAs across 16 strains — LTR/ERVK(IAP) + LINE/L1 dominant, more so in the wild-derived strains",fontsize=10.4,fontweight="bold",y=1.0)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_TE_private_families16.{e}",bbox_inches="tight")
print("wrote Fig_TE_private_families16.{png,pdf,svg} + source data")
