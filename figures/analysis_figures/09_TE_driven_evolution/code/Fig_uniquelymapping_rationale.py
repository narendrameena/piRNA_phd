#!/usr/bin/env python3
"""THEME 09 — WHY the TE-insertion test uses UNIQUELY-MAPPING (NH==1) piRNAs. Biological rationale, from data.
A strain-private piRNA SEQUENCE is born from a strain-private TE INSERTION (a new copy), but the TE family has many
near-identical ANCESTRAL copies genome-wide, so the piRNA read maps to the private insertion AND the ancestral copies.
The coordinate overlap test therefore cannot localize the PRODUCTION locus of a multi-mapping (NH>1) piRNA.
 A: per NH>1 strain-private candidate, its mapped loci split into private-insertion vs ancestral (a ~50/50 mix).
 B: only ~13-20% of NH>1 candidates are 'all-private' (localizable); most map to >=1 ancestral copy (ambiguous).
 C: consequence for fold-enrichment — the CLEAN, un-inflatable signal lives in NH==1 (uniquely-mapping); NH>1 collapses
    to ~the common-class 'maps-everywhere' artifact once the null is matched to mapping multiplicity.
Recomputed live from cand_self16 BAMs + ins16 private-insertion BEDs (klass5 strain-private vs expressed-elsewhere control)."""
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, pysam, subprocess, tempfile, os, numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
STRAINS=["CAST_EiJ","SPRET_EiJ","PWK_PhJ","WSB_EiJ"]   # the 4 wild-derived strains — where TE multi-mapping matters most
CLASSES={"unique: strain-private locus":"strain-private","expressed elsewhere (exact)":"common (control)"}
PRIV="#7a3b9a"; ANC="#c9b3d6"; CLEAN="#1b7837"; AMB="#d9a441"; NOPRIV="#bbbbbb"
WCOL={"CAST_EiJ":"#c51b8a","SPRET_EiJ":"#7a0177","PWK_PhJ":"#2c7fb8","WSB_EiJ":"#e6820a"}

def per_candidate(X, klass):
    """return dict id -> (n_total_loci, n_private_loci) for candidates of `klass` in strain X."""
    d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz"); d=d[d.strain==X]
    d=d.assign(id=X+"|"+d.timepoint.astype(str)+"|"+d.sequence)
    want=set(d.loc[d.klass5==klass,"id"])
    bam=pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam","rb")
    cb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG)
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.query_name not in want: continue
        cb.write(f"{a.reference_name}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
    cb.close(); bam.close()
    r=subprocess.run(f"sort -k1,1 -k2,2n {cb.name} | {BT} intersect -a - -b {PG}/ins16/{X}.ins_loci.bed -c",
                     shell=True,capture_output=True,text=True).stdout
    os.unlink(cb.name)
    tot=defaultdict(int); priv=defaultdict(int)
    for l in r.splitlines():
        f=l.split("\t"); q=f[3]; tot[q]+=1; priv[q]+=(1 if int(f[-1])>0 else 0)
    return {q:(tot[q],priv[q]) for q in tot}

def null_p(X):
    mrg=subprocess.run(f"sort -k1,1 -k2,2n {PG}/ins16/{X}.ins_loci.bed | {BT} merge",shell=True,capture_output=True,text=True).stdout
    insbp=sum(int(q[2])-int(q[1]) for q in (l.split('\t') for l in mrg.splitlines()) if len(q)>=3)
    gsize=sum(int(l) for c,l in (ln.split() for ln in open(f"{ROOT}/results/indexs/{X}/chrNameLength.txt")))
    return insbp/gsize

DATA={}
for X in STRAINS:
    DATA[X]={"p":null_p(X)}
    for kl,lab in CLASSES.items(): DATA[X][lab]=per_candidate(X,kl)

fig,(axA,axB,axC)=plt.subplots(1,3,figsize=(17,5.9),dpi=300,gridspec_kw=dict(width_ratios=[1.1,1.05,1.3],wspace=0.30))

# ---- A: fraction of a candidate's loci that are private (NH>1 strain-private) — violin per strain ----
fracs=[]
for X in STRAINS:
    sp=DATA[X]["strain-private"]; multi=[(kt,kp) for kt,kp in sp.values() if kt>1]
    fracs.append(np.array([kp/kt for kt,kp in multi]))
parts=axA.violinplot(fracs,positions=np.arange(len(STRAINS)),widths=0.82,showmedians=True,showextrema=False)
for i,pcb in enumerate(parts['bodies']): pcb.set_facecolor(WCOL[STRAINS[i]]); pcb.set_alpha(0.62); pcb.set_edgecolor("#444"); pcb.set_linewidth(0.6)
parts['cmedians'].set_color("#222"); parts['cmedians'].set_linewidth(1.5)
for i,fr in enumerate(fracs): axA.text(i,np.median(fr)+0.035,f"med {np.median(fr):.2f}",ha="center",fontsize=6.8,fontweight="bold",color="#222")
axA.axhline(1.0,color=CLEAN,lw=1.3,ls="--"); axA.text(len(STRAINS)-0.5,1.03,"all copies private = localizable",ha="right",va="bottom",fontsize=6.9,color=CLEAN)
axA.axhline(0.0,color="#999",lw=0.7,ls=":")
axA.set_xticks(np.arange(len(STRAINS))); axA.set_xticklabels([f"{s.replace('_','/')}\nn={len(fr):,}" for s,fr in zip(STRAINS,fracs)],fontsize=7.4)
axA.set_ylim(-0.10,1.16); axA.set_ylabel("fraction of a piRNA's mapped loci\nINSIDE a strain-private insertion",fontsize=8.6)
axA.set_title("A  Multi-mapping strain-private piRNAs map to a\nMIX of private + ancestral copies (median ~0.5)",fontsize=9.4,fontweight="bold",loc="left")
axA.spines[["top","right"]].set_visible(False)

# ---- B: localizability categories among NH>1 strain-private ----
cats=["all-private (localizable)","private + ancestral (ambiguous)","no private locus (elsewhere)"]; cc=[CLEAN,AMB,NOPRIV]
xb=np.arange(len(STRAINS)); bw=0.6; bottom=np.zeros(len(STRAINS)); vals={c:[] for c in cats}
for X in STRAINS:
    sp=DATA[X]["strain-private"]; multi=[(kt,kp) for kt,kp in sp.values() if kt>1]; n=len(multi)
    allp=sum(1 for kt,kp in multi if kp>0 and kp==kt); mix=sum(1 for kt,kp in multi if kp>0 and kp<kt); nop=sum(1 for kt,kp in multi if kp==0)
    vals[cats[0]].append(100*allp/n); vals[cats[1]].append(100*mix/n); vals[cats[2]].append(100*nop/n)
for c,col in zip(cats,cc):
    axB.bar(xb,vals[c],bw,bottom=bottom,color=col,edgecolor="white",lw=0.6,label=c)
    for xi,(v,b) in enumerate(zip(vals[c],bottom)):
        if v>4: axB.text(xi,b+v/2,f"{v:.0f}%",ha="center",va="center",fontsize=8,fontweight="bold",color="#222" if col!=NOPRIV else "#555")
    bottom+=np.array(vals[c])
axB.set_xticks(xb); axB.set_xticklabels([s.replace("_","/") for s in STRAINS],fontsize=9)
axB.set_ylabel("% of NH>1 strain-private piRNAs",fontsize=9); axB.set_ylim(0,100)
bpatch=[Patch(facecolor=col,label=c) for c,col in zip(cats,cc)]   # rendered as a figure-level legend below (avoids caption overlap)
axB.set_title("B  Only ~1/6 are cleanly localizable;\nmost map to >=1 ancestral copy",fontsize=9.6,fontweight="bold",loc="left")
axB.spines[["top","right"]].set_visible(False)

# ---- C: fold-enrichment decomposition — clean signal is in NH==1 ----
def folds(X):
    p=DATA[X]["p"]; out={}
    for lab in ["strain-private","common (control)"]:
        v=DATA[X][lab]
        for grp,sel in [("NH==1",lambda kt:kt==1),("NH>1",lambda kt:kt>1)]:
            sub=[(kt,kp) for kt,kp in v.values() if sel(kt)]
            if not sub: out[(lab,grp)]=np.nan; continue
            obs=np.mean([1 if kp>0 else 0 for kt,kp in sub])
            k=np.array([kt for kt,kp in sub]); expm=(1-(1-p)**k).mean()
            out[(lab,grp)]=obs/expm   # matched null (=single null for NH==1 since k=1)
    return out
F={X:folds(X) for X in STRAINS}
groups=[("strain-private","NH==1",CLEAN,"strain-private · NH==1 (clean)"),
        ("strain-private","NH>1",AMB,"strain-private · NH>1 (matched)"),
        ("common (control)","NH==1",NOPRIV,"common · NH==1"),
        ("common (control)","NH>1","#7f7f7f","common · NH>1 (matched)")]
xc=np.arange(len(STRAINS)); bw=0.2
for i,(lab,grp,col,leg) in enumerate(groups):
    vals=[F[X][(lab,grp)] for X in STRAINS]
    b=axC.bar(xc+(i-1.5)*bw,vals,bw,color=col,edgecolor="white",lw=0.4,label=leg)
    for xi,v in zip(xc+(i-1.5)*bw,vals):
        if np.isfinite(v): axC.text(xi,v+0.3,f"{v:.1f}×",ha="center",va="bottom",fontsize=6.6,fontweight="bold",color=col if col!=NOPRIV else "#555")
axC.axhline(1,color="#555",ls=":",lw=1); axC.text(len(STRAINS)-0.5,1.25,"chance (1×)",ha="right",fontsize=7,color="#555")
axC.set_xticks(xc); axC.set_xticklabels([s.replace("_","/") for s in STRAINS],fontsize=9)
axC.set_ylabel("fold-enrichment at private insertions\n(obs ÷ multiplicity-matched null)",fontsize=9)
_mx=np.nanmax([F[X][(g[0],g[1])] for X in STRAINS for g in groups]); axC.set_ylim(0,_mx*1.42)
axC.legend(fontsize=6.5,frameon=False,loc="upper center",ncol=2,columnspacing=1.1,handlelength=1.3,handletextpad=0.4); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  Clean signal is UNIQUELY in NH==1; NH>1 strain-private\ncollapses to the common-class 'maps-everywhere' artifact",fontsize=9.6,fontweight="bold",loc="left")

fig.suptitle("Why the TE-insertion test uses uniquely-mapping (NH==1) piRNAs — the discarded NH>1 fraction is un-localizable, not lost biology",fontsize=12,fontweight="bold",y=0.985)
fig.subplots_adjust(left=0.045,right=0.985,top=0.84,bottom=0.26,wspace=0.28)
fig.legend(handles=bpatch,loc="center",bbox_to_anchor=(0.5,0.145),ncol=3,fontsize=7.4,frameon=False,columnspacing=2.0,handlelength=1.4,handletextpad=0.5,title="Panel B — localizability of NH>1 strain-private piRNAs",title_fontsize=7.2)
fig.text(0.5,0.045,"A strain-private piRNA is born from a strain-private TE INSERTION, but the family's many near-identical ANCESTRAL copies mean a multi-mapping (NH>1) read hits the private insertion AND ancestral copies "
  "(A: median ~0.5 of loci private; B: most map to ≥1 ancestral copy → production locus AMBIGUOUS). NH==1 keeps only the unambiguous single-locus piRNAs — no localizable TE biology is lost; C: strain-private NH==1 enrichment has no "
  "common-control counterpart (~1×), while matched-null NH>1 strain-private ≈ the common-class 'maps-everywhere' artifact.",
  ha="center",va="center",fontsize=6.3,color="#555",wrap=True)
_SD=f"{ROOT}/figures/analysis_figures/09_TE_driven_evolution/data/source_data"; os.makedirs(_SD,exist_ok=True)
rows=[]
for X in STRAINS:
    for lab in ["strain-private","common (control)"]:
        for kt,kp in DATA[X][lab].values(): rows.append(dict(strain=X,klass=lab,n_loci=kt,n_private_loci=kp,frac_private=round(kp/kt,3),NH_group="NH==1" if kt==1 else "NH>1"))
pd.DataFrame(rows).to_csv(f"{_SD}/SourceData_uniquelymapping_rationale.csv",index=False)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_uniquelymapping_rationale.{e}",bbox_inches="tight")
print("wrote Fig_uniquelymapping_rationale.{png,pdf,svg} + source data")
for X in STRAINS:
    print(f"  {X}: strain-private NH==1 fold={F[X][('strain-private','NH==1')]:.1f}x  NH>1 matched={F[X][('strain-private','NH>1')]:.1f}x  |  common NH>1 matched={F[X][('common (control)','NH>1')]:.1f}x")
