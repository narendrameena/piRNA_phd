#!/usr/bin/env python3
"""THEME 20 - Fig 2: BIOLOGICAL DRIVERS of piRNA 3' raggedness. Ragged 3' ends (defined 5', variable 3') are NOT
noise: (A) they dominate the EXPRESSED pool (read-weighted >> sequence-weighted; ragged isoforms more abundant);
(B) universal across all 16 strains (developmental rise 8/8; wild = classical); (C) enriched in SECONDARY (ping-pong)
piRNAs at E16.5 (MILI/MIWI2, 10-nt overlap z=7.2); (D) TE<->genic origin FLIPS across the MILI->MIWI transition;
(E) ~90% come from pachytene clusters, the dominant/high-expression clusters being the most ragged. (F) Model,
BioMNI-verified (3/3): MIWI binds ~30nt, PNLDC1 trims, HENMT1 caps, A-MYB drives; Pnldc1-KO -> longer untrimmed.
Data: SourceData_{expression,perstrain,pingpong_e16,cluster_origin}.csv (+ TE/orientation recomputed from
clean_2read x sense_antisense). Deterministic (no randomness)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/20_pirna_3prime_length_heterogeneity"; DD=f"{T}/data"
TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#B2182B"}
ENT=[("16.5dpc",27),("12.5dpp",27),("12.5dpp",30),("20.5dpp",30)]; LAB=[f"{TPN[t]}\n{l} nt" for t,l in ENT]
fig,axs=plt.subplots(2,3,figsize=(16.5,9.7),dpi=300); (axA,axB,axC),(axD,axE,axF)=axs

# ---- A: expression — ragged isoforms dominate the read pool ----
ex=pd.read_csv(f"{DD}/SourceData_expression_ragged.csv"); ex=ex.set_index(["tp","L"]).loc[ENT].reset_index()
x=np.arange(len(ex)); w=0.38
axA.bar(x-w/2,ex.rag_pct_seq,w,color="#c4c4c4",edgecolor="white",label="by distinct sequence")
axA.bar(x+w/2,ex.rag_pct_reads,w,color=[TPCOL[t] for t in ex.tp],edgecolor="white",label="by READ abundance")
for i,r in ex.iterrows():
    axA.text(i-w/2,r.rag_pct_seq+1,f"{r.rag_pct_seq:.0f}",ha="center",va="bottom",fontsize=7,color="#777")
    axA.text(i+w/2,r.rag_pct_reads+1,f"{r.rag_pct_reads:.0f}",ha="center",va="bottom",fontsize=7.6,fontweight="bold")
axA.set_xticks(x); axA.set_xticklabels(LAB,fontsize=8); axA.set_ylabel("% ragged-3'",fontsize=9.5); axA.set_ylim(0,100)
axA.legend(fontsize=7.4,frameon=False,loc="upper left"); axA.spines[["top","right"]].set_visible(False)
axA.text(0.035,0.80,f"ragged piRNAs are MORE abundant\n(median {ex.med_ab_rag.iloc[-1]:.0f} vs {ex.med_ab_norag.iloc[-1]:.0f} reads, P20.5)",transform=axA.transAxes,ha="left",va="top",fontsize=7.2,color="#444",style="italic")
axA.set_title("A  Ragged isoforms DOMINATE the expressed pool\n(read-weighted >> sequence-weighted) - not low-count noise",fontsize=9.3,fontweight="bold",loc="left")

# ---- B: per-strain robustness ----
SHORTCOL={"E16.5":"#4393C3","P12.5":"#E8852B","P20.5":"#B2182B"}
ps=pd.read_csv(f"{DD}/SourceData_perstrain_ragged.csv"); ps["key"]=[f"{t}\n{l} nt" for t,l in zip(ps.tp,ps.L)]
for i,k in enumerate(LAB):
    sub=ps[ps.key==k].sort_values("rag_pct"); tp=sub.tp.iloc[0]; n=len(sub)
    xj=i+np.linspace(-0.16,0.16,n)
    axB.scatter(xj,sub.rag_pct,s=27,color=["#111" if w else SHORTCOL[tp] for w in sub.wild],edgecolor="white",lw=0.4,zorder=3)
    axB.plot([i-0.26,i+0.26],[sub.rag_pct.mean()]*2,color="#000",lw=2,zorder=4)
axB.set_xticks(range(4)); axB.set_xticklabels(LAB,fontsize=8); axB.set_ylabel("% ragged-3'  (one point = one strain)",fontsize=9.2)
axB.set_ylim(20,90); axB.set_xlim(-0.5,3.5); axB.spines[["top","right"]].set_visible(False)
axB.text(0.03,0.95,"every strain: E16.5 < P20.5  (8/8)\nwild (black) ≈ classical",transform=axB.transAxes,ha="left",va="top",fontsize=7.4,color="#444",style="italic")
axB.set_title("B  Universal across all 16 strains (bar = mean)\n- a conserved property of piRNA biogenesis",fontsize=9.3,fontweight="bold",loc="left")

# ---- C: E16.5 ping-pong (secondary) ----
pp=pd.read_csv(f"{DD}/SourceData_pingpong_e16.csv"); m=pp.groupby("pp").rag.mean()*100
vals=[m.get(False,np.nan),m.get(True,np.nan)]
axC.bar([0,1],vals,color=["#9ecae1","#08519c"],width=0.6,edgecolor="white")
for i,v in enumerate(vals): axC.text(i,v+1.5,f"{v:.0f}%",ha="center",va="bottom",fontsize=12,fontweight="bold")
axC.set_xticks([0,1]); axC.set_xticklabels(["primary\n(non-ping-pong)","secondary\n(ping-pong)"],fontsize=8.6)
axC.set_ylabel("% ragged-3'  (E16.5, 27 nt)",fontsize=9.5); axC.set_ylim(0,85); axC.spines[["top","right"]].set_visible(False)
axC.text(0.5,0.95,"10-nt 5'-overlap signature  z = 7.2  (p=3e-22)\nMILI (PIWIL2)  <->  MIWI2 (PIWIL4)",transform=axC.transAxes,ha="center",va="top",fontsize=7.5,color="#444",style="italic")
axC.set_title("C  Secondary (ping-pong) piRNAs are far more ragged\nat E16.5 - slicer-set 5', trimmed (ragged) 3'",fontsize=9.3,fontweight="bold",loc="left")

# ---- D: TE vs genic flip (recompute) ----
clean=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint"])
clean["L"]=clean.sequence.str.len(); clean["cand_id"]=clean.strain+"|"+clean.timepoint+"|"+clean.sequence
oc=pd.read_csv(f"{U}/sense_antisense/SourceData_sense_antisense16_percand.csv.gz"); TEm=dict(zip(oc.id,oc.family.notna()))
teres=[]
for tp,L in ENT:
    sub=clean[clean.timepoint==tp]; S={k:set(sub[sub.L==k].sequence) for k in range(23,35)}
    sset=S[L]; pre={k:{y[:L] for y in S[L+k]} for k in (1,2,3)}
    israg=lambda x,S=S,pre=pre,L=L: any(x[:L-k] in S[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3))
    seqte={}
    for r in sub[sub.L==L].itertuples(): seqte[r.sequence]=seqte.get(r.sequence,False) or TEm.get(r.cand_id,False)
    df=pd.DataFrame([(israg(s),seqte.get(s,False)) for s in sset],columns=["rag","te"]); g=df.groupby("te").rag.mean()*100
    teres.append((g.get(True,np.nan),g.get(False,np.nan)))
x=np.arange(len(ENT))
axD.bar(x-w/2,[t[0] for t in teres],w,color="#d6604d",edgecolor="white",label="TE-derived")
axD.bar(x+w/2,[t[1] for t in teres],w,color="#4393C3",edgecolor="white",label="genic / non-TE")
for i,t in enumerate(teres):
    axD.text(i-w/2,t[0]+1,f"{t[0]:.0f}",ha="center",va="bottom",fontsize=7,fontweight="bold",color="#b2182b")
    axD.text(i+w/2,t[1]+1,f"{t[1]:.0f}",ha="center",va="bottom",fontsize=7,fontweight="bold",color="#2166ac")
axD.set_xticks(x); axD.set_xticklabels(LAB,fontsize=8); axD.set_ylabel("% ragged-3'",fontsize=9.5); axD.set_ylim(0,85)
axD.legend(fontsize=7.4,frameon=False,loc="upper left"); axD.spines[["top","right"]].set_visible(False)
axD.text(0.035,0.80,"TE > genic at E16.5  ->  genic > TE at P20.5\n(across the MILI -> MIWI transition)",transform=axD.transAxes,ha="left",va="top",fontsize=7.3,color="#444",style="italic")
axD.set_title("D  Origin FLIPS with development: TE-silencing piRNAs\nragged early; genic/cluster piRNAs ragged in pachytene",fontsize=9.3,fontweight="bold",loc="left")

# ---- E: pachytene cluster origin ----
co=pd.read_csv(f"{DD}/SourceData_cluster_origin.csv"); mc=co.groupby("incluster").rag.mean()*100; frac=100*co.incluster.mean()
cl=co[co.incluster]; med=cl.fpm.median(); hi=cl[cl.fpm>=med].rag.mean()*100; lo=cl[cl.fpm<med].rag.mean()*100
bars=[(0,mc.get(False,np.nan),"#cccccc"),(1,mc.get(True,np.nan),"#B2182B"),(2.4,lo,"#fcae91"),(3.4,hi,"#a50f15")]
for xx,v,cc in bars: axE.bar(xx,v,width=0.6,color=cc,edgecolor="white"); axE.text(xx,v+1,f"{v:.0f}%",ha="center",va="bottom",fontsize=9.5,fontweight="bold")
axE.set_xticks([b[0] for b in bars]); axE.set_xticklabels(["non-\ncluster","pachytene\ncluster","low-expr\ncluster","high-expr\ncluster"],fontsize=7.8)
axE.set_ylabel("% ragged-3'  (P20.5, 30 nt)",fontsize=9.5); axE.set_ylim(0,85); axE.spines[["top","right"]].set_visible(False); axE.axvline(1.7,color="#e5e5e5",lw=0.9)
axE.text(0.5,0.95,f"{frac:.0f}% of pachytene piRNAs are cluster-derived;\nthe dominant (high-expression) clusters are raggedest",transform=axE.transAxes,ha="center",va="top",fontsize=7.3,color="#444",style="italic")
axE.set_title("E  ~90% come from pachytene clusters; the dominant\nclusters produce the MOST ragged piRNAs",fontsize=9.3,fontweight="bold",loc="left")

# ---- F: biology synthesis schematic ----
axF.axis("off"); axF.set_xlim(0,1); axF.set_ylim(0,1)
axF.text(0.5,0.965,"Model: defined 5' end + imprecise 3' trim",ha="center",fontsize=9.5,fontweight="bold")
y0=0.80
for i,(frac_len,col) in enumerate([(0.30,"#4393C3"),(0.35,"#4393C3"),(0.41,"#B2182B"),(0.45,"#B2182B")]):
    yy=y0-i*0.06; axF.plot([0.14,0.14+frac_len],[yy,yy],color=col,lw=5,solid_capstyle="butt")
    axF.plot([0.14+frac_len],[yy],marker="o",ms=3,color=col)
axF.plot([0.14,0.14],[y0+0.05,y0-0.22],color="#000",lw=1.3,ls=":")
axF.text(0.14,y0+0.085,"5' end FIXED\n(PIWI footprint, 1U)",ha="center",fontsize=7,color="#222")
axF.annotate("3' end RAGGED",xy=(0.59,y0-0.09),xytext=(0.62,y0-0.02),fontsize=7.4,color="#222",fontweight="bold",arrowprops=dict(arrowstyle="->",color="#666",lw=0.9))
axF.text(0.62,y0-0.20,"PNLDC1 trims 3'->5';\nHENMT1 2'-O-Me caps",ha="left",fontsize=6.8,color="#555")
axF.text(0.5,0.40,"MILI ~27 nt (E16.5, pre-pachytene)  ->  MIWI ~30 nt (P20.5, pachytene)",ha="center",fontsize=7.7,color="#222",fontweight="bold")
axF.text(0.5,0.31,"pachytene MORE ragged (68% vs 39%) = more 3' trimming imprecision",ha="center",fontsize=7.3,color="#555")
axF.text(0.5,0.205,"Pnldc1-KO -> LONGER untrimmed (31-35 nt) piRNAs",ha="center",fontsize=7.3,color="#555",style="italic")
axF.text(0.5,0.115,"BioMNI 3/3 verified: MIWI binds . PNLDC1 trims . HENMT1 caps . A-MYB drives",ha="center",fontsize=6.9,color="#777")
axF.text(0.5,0.03,"genomically validated: 99-100% of isoforms = SAME 5' locus, exact prefix, no SNP",ha="center",fontsize=6.7,color="#999")
axF.set_title("F  Biology synthesis (BioMNI-verified)",fontsize=9.3,fontweight="bold",loc="left")

fig.suptitle("Ragged 3' ends are a BIOLOGICAL SIGNATURE of piRNA processing - abundant, strain-universal, secondary/ping-pong-linked, and pachytene-cluster-derived",fontsize=11.6,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.965])
out=f"{T}/figures/Fig_pirna_ragged_3prime_drivers"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pd.DataFrame([dict(panel="C_pingpong",ragged_pp=round(m.get(True),1),ragged_nonpp=round(m.get(False),1)),
             dict(panel="E_cluster",cluster_frac=round(frac,1),ragged_cluster=round(mc.get(True),1),ragged_noncluster=round(mc.get(False),1),ragged_hi_expr=round(hi,1),ragged_lo_expr=round(lo,1))]).to_csv(f"{DD}/SourceData_Fig_drivers_summary.csv",index=False)
print("wrote",out); print("ping-pong:",round(m.get(True),1),"vs",round(m.get(False),1)); print("cluster:",round(mc.get(True),1),"vs",round(mc.get(False),1),"| hi/lo expr:",round(hi,1),round(lo,1)); print("TE flip:",[(round(a,1),round(b,1)) for a,b in teres])
