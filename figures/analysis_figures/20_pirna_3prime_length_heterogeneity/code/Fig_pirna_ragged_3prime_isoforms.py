#!/usr/bin/env python3
"""THEME 20 — piRNA 3' length heterogeneity: imprecise 3'-end trimming gives RAGGED 3' ends (a defined 5' end but
VARIABLE 3' length), so a shorter read is the EXACT 5' prefix of a longer piRNA isoform. piRNAs are NOT cut at a
single exact length; they spread around the per-tp stage-characteristic window (E16.5 27 nt pre-pachytene -> P20.5
30 nt pachytene), and the raggedness INCREASES across development, peaking in pachytene 30-nt piRNAs.
Biology: length is set by the bound PIWI footprint (pre-pachytene MILI ~26-27 nt; pachytene MIWI ~29-30 nt) plus
3'->5' trimming by PNLDC1, terminated by HENMT1 2'-O-methylation -- trimming is not single-nt precise -> ragged 3'.
The E16.5 pre-pachytene pool ALSO includes the fetal nuclear PIWI MIWI2 (~28 nt; de-novo TE methylation; MILI/MIWI2
ping-pong, see drivers panel C). NOTE: MILI vs MIWI2 are NOT separable by read length (overlapping footprints + 3'
raggedness), so the windows here are READ-LENGTH classes, not protein assignments -- 28 nt is acknowledged, not analysed.
Data: unique_pirna/unique16/final_classified_clean_2read.csv.gz (>=2-read, 25-32 nt, within-tp classified piRNAs).
Window (per-tp, data-verified): E16.5 -> 27; P12.5 -> 27 & 30; P20.5 -> 30 (make_stage_peak_unique.py). 24 nt excluded
by the 1U-peak test (piRNA range 25-32 nt)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/20_pirna_3prime_length_heterogeneity"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#B2182B"}
# windows are READ-LENGTH classes (the stage-characteristic peak); MILI vs MIWI2 are NOT separable by length
# (overlapping footprints + 3' raggedness), so MIWI2 (~28 nt, fetal) is acknowledged but not analysed as its own window
WIN={"16.5dpc":[27],"12.5dpp":[27,30],"20.5dpp":[30]}; LMIN,LMAX=25,32
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","timepoint"]); d["L"]=d.sequence.str.len()
# per-tp distinct-sequence sets by length
SETS={tp:{L:set(d[(d.timepoint==tp)&(d.L==L)].sequence) for L in range(23,35)} for tp in TPS}
def ragged(tp,L):
    s=SETS[tp]; sset=s[L]; n=len(sset) or 1; off={}
    for k in (1,2,3): off[-k]=np.mean([x[:L-k] in s[L-k] for x in sset])
    pre={k:{x[:L] for x in s[L+k]} for k in (1,2,3)}
    for k in (1,2,3): off[k]=np.mean([x in pre[k] for x in sset])
    anyiso=np.mean([any((x[:L-k] in s[L-k]) for k in (1,2,3)) or any((x in pre[k]) for k in (1,2,3)) for x in sset])
    return n,off,anyiso
ENTRIES=[(tp,L) for tp in TPS for L in WIN[tp]]; rag={(tp,L):ragged(tp,L) for tp,L in ENTRIES}
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(13,9.6),dpi=300)
# A: per-tp length distribution (developmental shift + broad/imprecise)
PROF={}
for tp in TPS:
    prof=np.array([len(SETS[tp][L]) for L in range(LMIN,LMAX+1)],float); prof/=prof.sum(); PROF[tp]=100*prof
    axA.plot(range(LMIN,LMAX+1),PROF[tp],"-o",color=TPCOL[tp],ms=4,label=TPN[tp])
    for L in WIN[tp]: axA.scatter([L],[PROF[tp][L-LMIN]],s=130,facecolor="none",edgecolor=TPCOL[tp],lw=1.8,zorder=5)
top=max(PROF[tp].max() for tp in TPS)*1.17; axA.set_ylim(0,top)
axA.set_xlabel("piRNA length (nt)",fontsize=9.5); axA.set_ylabel("% of distinct piRNAs",fontsize=9.5); axA.set_xticks(range(LMIN,LMAX+1))
axA.legend(fontsize=8,frameon=False,title="timepoint",title_fontsize=8,loc="upper right"); axA.spines[["top","right"]].set_visible(False)
y27=max(PROF[tp][27-LMIN] for tp in TPS); y30=PROF["20.5dpp"][30-LMIN]
axA.annotate("27 nt pre-pachytene\n(MILI footprint)",xy=(27,y27),xytext=(25.6,top*0.99),ha="center",va="top",fontsize=6.8,color="#4393C3",fontweight="bold",arrowprops=dict(arrowstyle="->",color="#4393C3",lw=1,shrinkB=4))
axA.annotate("30 nt pachytene\n(MIWI footprint)",xy=(30,y30),xytext=(28.3,top*0.87),ha="center",va="top",fontsize=6.8,color="#B2182B",fontweight="bold",arrowprops=dict(arrowstyle="->",color="#B2182B",lw=1,shrinkB=4))
axA.set_title("A  piRNA length is NOT a single exact cut — a broad spread that\nSHIFTS 27 nt (pre-pachytene) -> 30 nt (pachytene) across development",fontsize=9.6,fontweight="bold",loc="left")
# B: any ragged-3' isoform fraction per tp x window (developmental increase)
labs=[f"{TPN[tp]}\n{L} nt" for tp,L in ENTRIES]; vals=[100*rag[e][2] for e in ENTRIES]; cols=[TPCOL[tp] for tp,_ in ENTRIES]
axB.bar(range(len(ENTRIES)),vals,color=cols,edgecolor="white",zorder=2)
_psB=pd.read_csv(f"{T}/data/SourceData_perstrain_ragged.csv")   # per-strain ragged% (16 strains) for the spread
for i,(tp,L) in enumerate(ENTRIES):
    pv=_psB[(_psB.tp==TPN[tp])&(_psB.L==L)].rag_pct.values
    axB.scatter(i+np.linspace(-0.17,0.17,len(pv)),pv,s=11,color="#333",alpha=0.55,edgecolor="white",lw=0.3,zorder=4)
    axB.text(i,max(vals[i],pv.max() if len(pv) else vals[i])+2,f"{vals[i]:.0f}%\n(n={rag[(tp,L)][0]:,})",ha="center",va="bottom",fontsize=7.1,fontweight="bold")
axB.set_xticks(range(len(ENTRIES))); axB.set_xticklabels(labs,fontsize=8.3); axB.set_ylabel("% with a ragged-3' isoform  (bar = pooled · dots = 16 strains)",fontsize=8.3)
axB.set_ylim(0,max(max(vals),_psB.rag_pct.max())*1.16); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  A shorter read is the EXACT 5' prefix of a longer piRNA\n(ragged 3' ends) — rising with development, peaking in pachytene",fontsize=9.6,fontweight="bold",loc="left")
# C: 3'-end offset MAGNITUDE profile (1/2/3 nt shorter <- | -> longer)
offs=[-3,-2,-1,1,2,3]
for tp,L in ENTRIES:
    _,off,_=rag[(tp,L)]; axC.plot(offs,[100*off[o] for o in offs],"-o",ms=4.5,lw=1.6,color=TPCOL[tp],alpha=0.55+0.30*(L==30),label=f"{TPN[tp]} · {L}nt")
axC.axvline(0,color="#aaa",lw=0.9,ls=":"); axC.set_xticks(offs); axC.set_xticklabels(["-3","-2","-1","+1","+2","+3"])
axC.set_xlabel("3' length offset of isoform (nt)   ·   shorter $\\leftarrow$  |  $\\rightarrow$  longer",fontsize=9.2)
axC.set_ylabel("% of window piRNAs with this isoform",fontsize=9.2); axC.set_ylim(0,axC.get_ylim()[1]*1.34); axC.legend(fontsize=6.8,frameon=False,ncol=1,loc="upper right"); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  How big is the cut imprecision? Mostly +/-1 nt, decaying to +/-3 nt\n(defined 5' end, single-nt-ragged 3' end)",fontsize=9.6,fontweight="bold",loc="left")
# D: developmental TRAJECTORY of raggedness — two protein-class tracks (MILI 27nt, MIWI 30nt); P12.5 = handoff (both)
e27=("16.5dpc",27); e30=("20.5dpp",30); p27=("12.5dpp",27); p30=("12.5dpp",30); TX={"16.5dpc":0,"12.5dpp":1,"20.5dpp":2}
mili=[e27,p27]; miwi=[p30,e30]
axD.plot([TX[t] for t,_ in mili],[100*rag[e][2] for e in mili],"-o",color="#2166ac",lw=2.6,ms=12,label="MILI · 27 nt (pre-pachytene)",zorder=3)
axD.plot([TX[t] for t,_ in miwi],[100*rag[e][2] for e in miwi],"-o",color="#b2182b",lw=2.6,ms=12,label="MIWI · 30 nt (pachytene)",zorder=3)
LO={e27:(0,3.4,"center"),p27:(-0.14,-3.2,"right"),p30:(-0.14,3.0,"right"),e30:(0,3.4,"center")}
for e in mili+miwi:
    dx,dy,ha=LO[e]; axD.text(TX[e[0]]+dx,100*rag[e][2]+dy,f"{100*rag[e][2]:.0f}%",ha=ha,va="center",fontsize=9.5,fontweight="bold",color=("#2166ac" if e[1]==27 else "#b2182b"))
axD.annotate("",xy=(1,100*rag[p30][2]-1.5),xytext=(1,100*rag[p27][2]+1.5),arrowprops=dict(arrowstyle="-|>",color="#999",lw=1.7))
axD.text(1.08,(100*rag[p27][2]+100*rag[p30][2])/2,"MILI->MIWI\nhandoff",ha="left",va="center",fontsize=6.9,color="#666",style="italic")
axD.set_xticks([0,1,2]); axD.set_xticklabels(["E16.5","P12.5","P20.5"],fontsize=9.5); axD.set_xlim(-0.4,2.5)
axD.set_ylabel("% ragged-3' piRNAs",fontsize=9.5); axD.set_ylim(30,100*rag[e30][2]*1.3); axD.spines[["top","right"]].set_visible(False)
axD.legend(fontsize=7.3,frameon=False,loc="upper left")
axD.text(0.035,0.56,"E16.5 pre-pachytene also loads MIWI2\n(~28 nt, fetal nuclear PIWI, de novo TE\nmethylation) — not separable from MILI\nby read length alone",transform=axD.transAxes,fontsize=6.2,color="#2166ac",va="top",style="italic")
axD.text(0.97,0.05,f"raggedness climbs along BOTH tracks and jumps at the\nMILI->MIWI switch ({100*rag[e27][2]:.0f}% -> {100*rag[e30][2]:.0f}%, {rag[e30][2]/rag[e27][2]:.1f}x overall)",transform=axD.transAxes,ha="right",va="bottom",fontsize=7,color="#444",style="italic")
axD.set_title("D  Biology: developmental MILI(27)->MIWI(30) handoff —\npachytene piRNAs MORE ragged; P12.5 carries BOTH classes",fontsize=9.6,fontweight="bold",loc="left")
fig.suptitle("piRNA 3' ends are RAGGED, not cut at an exact length — imprecise PNLDC1 trimming around the per-stage window (27 nt pre-pachytene -> 30 nt pachytene)",fontsize=11.5,fontweight="bold",y=1.0)
fig.text(0.5,0.012,"3'-isoform = a shorter read whose sequence is the EXACT 5' prefix of a longer read (no mismatch / no SNP); genomically validated — 99-100% map to the SAME 5' locus (one molecule, ragged 3' trim)",ha="center",fontsize=7,color="#555",style="italic")
fig.tight_layout(rect=[0,0.028,1,0.96])
out=f"{T}/figures/Fig_pirna_ragged_3prime_isoforms"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
src=[]
for tp,L in ENTRIES:
    n,off,anyiso=rag[(tp,L)]; src.append(dict(timepoint=TPN[tp],window_len=L,n=n,ragged_any_pct=round(100*anyiso,1),**{f"iso_{o:+d}nt_pct":round(100*off[o],1) for o in offs}))
pd.DataFrame(src).to_csv(f"{T}/data/source_data/SourceData_Fig_pirna_ragged_3prime_isoforms.csv",index=False)
print("wrote",out)
for tp,L in ENTRIES: print(f"  {TPN[tp]} {L}nt: n={rag[(tp,L)][0]:,} ragged-any={100*rag[(tp,L)][2]:.0f}%")
