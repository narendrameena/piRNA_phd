#!/usr/bin/env python3
"""Comprehensive locus figure: EXAMPLE + ANATOMY + NUCLEOTIDE RESOLUTION (real read depth).
(1) cross-strain status strip (strain-private); (2) TE content; (3) strand-resolved piRNA coverage
(antisense-to-TE=silencing red, sense grey) + 1U inset; (4) a REAL piRNA pair at base resolution
(ping-pong: sense+antisense 10-nt 5' overlap with 1U/10A; or the top piRNA if no pair).
Usage: make_locus_full.py <tag> <chrom> <start> <end> <te_label> <te_strand>"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER
import pav_clusters as pc   # shared styling: pbadge panel headers + semantic colour palette + genes_at gene track
import numpy as np, pandas as pd, pysam
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
RB="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise"
tag,CHROM,S,E,TELAB,TEST=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],sys.argv[6]
import glob as _glob, os as _os
STRAIN=sys.argv[7] if len(sys.argv)>7 else "SPRET_EiJ"                              # strain to render (default SPRET for back-compat)
BAMS=sorted(_glob.glob(f"{RB}/{STRAIN}/{STRAIN}-*/Aligned.sortedByCoord.out.bam"))  # POOL all tp/rep so embryonic (16.5dpc) source loci are not missed
N=E-S; comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))
nb=200; fwd=np.zeros(nb); rev=np.zeros(nb); first=[]
sense=defaultdict(Counter); anti=defaultdict(Counter)
for BAM in BAMS:
    bam=pysam.AlignmentFile(BAM,"rb")
    for a in bam.fetch(CHROM,S,E):
        if a.is_unmapped: continue
        L=a.reference_end-a.reference_start
        if not(24<=L<=32) or not a.query_sequence: continue
        b0=int((a.reference_start-S)/N*nb); b1=int((a.reference_end-S)/N*nb)
        for b in range(max(0,b0),min(nb,b1+1)): (rev if a.is_reverse else fwd)[b]+=1
        if a.is_reverse: anti[a.reference_end-1][rc(a.query_sequence)]+=1; first.append(comp.get(a.query_sequence[-1],"N"))
        else: sense[a.reference_start][a.query_sequence]+=1; first.append(a.query_sequence[0])
    bam.close()
fc=pd.Series(first).value_counts(normalize=True)*100
antiC,senseC=(fwd,rev) if TEST=="-" else (rev,fwd)
# find a real ping-pong pair (sense@i, antisense@i+9), prefer 1U(anti)+10A(sense)
pair=None
for i in sense:
    if i+9 in anti:
        ss,sc=sense[i].most_common(1)[0]; aq,acn=anti[i+9].most_common(1)[0]
        if len(ss)>=10 and len(aq)>=10:
            sgn=(sc+acn)*(3 if(aq[0]=="T" and ss[9]=="A") else 1)
            if pair is None or sgn>pair[0]: pair=(sgn,i,ss,sc,aq,acn)
_tebed=f"{U}/sense_antisense/{STRAIN}.TE_stranded.bed"
if _os.path.exists(_tebed):
    TE=pd.read_csv(_tebed,sep="\t",header=None,names=["c","s","e","fam","sc","st"]); TE=TE[(TE.c==CHROM)&(TE.s<E)&(TE.e>S)]
else:                                                                              # any strain w/o a bed: RepeatMasker .out via te_at (own coords, PanSN naming)
    TE=pd.DataFrame(pc.te_at(STRAIN,CHROM.split("chr")[-1],S,E),columns=["s","e","st","fam"]); TE["c"]=CHROM
TECOL={"LINE/L1":"#E69F00","LTR/ERVK":"#6a3d9a","LTR/ERVL-MaLR":"#b15928","LTR/ERV1":"#cab2d6","SINE/Alu":"#33a02c","SINE/B2":"#1f78b4","LTR/ERVL":"#a6cee3"}
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(11,8.8),dpi=300)
gs=fig.add_gridspec(4,4,height_ratios=[0.95,0.8,1.9,1.5],hspace=0.58,wspace=0.32)
axS=fig.add_subplot(gs[0,:]); axT=fig.add_subplot(gs[1,:3]); axC=fig.add_subplot(gs[2,:3]); ax1=fig.add_subplot(gs[2,3]); axN=fig.add_subplot(gs[3,:])
# strip
axS.set_xlim(0,16); axS.set_ylim(0,1.7); axS.axis("off"); order=[s for s in STRAIN_ORDER if s!="C57BL_6"]; _WILD={"WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"}
for i,s in enumerate(order):
    has=s==STRAIN; axS.add_patch(Rectangle((i+0.12,0.5),0.76,0.62,fc=(pc.C_SILENCE if has else "white"),ec=("#333" if has else "#ccc"),lw=0.7))
    if has: axS.text(i+0.5,0.81,"●",ha="center",va="center",fontsize=8,color="white")
    axS.text(i+0.5,0.36,s.replace("_","/"),ha="right",va="top",fontsize=5.2,rotation=45,color=(pc.C_WILD if s in _WILD else "#333"),fontweight=("bold" if has else "normal"))
pc.pbadge(axS,"A",f"Insertion + piRNA SOURCE LOCUS present in {STRAIN.replace('_','/')} only (1/16) — strain-private (individual piRNAs, NOT a PICB cluster)",fs=8.0,y=1.32)
# TE
axT.set_xlim(S,E); axT.set_ylim(0,1.18); axT.axis("off")
axT.add_patch(Rectangle((S,0.58),N,0.34,fc="#f6f6f6",ec="#e3e3e3",lw=0.4))                                # TE lane (framed, like the locus figures)
for _,r in TE.iterrows(): axT.add_patch(Rectangle((max(r.s,S),0.61),min(r.e,E)-max(r.s,S),0.28,fc=TECOL.get(r.fam.split("|")[-1],"#ddd"),ec="white",lw=0.2))
axT.text(S-N*0.013,0.75,"TEs",ha="right",va="center",fontsize=6.4,fontweight="bold",color="#333",bbox=dict(boxstyle="round,pad=0.18",fc="#ececec",ec="#c9c9c9",lw=0.5))
axT.add_patch(Rectangle((S,0.14),N,0.34,fc="#edf2f9",ec="#d8e0ec",lw=0.4))                                # gene lane (framed)
_GEN=pc.genes_at(STRAIN,CHROM.split("chr")[-1],S,E); _gw=[min(E,ge)-max(S,gs) for (gs,ge,gx,gn,gb) in _GEN]
_sym=lambda nm: not str(nm).startswith(("ENSMUS","ENSMSP","ENSRNO"))                                       # real gene symbol vs bare Ensembl ID
_glbl=set(sorted([i for i in range(len(_GEN)) if _sym(_GEN[i][3]) or _gw[i]>=0.04*N],key=lambda i:(-int(_sym(_GEN[i][3])),-_gw[i]))[:5])   # NAME real symbols always; predicted IDs only if wide (not every gene)
for i,(gs,ge,gx,gn,gb) in enumerate(_GEN):
    a=max(S,gs); b=min(E,ge)
    if b<=a: continue
    axT.add_patch(Rectangle((a,0.17),b-a,0.28,fc="#c9d6ea",ec=pc.C_GENE,lw=0.7))
    if i in _glbl: axT.text((a+b)/2,0.085,gn,ha="center",va="top",fontsize=4.7,style="italic",fontweight="bold",color=pc.C_GENE)
axT.text(S-N*0.013,0.31,"genes",ha="right",va="center",fontsize=6.4,fontweight="bold",color="white",bbox=dict(boxstyle="round,pad=0.18",fc=pc.C_GENE,ec=pc.C_GENE,lw=0.5))
hx=0.13                                                                                                    # TE family key (top) — only families present at this locus
_pf=set(TE.fam.str.split("|").str[-1]); _key=[(f,TECOL[f]) for f in TECOL if f in _pf][:6]
for cl,col in _key: axT.add_patch(Rectangle((S+N*hx,1.02),N*0.011,0.11,fc=col,ec="none",clip_on=False)); axT.text(S+N*(hx+0.014),1.075,cl,fontsize=4.7,va="center"); hx+=0.155
# coverage
x=np.linspace(S,E,nb)
axC.fill_between(x,0,antiC,step="mid",color=pc.C_SILENCE,alpha=0.9,label="antisense to TE (silencing)")
axC.fill_between(x,0,-senseC,step="mid",color="#9a9a9a",alpha=0.85,label="sense to TE")
axC.axhline(0,color="#333",lw=0.8); axC.set_xlim(S,E); axC.set_ylabel("piRNA coverage",fontsize=8); axC.set_xlabel(f"{CHROM.split('#')[-1]} position (bp)",fontsize=8)
axC.ticklabel_format(axis="x",style="plain"); axC.tick_params(labelsize=6.5)
_h = nb // 2; _q = {"upper left": antiC[:_h].sum(), "upper right": antiC[_h:].sum(), "lower left": senseC[:_h].sum(), "lower right": senseC[_h:].sum()}
axC.legend(fontsize=6.2, frameon=False, loc=min(_q, key=_q.get))   # legend in the EMPTIEST corner — dynamic per the Panel B coverage (antisense=top, sense=bottom; left/right halves), so it never overlaps the bars
axC.spines[['top','right']].set_visible(False)
axC.set_ylim(-max(senseC.max(),1)*1.18, max(antiC.max(),1)*1.42)   # top headroom so the ping-pong-pair callout clears the peak
pc.pbadge(axC,"B",f"Strand-resolved piRNA coverage (24–32 nt) — antisense-to-TE = silencing (red) · sense = grey · TE {TELAB}",fs=7.5,y=1.17)   # raised clear of the ping-pong callout (at axes-frac y≈1.01)
# 1U
vals=[fc.get(n,0) for n in ["A","C","G","T"]]; ax1.bar(["A","C","G","U"],vals,color=["#cfcfcf","#cfcfcf","#cfcfcf",pc.C_1U]); ax1.set_ylim(0,max(vals+[1])*1.28)
for i,v in enumerate(vals): ax1.text(i,v+1,f"{v:.0f}",ha="center",fontsize=6,fontweight="bold")
ax1.set_title("5′ nt (1U)",fontsize=7.4,fontweight="bold",color=pc.C_1U); ax1.set_ylabel("%",fontsize=7); ax1.tick_params(labelsize=7); ax1.spines[['top','right']].set_visible(False)
# nucleotide panel
axN.axis("off"); bw=0.42
if pair:
    _,i,ss,sc,aq,acn=pair
    for k,b in enumerate(ss):
        axN.add_patch(Rectangle((k*bw,1.5),bw*0.9,0.5,fc=NT.get(b,"#999"),ec="white",lw=0.4)); axN.text(k*bw+bw*0.45,1.75,b,ha="center",va="center",fontsize=7,color="white",fontweight="bold")
    axN.text(-0.35,1.75,"5′",ha="right",fontsize=8); axN.text(len(ss)*bw+0.1,1.75,"3′ sense",ha="left",fontsize=7.5,color="#555")
    x0=9*bw
    for k,b in enumerate(aq):
        xx=x0-k*bw; axN.add_patch(Rectangle((xx+bw*0.1,0.7),bw*0.9,0.5,fc=NT.get(b,"#999"),ec="white",lw=0.4)); axN.text(xx+bw*0.55,0.95,b,ha="center",va="center",fontsize=7,color="white",fontweight="bold")
    axN.text(x0+bw+0.1,0.95,"5′",ha="left",fontsize=8); axN.text(x0-len(aq)*bw,0.95,"3′ antisense",ha="right",fontsize=7.5,color="#555")
    axN.add_patch(Rectangle((0,0.65),10*bw,1.45,fc="#fff3cd",ec="#e0a800",lw=1,zorder=0))
    for k in range(10): axN.plot([k*bw+bw*0.45]*2,[1.45,1.25],color="#888",lw=0.6)
    axN.annotate("1U",xy=(x0+bw*0.55,1.2),xytext=(x0+2,0.3),fontsize=8,color="#e31a1c",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#e31a1c",lw=1))
    axN.annotate("10A",xy=(9*bw+bw*0.45,2.0),xytext=(9*bw+1.6,2.5),fontsize=8,color="#1f78b4",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#1f78b4",lw=1))
    pc.pbadge(axN,"C","Nucleotide resolution — a real ping-pong pair: 10-nt complementary 5′ overlap, 1U (antisense) + 10A (sense)",fs=7.5,y=1.02)
    # genomic coordinate ruler aligned to the nucleotides (genomic = i + x/bw; both strands share the frame)
    gmin=i+10-len(aq); gmax=i+len(ss)-1; yr=0.45
    axN.plot([(gmin-i)*bw,(gmax-i)*bw],[yr,yr],color="#555",lw=0.8)
    g=int(np.ceil(gmin/5.0)*5)
    while g<=gmax:
        xx=(g-i)*bw; axN.plot([xx,xx],[yr,yr-0.05],color="#555",lw=0.8)
        axN.text(xx,yr-0.08,f"{g:,}",ha="right",va="top",fontsize=5,rotation=90,color="#555"); g+=5
    axN.text(((gmin+gmax)/2-i)*bw,yr-0.66,f"{CHROM.split('#')[-1]} position",ha="center",va="top",fontsize=7,color="#333")
    axN.set_xlim((10-len(aq))*bw-2.0,len(ss)*bw+3); axN.set_ylim(-0.62,2.7)
else:   # no sense+antisense pair (single-strand-dominant locus) → show the top piRNAs at base resolution instead of an empty panel
    _sdom = senseC.sum() >= antiC.sum()                          # sense-to-TE dominant? (senseC/antiC already account for TEST)
    _dom = (sense if TEST == "+" else anti) if _sdom else (anti if TEST == "+" else sense)   # sense/anti dicts are by genomic strand (fwd/rev); map to TE-relative via TEST
    _lab = "sense" if _sdom else "antisense"
    _top = sorted(((cnt, seq) for d in _dom.values() for seq, cnt in d.items()), reverse=True)[:7]
    _y = len(_top)
    for _cnt, _seq in _top:
        for k, b in enumerate(_seq):
            axN.add_patch(Rectangle((k*bw, _y-0.4), bw*0.9, 0.8, fc=NT.get(b, "#999"), ec="white", lw=0.4)); axN.text(k*bw+bw*0.45, _y, b, ha="center", va="center", fontsize=5.6, color="white", fontweight="bold")
        if _seq and _seq[0] == "T": axN.add_patch(Rectangle((0, _y-0.45), bw*0.9, 0.9, fc="none", ec="#e31a1c", lw=1.4))
        axN.text(len(_seq)*bw+0.2, _y, f"×{_cnt}", ha="left", va="center", fontsize=5.4, color="#666"); _y -= 1
    axN.text(-0.4, len(_top), "5′", ha="right", va="center", fontsize=7, color="#555")
    axN.set_xlim(-1.6, (max(len(s) for c, s in _top) if _top else 30)*bw+4.5); axN.set_ylim(-0.3, len(_top)+0.8)
    pc.pbadge(axN, "C", f"Nucleotide resolution — NO ping-pong pair (single-strand-dominant: {_lab}-to-TE); top piRNAs, 5′-U (1U) outlined red", fs=7.2, y=1.0)
if pair:   # zoom-in callout: tie the nucleotide panel to its true position on the coverage track (different scales)
    pc=(gmin+gmax)/2.0
    axC.axvspan(gmin,gmax,color="#e0a800",alpha=0.25,zorder=1)
    axC.text(pc,1.012,f"ping-pong pair  {CHROM.split('#')[-1]}:{gmin:,}  (zoom ↓)",ha="center",va="bottom",fontsize=5.2,color="#b8860b",fontweight="bold",transform=axC.get_xaxis_transform(),clip_on=False)   # above the axes (axes-frac y) so it never collides with bars on sense-dominant loci
    yC=axC.get_ylim()[0]; xL=(10-len(aq))*bw-0.3; xR=(len(ss)-1)*bw+0.3; yNt=2.12
    for xa,xb in [(gmin,xL),(gmax,xR)]:
        fig.add_artist(ConnectionPatch(xyA=(xa,yC),coordsA=axC.transData,xyB=(xb,yNt),coordsB=axN.transData,color="#e0a800",lw=0.9,ls=(0,(4,2)),alpha=0.9,zorder=20))
fig.suptitle(f"Strain-private {TELAB} → piRNA SOURCE LOCUS (individual strain-private piRNA sequences, not a PICB cluster): example + anatomy + nucleotide resolution ({STRAIN.replace('_','/')} {CHROM.split('#')[-1]}:{S:,}-{E:,})",fontsize=8.8,fontweight="bold",y=0.995)
for ext in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_locus_full_{tag}.{ext}",bbox_inches="tight")
print(f"[{tag}] 1U={fc.get('T',0):.0f}% anti={antiC.sum():.0f} sense={senseC.sum():.0f} TEs={len(TE)} pair={'yes' if pair else 'no'}")
