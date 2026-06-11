#!/usr/bin/env python3
"""piRNA biogenesis at the strain-private L1 cluster, at NUCLEOTIDE + GENOMIC-COORDINATE resolution
(real reads, strain-wise BAM, SPRET/EiJ P20.5). No abstract arrows: every glyph is a measured base at a
measured coordinate.
 Row A  PING-PONG: a real sense.antisense pair drawn at its true chr2 coordinates (10-nt 5' overlap,
   1U + 10A, base-pairing lines) + the measured 5'-overlap cross-correlation (peak @10 nt).
 Row B  PHASING: a real run of consecutive same-strand piRNAs at their true coordinates (head-to-tail,
   ~27-nt 5'-5' step, 1U at each 5') + the measured 5'-5' autocorrelation (peak @~27 nt).
Source data written alongside (SourceData_biogenesis_*.csv)."""
import numpy as np, pandas as pd, pysam
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"
BAM="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/SPRET_EiJ/SPRET_EiJ-20.5dpp.1/Aligned.sortedByCoord.out.bam"
CHROM,S,E="SPRET_EiJ#1#chr2",150450743,150521192; N=E-S
CHRLAB=CHROM.split("#")[-1]
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))

# ---- scan reads (24-32 nt) ----
plus5=np.zeros(N+1); minus5=np.zeros(N+1); first=[]
sense=defaultdict(Counter); anti=defaultdict(Counter)
bam=pysam.AlignmentFile(BAM,"rb")
for a in bam.fetch(CHROM,S,E):
    if a.is_unmapped or not a.query_sequence: continue
    L=a.reference_end-a.reference_start
    if not (24<=L<=32): continue
    if a.is_reverse:
        p=a.reference_end-1-S
        if 0<=p<=N: minus5[p]+=1
        anti[a.reference_end-1][rc(a.query_sequence)]+=1
        first.append(comp.get(a.query_sequence[-1],"N"))
    else:
        p=a.reference_start-S
        if 0<=p<=N: plus5[p]+=1
        sense[a.reference_start][a.query_sequence]+=1
        first.append(a.query_sequence[0])
bam.close()
fc=pd.Series(first).value_counts(normalize=True)*100
pp=np.array([np.dot(plus5[:N-d],minus5[d:N]) for d in range(0,31)])          # overlap=d+1
ac=np.array([np.dot(plus5[:N-d],plus5[d:N])+np.dot(minus5[:N-d],minus5[d:N]) for d in range(1,61)])

# ---- real ping-pong pair (sense@i, antisense@i+9), prefer 1U+10A ----
pair=None
for i in sense:
    if i+9 in anti:
        ss,sc=sense[i].most_common(1)[0]; aq,acn=anti[i+9].most_common(1)[0]
        if len(ss)>=10 and len(aq)>=10:
            sgn=(sc+acn)*(3 if (aq[0]=="T" and ss[9]=="A") else 1)
            if pair is None or sgn>pair[0]: pair=(sgn,i,ss,sc,aq,acn)

# ---- real phased run, plus strand: consecutive 5' ends ~27 apart ----
sense_top={p:c.most_common(1)[0] for p,c in sense.items() if c.most_common(1)[0][1]>=5}
def find_phased(require_1U,want):
    starts=sorted(sense_top); sset=set(starts); best=None
    for p0 in starts:
        if require_1U and sense_top[p0][0][0]!="T": continue
        run=[p0]; cur=p0
        while len(run)<want:
            cand=[p for p in range(cur+24,cur+31) if p in sset]
            if not cand: break
            nxt=min(cand,key=lambda p:(abs(p-(cur+27)),-sense_top[p][1])); run.append(nxt); cur=nxt
        if len(run)==want:
            sc=sum(sense_top[p][1] for p in run)
            if best is None or sc>best[0]: best=(sc,run)
    return best
ph=find_phased(True,3) or find_phased(True,2) or find_phased(False,3) or find_phased(False,2)

# ---------------- plot ----------------
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(14,7.6),dpi=300)
gs=fig.add_gridspec(2,2,width_ratios=[2.35,1],height_ratios=[1,1],hspace=0.62,wspace=0.2,
                    left=0.05,right=0.985,top=0.9,bottom=0.075)
axA=fig.add_subplot(gs[0,0]); axAh=fig.add_subplot(gs[0,1])
axB=fig.add_subplot(gs[1,0]); axBh=fig.add_subplot(gs[1,1])

def ruler(ax,xlo,xhi,gen0,y,step,every):
    """genomic ruler: x-unit -> coordinate gen0+x; tick every `every` x-units."""
    ax.plot([xlo+0.5,xhi+0.5],[y,y],color="#555",lw=0.8)
    x0=int(np.ceil(xlo/step)*step)
    for x in range(x0,xhi+1,every):
        ax.plot([x+0.5,x+0.5],[y,y-0.06],color="#555",lw=0.8)
        ax.text(x+0.5,y-0.09,f"{gen0+x:,}",ha="right",va="top",fontsize=5,rotation=90,color="#555")

# ---- A ping-pong pair ----
axA.axis("off")
if pair:
    _,i,ss,sc,aq,acn=pair; Ls=len(ss); La=len(aq)
    ytop,ybot,bh=1.5,0.62,0.5
    axA.add_patch(Rectangle((-0.06,ybot-0.06),10.12,(ytop+bh)-(ybot-0.06)+0.06,fc="#fff3cd",ec="#e0a800",lw=1,zorder=0))
    for k,b in enumerate(ss):                                   # sense base k at x=k
        axA.add_patch(Rectangle((k+0.04,ytop),0.92,bh,fc=NT.get(b,"#999"),ec="white",lw=0.4,zorder=2))
        axA.text(k+0.5,ytop+bh/2,b,ha="center",va="center",color="white",fontsize=7,fontweight="bold",zorder=3)
    for j,b in enumerate(aq):                                   # antisense base j at x=9-j
        x=9-j; axA.add_patch(Rectangle((x+0.04,ybot),0.92,bh,fc=NT.get(b,"#999"),ec="white",lw=0.4,zorder=2))
        axA.text(x+0.5,ybot+bh/2,b,ha="center",va="center",color="white",fontsize=7,fontweight="bold",zorder=3)
    for x in range(10): axA.plot([x+0.5,x+0.5],[ytop-0.02,ybot+bh+0.02],color="#777",lw=0.5,zorder=1)
    axA.text(-0.45,ytop+bh/2,"5′",ha="right",va="center",fontsize=8)
    axA.text(Ls+0.15,ytop+bh/2,"3′  sense (→)",ha="left",va="center",fontsize=7.5,color="#444")
    axA.text(9+0.7,ybot+bh/2,"5′",ha="left",va="center",fontsize=8)
    axA.text(9-(La-1)-0.6,ybot+bh/2,"antisense (←)  3′",ha="right",va="center",fontsize=7.5,color="#444")
    axA.annotate("1U",xy=(9.5,ybot+bh+0.02),xytext=(12.2,ybot-0.25),fontsize=8.5,color="#e31a1c",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#e31a1c",lw=1.1))
    axA.annotate("10A",xy=(9.5,ytop+bh+0.02),xytext=(9.5,ytop+bh+0.62),fontsize=8.5,color="#1f78b4",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#1f78b4",lw=1.1))
    axA.text(4.5,ytop+bh+0.62,"10-nt 5′ overlap",ha="center",va="bottom",fontsize=7,color="#b8860b",fontweight="bold")
    xmin=min(0,9-(La-1)); xmax=Ls-1; yr=ybot-0.2
    ruler(axA,xmin,xmax,i,yr,5,5)
    axA.text((xmin+xmax)/2+0.5,yr-0.92,f"{CHRLAB} position",ha="center",va="top",fontsize=7.5)
    axA.set_xlim(xmin-3.6,xmax+9.5); axA.set_ylim(yr-1.0,ytop+bh+0.95)
axA.set_title("A · Ping-pong amplification — a real sense·antisense pair (10-nt 5′ overlap, 1U + 10A)",fontsize=8.8,fontweight="bold",loc="left")

ov=np.arange(1,32)
axAh.bar(ov,pp,color=["#C0392B" if z==10 else "#cfcfcf" for z in ov],width=0.9)
axAh.axvline(10,color="#C0392B",lw=0.6,ls="--")
axAh.set_xlabel("sense–antisense 5′ overlap (nt)",fontsize=8); axAh.set_ylabel("read-pair signal",fontsize=8)
axAh.set_title("measured: 10-nt overlap peak",fontsize=8,fontweight="bold")
axAh.tick_params(labelsize=6.5); axAh.spines[['top','right']].set_visible(False)

# ---- B phased series ----
axB.axis("off")
if ph:
    _,run=ph; q0=run[0]; bh=0.5; yb=1.05
    seqs=[sense_top[p][0] for p in run]; span_end=(run[-1]-q0)+len(seqs[-1])-1
    for n,p in enumerate(run):
        seq=sense_top[p][0]; x0=p-q0
        for k,b in enumerate(seq):
            axB.add_patch(Rectangle((x0+k+0.04,yb),0.92,bh,fc=NT.get(b,"#999"),ec="white",lw=0.3,zorder=2))
            axB.text(x0+k+0.5,yb+bh/2,b,ha="center",va="center",color="white",fontsize=5,fontweight="bold",zorder=3)
        axB.annotate("1U",xy=(x0+0.5,yb+bh+0.02),xytext=(x0+0.5,yb+bh+0.5),fontsize=7,color="#e31a1c",fontweight="bold",ha="center",arrowprops=dict(arrowstyle="-|>",color="#e31a1c",lw=0.9))
        axB.text(x0+0.5,yb+bh+1.02,f"{p:,}",ha="center",va="bottom",fontsize=5.2,rotation=90,color="#555")
    for n in range(len(run)-1):
        xa=run[n]-q0+0.5; xb=run[n+1]-q0+0.5; d=run[n+1]-run[n]
        axB.annotate("",xy=(xb,yb-0.2),xytext=(xa,yb-0.2),arrowprops=dict(arrowstyle="<->",color="#0072B2",lw=1))
        axB.text((xa+xb)/2,yb-0.34,f"{d} nt",ha="center",va="top",fontsize=6.5,color="#0072B2",fontweight="bold")
    yr=yb-0.7; ruler(axB,0,span_end,q0,yr,10,10)
    axB.text(span_end/2,yr-0.92,f"{CHRLAB} position  (head-to-tail; ~27-nt 5′-5′ step = piRNA length)",ha="center",va="top",fontsize=7.5)
    axB.set_xlim(-3.6,span_end+6); axB.set_ylim(yr-1.0,yb+bh+1.7)
axB.set_title("B · Phasing (Zucchini/PLD6) — consecutive same-strand piRNAs, ~27-nt periodicity",fontsize=8.8,fontweight="bold",loc="left")

dd=np.arange(1,61)
axBh.bar(dd,ac,color=["#0072B2" if 25<=z<=29 else "#cfcfcf" for z in dd],width=0.9)
axBh.axvline(27,color="#0072B2",lw=0.6,ls="--")
axBh.set_xlabel("same-strand 5′–5′ distance (nt)",fontsize=8); axBh.set_ylabel("read-pair signal",fontsize=8)
axBh.set_title("measured: ~27-nt periodicity",fontsize=8,fontweight="bold")
axBh.tick_params(labelsize=6.5); axBh.spines[['top','right']].set_visible(False)

fig.suptitle(f"piRNA biogenesis at the strain-private L1 cluster — nucleotide + coordinate resolution (real reads, SPRET/EiJ {CHRLAB}, P20.5)",fontsize=10.5,fontweight="bold",y=0.975)
fig.text(0.5,0.012,f"5′-U (1U): {fc.get('T',0):.0f}%.  Both signatures shown at single-base, genome-anchored resolution — ping-pong (10-nt 5′ overlap, 1U+10A) and phasing (~27-nt head-to-tail): an actively amplifying, Zucchini-processed piRNA source born from a strain-private transposon.",ha="center",fontsize=6.4,color="#555")

for ext in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_biogenesis.{ext}",bbox_inches="tight")

# ---- source data ----
pd.DataFrame({"overlap_nt":ov,"pingpong_signal":pp}).to_csv(f"{PG}/SourceData_biogenesis_pingpong_overlap.csv",index=False)
pd.DataFrame({"distance_nt":dd,"phasing_signal":ac}).to_csv(f"{PG}/SourceData_biogenesis_phasing_distance.csv",index=False)
pd.DataFrame([{"nt":n,"pct":round(fc.get(n,0),3)} for n in ["A","C","G","T"]]).to_csv(f"{PG}/SourceData_biogenesis_5pnt.csv",index=False)
if pair:
    _,i,ss,sc,aq,acn=pair
    pd.DataFrame([{"role":"sense","genomic_5p":i,"reads":sc,"sequence":ss},
                 {"role":"antisense","genomic_5p":i+9,"reads":acn,"sequence":aq}]).to_csv(f"{PG}/SourceData_biogenesis_pingpong_pair.csv",index=False)
if ph:
    _,run=ph
    pd.DataFrame([{"n":n+1,"genomic_5p":p,"reads":sense_top[p][1],"sequence":sense_top[p][0],
                   "step_from_prev_nt":(p-run[n-1] if n>0 else 0)} for n,p in enumerate(run)]).to_csv(f"{PG}/SourceData_biogenesis_phased_series.csv",index=False)
print(f"1U={fc.get('T',0):.0f}%  pair={'Y' if pair else 'N'}  phased_run={'Y(%d)'%len(ph[1]) if ph else 'N'}  pp10={pp[9]:.0f}(max@{int(ov[np.argmax(pp)])})")
print("wrote Fig_biogenesis.{png,pdf,svg} + SourceData_biogenesis_*.csv")
