#!/usr/bin/env python3
"""Nucleotide-resolution SNP-variant (route 2 = NOT unique): a REAL SPRET strain-specific piRNA aligned
to the CONSERVED piRNA it differs from (the reference sequence at its homologous locus in CAST/C57,
which those strains also express). Actual sequences; SNP positions highlighted. Shows why a 1-3 mm
variant is a strain ALLELE, not a novel piRNA."""
import pysam, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))
d=pd.read_csv(f"{U}/step4/SPRET_EiJ.step4_classified.csv.gz")
sv=set(d.loc[d.klass=="SNP-variant of expressed (1-3mm)","id"])
bam=pysam.AlignmentFile(f"{U}/step4/SPRET_EiJ.cand_to_CAST_EiJ.Aligned.sortedByCoord.out.bam","rb")
pick=None
for a in bam.fetch(until_eof=True):
    if a.is_unmapped or a.query_name not in sv: continue
    nm=int(a.get_tag("NM"));
    if nm not in (2,3): continue                 # 2-3 SNPs = clear illustration
    try: ref=a.get_reference_sequence().upper()
    except Exception: continue
    q=a.query_sequence
    if len(q)!=len(ref) or not(25<=len(q)<=32): continue
    spret = rc(q) if a.is_reverse else q          # SPRET piRNA 5'->3'
    cons  = rc(ref) if a.is_reverse else ref       # conserved piRNA 5'->3'
    info=(a.reference_name,a.reference_start,a.reference_end,a.is_reverse)
    if spret[0]=="T":                              # prefer a 1U example
        pick=(spret,cons,nm,info); break
    if pick is None: pick=(spret,cons,nm,info)
bam.close()
spret,cons,nm,info=pick
refname,refstart,refend,isrev=info
snps=[i for i in range(len(spret)) if spret[i]!=cons[i]]
print(f"SPRET piRNA : {spret}\nconserved   : {cons}\nSNPs={len(snps)} at {snps}")
plt.rcParams.update({"font.family":"DejaVu Sans Mono"})
fig,ax=plt.subplots(figsize=(11,3.4),dpi=300); ax.axis("off")
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
bw=0.5
def row(seq,y,label,other=None):
    for k,b in enumerate(seq):
        snp = other is not None and b!=other[k]
        ax.add_patch(Rectangle((k*bw,y),bw*0.9,0.62,fc=NT.get(b,"#999"),ec=("black" if snp else "white"),lw=(1.8 if snp else 0.5),zorder=3))
        ax.text(k*bw+bw*0.45,y+0.31,b,ha="center",va="center",fontsize=8.5,color="white",fontweight="bold")
    ax.text(-0.3,y+0.31,"5′",ha="right",va="center",fontsize=9); ax.text(len(seq)*bw+0.1,y+0.31,f"3′  {label}",ha="left",va="center",fontsize=8.5,color="#333")
row(spret,1.5,"SPRET/EiJ piRNA (strain-specific by sequence)",other=cons)
row(cons,0.4,"conserved piRNA — also expressed in CAST/EiJ & C57BL/6NJ",other=spret)
for k in range(len(spret)):
    if spret[k]==cons[k]: ax.plot([k*bw+bw*0.45]*2,[1.45,1.05],color="#bbb",lw=0.7)
    else: ax.text(k*bw+bw*0.45,1.05,"✗",ha="center",va="center",fontsize=9,color="black",fontweight="bold")
# genomic coordinate ruler: anchor the conserved piRNA to its real homologous locus in CAST/EiJ
chrlab=refname.split("#")[-1]; yr=0.1
ax.plot([0,(len(cons)-1)*bw],[yr,yr],color="#555",lw=0.7)
for k in range(0,len(cons),5):
    g=(refstart+k) if not isrev else (refend-1-k)
    ax.plot([k*bw,k*bw],[yr,yr-0.05],color="#555",lw=0.7)
    ax.text(k*bw,yr-0.07,f"{g:,}",ha="right",va="top",fontsize=5,rotation=90,color="#555")
ax.text((len(cons)-1)*bw/2,yr-0.52,f"homologous conserved-piRNA locus — CAST/EiJ {chrlab}",ha="center",va="top",fontsize=6.3,color="#333")
ax.set_xlim(-2,len(spret)*bw+6); ax.set_ylim(-0.6,2.6)
ax.set_title(f"Route 2 — SNP-variant of a conserved piRNA ({len(snps)} SNPs) → a strain ALLELE, NOT a novel piRNA",fontsize=10,fontweight="bold")
ax.text(len(spret)*bw+0.4,1.27,f"{len(snps)} mismatches (✗)\n≤3 = data-driven\nPoisson cutoff\n→ NOT unique",fontsize=7.5,va="center",color="#b00")
fig.text(0.5,0.02,"The SPRET piRNA looked strain-specific only because its sequence carries a few SNPs; the SAME piRNA (conserved sequence) is expressed in the other strains. "
  "The STAR genome-anchored test catches this (≤3-mm match to an expressed sequence) and correctly calls it NOT unique.",ha="center",fontsize=6.4,color="#555")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_snp_variant_nucleotide.{e}",bbox_inches="tight")
print("wrote Fig_snp_variant_nucleotide.{png,pdf,svg}")
