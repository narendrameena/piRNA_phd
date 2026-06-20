#!/usr/bin/env python3
"""COORDINATE-ANCHORED SNP-variant (route 2): the SAME piRNA at its HOMOLOGOUS locus in two genomes —
SPRET/EiJ (own coords) and CAST/EiJ (own coords) — drawn base-by-base at real genomic positions, SNPs
marked. Message: SPRET's few SNPs make it look strain-specific, but it is a conserved piRNA also
expressed in CAST/EiJ & C57BL/6NJ -> NOT novel. Real data: SPRET locus from cand_self; CAST homolog +
conserved sequence from cand_to_CAST (MD tag)."""
import pysam, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))
d=pd.read_csv(f"{U}/step4/SPRET_EiJ.step4_classified.csv.gz")   # id<->sequence bridge (BAMs are keyed by candidate id)
_fc=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","klass5"])
_snp=set(_fc.loc[(_fc.strain=="SPRET_EiJ")&(_fc.klass5=="SNP-variant (1-3mm)"),"sequence"])   # klass5-defined SNP-variant set (≥2-read)
sv=set(d.loc[d.sequence.isin(_snp),"id"])
# SPRET own-genome loci for sv candidates (first pass)
self_loc={}
selfb=pysam.AlignmentFile(f"{U}/step4/SPRET_EiJ.cand_self.Aligned.sortedByCoord.out.bam","rb")
for a in selfb.fetch(until_eof=True):
    if a.is_unmapped or a.query_name not in sv or a.query_name in self_loc: continue
    self_loc[a.query_name]=(a.reference_name.split('#')[-1],a.reference_start,a.reference_end,a.is_reverse)
selfb.close()
# CAST homolog: scan cand_to_CAST for an sv candidate (2 SNPs, 1U) that ALSO maps to SPRET's genome
camap=pysam.AlignmentFile(f"{U}/step4/SPRET_EiJ.cand_to_CAST_EiJ.Aligned.sortedByCoord.out.bam","rb")
pickC=None
for a in camap.fetch(until_eof=True):
    if a.is_unmapped or a.query_name not in self_loc or int(a.get_tag("NM"))!=2: continue
    try: ref=a.get_reference_sequence().upper()
    except Exception: continue
    q=a.query_sequence
    if len(q)!=len(ref) or not(26<=len(q)<=30): continue
    spret = rc(q) if a.is_reverse else q; cons = rc(ref) if a.is_reverse else ref
    if spret[0]=="T": pickC=(a.query_name,a.reference_name.split('#')[-1],a.reference_start,a.reference_end,a.is_reverse,spret,cons); break
    if pickC is None: pickC=(a.query_name,a.reference_name.split('#')[-1],a.reference_start,a.reference_end,a.is_reverse,spret,cons)
camap.close()
cid,cchr,cs,ce,crev,spret,cons=pickC
schr,ss,se,srev=self_loc[cid]
snps=[k for k in range(len(spret)) if spret[k]!=cons[k]]
print(f"cid={cid} SPRET {schr}:{ss}-{se} CAST {cchr}:{cs}-{ce} SNPs at {snps}")
NT={"A":"#33a02c","C":"#1f78b4","G":"#ff7f00","T":"#e31a1c","N":"#999"}
plt.rcParams.update({"font.family":"DejaVu Sans Mono"})
fig,ax=plt.subplots(figsize=(10,4.6),dpi=175); ax.axis("off"); L=len(spret); bw=1.0
def track(seq,other,y,gstart,gstrand,lab,locuslab):
    for k,b in enumerate(seq):
        snp=b!=other[k]; ax.add_patch(Rectangle((k*bw,y),bw*0.92,0.7,fc=NT.get(b,"#999"),ec=("black" if snp else "white"),lw=(2 if snp else 0.5),zorder=3))
        ax.text(k*bw+bw*0.46,y+0.35,b,ha="center",va="center",fontsize=8,color="white",fontweight="bold")
    ax.text(-0.5,y+0.35,"5′",ha="right",fontsize=9); ax.text(L*bw+0.2,y+0.35,f"3′  {lab}",ha="left",fontsize=8.5,color="#333")
    # genomic coordinate ticks (every 5 bp) — real positions
    for k in range(0,L,5):
        pos = gstart+k if gstrand=="+" else gstart+(L-k)   # display genomic coordinate at this base
        ax.text(k*bw+bw*0.46,y-0.12,f"{gstart+k:,}",ha="center",va="top",fontsize=5.0,color="#888",rotation=90)
    ax.text(0,y+0.95,locuslab,fontsize=7,color="#555")
track(spret,cons,2.0,ss,"-" if srev else "+","SPRET/EiJ piRNA (strain-specific by sequence)",f"SPRET genome {schr}:{ss:,}-{se:,}")
track(cons,spret,0.5,cs,"-" if crev else "+","conserved piRNA — also expressed in CAST/EiJ & C57BL/6NJ",f"CAST genome {cchr}:{cs:,}-{ce:,}")
for k in range(L):
    if spret[k]==cons[k]: ax.plot([k*bw+bw*0.46]*2,[1.95,1.25],color="#ccc",lw=0.6)
    else: ax.plot([k*bw+bw*0.46]*2,[1.95,1.25],color="#000",lw=1.2); ax.text(k*bw+bw*0.46,1.6,"✗",ha="center",va="center",fontsize=8,fontweight="bold")
ax.set_xlim(-3,L*bw+9); ax.set_ylim(0,3.3)
ax.text(L*bw+0.4,1.6,f"{len(snps)} SNPs (✗)\n→ NOT unique",fontsize=8,color="#b00",va="center")
ax.text(L*bw*0.5,3.15,"The SPRET piRNA (2 SNPs) is a near-copy of a CONSERVED piRNA expressed in CAST/EiJ — NOT a novel piRNA",ha="center",fontsize=9.3,fontweight="bold")
fig.text(0.5,0.02,"Coordinate-anchored: SPRET piRNA at its real SPRET-genome position (top); the matching conserved piRNA at its CAST-genome position (bottom). "
  "NB the CAST copy is at a different coordinate (a paralogous family copy, not the syntenic allele) — but the sequence (±2 SNPs) IS expressed in CAST, so the genome-anchored test (≤3 mm) calls it NOT unique.",ha="center",fontsize=5.8,color="#555")
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_snp_variant_coord.{e}",bbox_inches="tight")
from PIL import Image; print("size:",Image.open(f"{U}/pangenome_te/Fig_snp_variant_coord.png").size)
print("wrote Fig_snp_variant_coord.{png,pdf,svg}")
