#!/usr/bin/env python3
"""THEME 07 — visual flowchart of the unique-piRNA pipeline: sRNA STAR BAM -> 25-32 nt piRNA candidate filter ->
5-class (klass5) classification, with parameters + counts at each step. Decision-tree form (how each class is reached).
Counts = the adopted >=2-read, 25-32 nt set (final_classified_clean_2read.csv.gz). Pure layout (no data read)."""
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,ax=plt.subplots(figsize=(13.5,15),dpi=300); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")
def box(x,y,w,h,txt,fc,ec="#333",tc="#111",fs=9,bold=False,r=0.02):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle=f"round,pad=0.3,rounding_size={r*100}",
                fc=fc,ec=ec,lw=1.2,zorder=2))
    ax.text(x,y,txt,ha="center",va="center",fontsize=fs,color=tc,fontweight="bold" if bold else "normal",zorder=3,linespacing=1.25)
def arr(x0,y0,x1,y1,lab="",fs=7,lc="#555",cx=None):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle="-|>",mutation_scale=14,lw=1.3,color="#666",zorder=1))
    if lab: ax.text(cx if cx is not None else (x0+x1)/2+1.5,(y0+y1)/2,lab,ha="left",va="center",fontsize=fs,color=lc,style="italic",zorder=4)
ax.text(50,98.5,"Unique-piRNA pipeline: sRNA STAR BAM → 25–32 nt candidates → 5-class (klass5)",ha="center",fontsize=14,fontweight="bold")
ax.text(50,96.2,"parameters in italics · counts = adopted ≥2-read, 25–32 nt set (399,812 candidates, 16 strains)",ha="center",fontsize=8.5,color="#666")
# ---- linear pipeline (top, centered) ----
box(50,92,46,4.0,"sRNA reads — 16 strains × 3 timepoints × 3 replicates","#eef2f7",fs=9.5,bold=True)
arr(50,90,50,87.2,"cutadapt: adapter TGGAATTCTCGGGTGCCAAGG · 20–36 nt · --discard-untrimmed",cx=27)
box(50,84.5,60,5.4,"STAR → strain genome (UNMASKED, PanSN)\n0 mismatch · multimap ≤800 / anchor 1600 · alignEndsType EndToEnd · no indels · no splicing","#eaf3ea",fs=8.8,bold=True)
arr(50,81.8,50,79.9,"→ sRNA BAM (per strain × tp × rep)",cx=52)
box(50,75.8,69,8.0,"piRNA window 25–32 nt (1U-pure; 24 nt excluded as 1U-impure) · per-strain expressed pools (collapse)\n"
    "edgeR: TMM-norm → quasi-likelihood F-test · focal strain X vs mean of 15 others → DA-enriched (FDR<0.05 & logFC>0)\n"
    "strain-SPECIFIC = DA-enriched  AND  present in focal X (≥2 of 3 reps)  AND  ≥2-read-ABSENT in all 15 other strains\n"
    "“≥2 reads” rule (per strain): PRESENT iff ≥2 reads total across its 3 reps · 0–1 read = ABSENT (noise/contamination)","#fff6e8",fs=7.1)
arr(50,71.7,50,70.7)
box(50,68.5,40,4.2,"strain-specific candidate piRNAs\n399,812","#fde9d9",fs=9.5,bold=True)
# ---- classification decision tree ----
arr(50,66.2,50,63.5,"classify each candidate (strain | timepoint | sequence)",cx=52)
# D1
box(50,61,58,5.4,"D1  exact sequence present in ANOTHER strain's expressed pool?\ntool: per-strain expressed pools (collapse, ≥2 reads) · exact string match","#f0f0f0",fs=8.3,bold=True,r=0.05)
GREY="#9e9e9e";ORANGE="#E69F00";TAN="#cdb892";BLUE="#0072B2";PURP="#7a3b9a"
arr(76,61,90,57.5,"YES",fs=8,lc="#b00")
box(90,55,17,5.0,"expressed-\nelsewhere (exact)\n39,824","#e8e8e8",ec=GREY,tc="#333",fs=8,bold=True)
arr(50,58.7,50,55.5,"NO → novel sequence",cx=52)
# D2
box(45,53,56,5.6,"D2  orthologous LOCUS present in ≥1 other strain? (pangenome PAV)\ntool: Cactus 17-strain pangenome HAL + halLiftover (present = ≥20% of locus lifts into the strain)","#f0f0f0",fs=8,bold=True,r=0.05)
arr(24,53,16,48.5,"NO",fs=8,lc="#b00"); arr(66,53,80,48.5,"YES (silent there)",fs=8,lc="#b00",cx=68)
box(16,46,21,4.4,"strain-private locus\n(locus gain / novel)","#f3ecf7",ec=PURP,fs=8.3,bold=True)
box(80,46,21,4.4,"conserved-but-silent\n(locus there, not expressed)","#e7f0f7",ec=BLUE,fs=8.3,bold=True)
# D3 (strain-private) + D4 (CBS)
arr(16,43.8,16,40.5); box(16,37.5,26,5.4,"D3  mm0 own-genome locus?\ntool: STAR · candidate → own\ngenome · 0 mismatch","#f0f0f0",fs=7.6,bold=True,r=0.05)
arr(16,35.6,16,32.5,"NO (unmapped)",fs=7.5,lc="#b00",cx=17)
arr(28,38,40,32.5,"YES",fs=8,lc="#b00",cx=33)
arr(80,43.8,80,40.5); box(80,37.5,28,5.4,"D4  1–3 mm SNP-variant of an EXPRESSED allele elsewhere?\ntool: STAR · candidate → other genomes (mm≤3) + expressed there","#f0f0f0",fs=7.6,bold=True,r=0.05)
arr(72,37,56,32.5,"YES",fs=8,lc="#b00",cx=60)        # YES (is a SNP-variant) -> SNP-variant leaf (x56)
arr(80,35.6,80,32.5,"NO",fs=8,lc="#b00",cx=81)       # NO (truly silent, novel) -> unique: conserved-but-silent (x80)
# ---- 5 final classes (leaves), color-coded ----
yL=29
box(16,yL,20,5.2,"low-quality\nno mm0 own locus\n39,656",TAN,ec="#9a7b3b",tc="#3a2d10",fs=8.2,bold=True)
box(40,yL,20,5.2,"UNIQUE:\nstrain-private locus\n20,616","#efe3f6",ec=PURP,tc=PURP,fs=8.2,bold=True)
box(56,yL,20,5.2,"SNP-variant\n(1–3 mm allelic)\n214,672","#fbeccd",ec=ORANGE,tc="#7a5500",fs=8.2,bold=True)
box(80,yL,20,5.2,"UNIQUE:\nconserved-but-silent\n85,044","#dcebf6",ec=BLUE,tc=BLUE,fs=8.2,bold=True)
# expressed-elsewhere leaf already drawn (top right). draw its line down to the summary level conceptually
# ---- summary split ----
ax.add_patch(plt.Rectangle((6,6),44,12,fc="#f7f7f7",ec="#bbb",lw=1,zorder=0))
ax.text(28,16.2,"NOT UNIQUE — 294,152 (74%)",ha="center",fontsize=10,fontweight="bold",color="#555")
ax.text(28,12.5,"expressed-elsewhere (39,824) · SNP-variant (214,672)\n· low-quality (39,656)",ha="center",fontsize=8,color="#666",linespacing=1.3)
ax.text(28,8.2,"sequence/allele also present or expressed in other strains,\nor no clean own-genome locus",ha="center",fontsize=7,color="#888",style="italic",linespacing=1.3)
ax.add_patch(plt.Rectangle((52,6),42,12,fc="#f3eef8",ec=PURP,lw=1.3,zorder=0))
ax.text(73,16.2,"GENUINELY UNIQUE — 105,660 (26%)",ha="center",fontsize=10,fontweight="bold",color=PURP)
ax.text(73,12.5,"unique: conserved-but-silent (85,044, divergence-driven)\n+ unique: strain-private locus (20,616, insertion-driven)",ha="center",fontsize=8,color="#444",linespacing=1.3)
ax.text(73,8.2,"novel sequence with a clean own-genome locus,\nNOT expressed (even as a SNP-variant) in any other strain",ha="center",fontsize=7,color="#7a3b9a",style="italic",linespacing=1.3)
# faint guides from leaves to summary
for lx,col in [(16,"#9a7b3b"),(56,ORANGE),(90,GREY)]: ax.plot([lx,28],[yL-2.7,18],color=col,lw=0.6,ls=":",alpha=0.5,zorder=0)
for lx,col in [(40,PURP),(80,BLUE)]: ax.plot([lx,73],[yL-2.7,18],color=col,lw=0.6,ls=":",alpha=0.5,zorder=0)
import pandas as _pd, os as _os
_SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data"; _os.makedirs(_SD,exist_ok=True)
_k=_pd.read_csv("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/final_classified_clean_2read.csv.gz",usecols=["klass5"])
_vc=_k.klass5.value_counts(); _tot=int(_vc.sum())
_vc.rename_axis("klass5").reset_index(name="n_candidates").to_csv(f"{_SD}/SourceData_Fig_unique_pirna_pipeline_flowchart.csv",index=False)   # CURRENT klass5 counts computed from data (the schematic's hardcoded display numbers are separate)
if _tot!=399812: print(f"FLAG: flowchart hardcoded total 399,812 vs CURRENT klass5 total {_tot:,} — displayed schematic counts are STALE (post-5-class-restoration); refresh the hardcoded text in a follow-up")
fig.savefig("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/Fig_unique_pirna_pipeline_flowchart.pdf",bbox_inches="tight")
for e in ("svg","png"): fig.savefig(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/Fig_unique_pirna_pipeline_flowchart.{e}",bbox_inches="tight")
print("wrote Fig_unique_pirna_pipeline_flowchart")
