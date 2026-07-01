#!/usr/bin/env python3
"""
Characterise the cross-strain DIFFERENCE at Zamore pachytene piRNA loci by
TE content, sequence divergence, and genomic region.

Per locus:
  region_class : GRCm39 gene overlap (protein_coding > lncRNA > other_ncRNA > intergenic)
  te_fraction  : fraction of the locus span that is TE (project GRCm39->C57BL_6NJ via
                 cactus halLiftover [pangenome], intersect C57BL_6NJ RepeatMasker)
  dom_te_class : dominant TE class (LINE/SINE/LTR/DNA/...)
  divergence   : 1 - mean(fraction_lifted across strains)   [pangenome retention]
  sv_rate      : fraction of strains with a direct SV at the locus
Output: combined_rebuild/zamore_locus_annotation.csv
"""
import os, subprocess
import pandas as pd, numpy as np

BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CR   = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
SIF  = f"{BASE}/cactus_v2.9.3.sif"
HAL  = f"{BASE}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
GFF  = f"{BASE}/resources/black6/genome/Mus_musculus.GRCm39.115.chr.gff3"
RM   = f"{BASE}/resources/repeatMasker/C57BL_6NJ_repeatmasker.bed"
ZAM  = f"{CR}/_zamore_loci_noprefix.bed"          # GRCm39 chr/start/end/locus (bare-numeric)
W    = f"{CR}/_annot"; os.makedirs(W, exist_ok=True)

ref = pd.read_csv(ZAM, sep="\t", header=None, names=["chr","start","end","locus"])
ref["chr"]=ref["chr"].astype(str); ref["span"]=ref["end"]-ref["start"]
zam_sorted=f"{W}/zam.sorted.bed"
ref.sort_values(["chr","start"]).to_csv(zam_sorted,sep="\t",header=False,index=False)

# ── 1. genomic region (GRCm39 gff3) ───────────────────────────────────────────
genebed=f"{W}/genes.bed"
if not os.path.exists(genebed) or os.path.getsize(genebed)==0:
    subprocess.run(
        "grep -v '^#' %s | awk -F'\\t' '$3==\"gene\"||$3==\"ncRNA_gene\"{m=$9; "
        "i=index(m,\"biotype=\"); bt=substr(m,i+8); j=index(bt,\";\"); if(j>0)bt=substr(bt,1,j-1); "
        "print $1\"\\t\"$4-1\"\\t\"$5\"\\t\"bt}' | sort -k1,1 -k2,2n > %s" % (GFF, genebed),
        shell=True, check=True)
r=subprocess.run(["bedtools","intersect","-a",zam_sorted,"-b",genebed,"-wa","-wb"],
                 capture_output=True,text=True)
ov={}
for ln in r.stdout.strip().split("\n"):
    if not ln: continue
    f=ln.split("\t"); loc=f[3]; bt=f[7]
    ov.setdefault(loc,set()).add(bt)
def region_class(loc):
    bts=ov.get(loc,set())
    if "protein_coding" in bts: return "protein_coding"
    if "lncRNA" in bts: return "lncRNA"
    if bts: return "other_ncRNA"
    return "intergenic"
ref["region_class"]=ref["locus"].map(region_class)

# ── 2. TE content: halLiftover GRCm39 -> C57BL_6NJ, intersect C57BL_6NJ RM ─────
lifted=f"{W}/zam_in_C57.bed"
subprocess.run(["singularity","exec","--bind","/mnt",SIF,"halLiftover",HAL,"GRCm39",ZAM,"C57BL_6NJ",lifted],
               capture_output=True)
rmb=f"{W}/rm_bare.bed"
if not os.path.exists(rmb) or os.path.getsize(rmb)==0:
    subprocess.run(
        "awk -F'\\t' '$1~/^chr([0-9]+|X|Y)$/{sub(/^chr/,\"\",$1); n=split($4,a,\"|\"); "
        "cls=a[2]; print $1\"\\t\"$2\"\\t\"$3\"\\t\"cls}' %s | sort -k1,1 -k2,2n > %s" % (RM, rmb),
        shell=True, check=True)
lift_sorted=f"{W}/zam_in_C57.sorted.bed"
subprocess.run(f"sort -k1,1 -k2,2n {lifted} > {lift_sorted}", shell=True, check=True)
# projected span per locus (merge segments per locus)
lf=pd.read_csv(lifted,sep="\t",header=None,names=["chr","start","end","locus"])
proj_span=lf.assign(l=lf.end-lf.start).groupby("locus")["l"].sum().to_dict()
# TE bp per locus per class
r=subprocess.run(["bedtools","intersect","-a",lift_sorted,"-b",rmb,"-wo"],capture_output=True,text=True)
te_bp={}; te_cls={}
for ln in r.stdout.strip().split("\n"):
    if not ln: continue
    f=ln.split("\t"); loc=f[3]; cls=f[7]; ov_bp=int(f[-1])
    # collapse to broad TE class
    broad = cls.split("/")[0] if "/" in cls else cls
    if broad in ("LINE","SINE","LTR","DNA"):
        te_bp[loc]=te_bp.get(loc,0)+ov_bp
        te_cls.setdefault(loc,{}); te_cls[loc][broad]=te_cls[loc].get(broad,0)+ov_bp
def te_fraction(loc):
    ps=proj_span.get(loc,0)
    return (te_bp.get(loc,0)/ps) if ps>0 else np.nan
def dom_te(loc):
    d=te_cls.get(loc,{})
    return max(d,key=d.get) if d else "none"
ref["proj_span_C57"]=ref["locus"].map(lambda l: proj_span.get(l,0))
ref["te_fraction"]=ref["locus"].map(te_fraction)
ref["dom_te_class"]=ref["locus"].map(dom_te)

# ── 3. divergence (1 - mean retention) ────────────────────────────────────────
frac=pd.read_csv(f"{CR}/zamore_fraction_lifted.csv")
div=(1-frac.groupby("locus")["fraction_lifted"].mean())
ref["divergence"]=ref["locus"].map(div)

# ── 4. SV rate (fraction of strains with a direct SV) ─────────────────────────
sv=pd.read_csv(f"{CR}/all_strains_SV_matrix.csv")
svd=sv[sv.window=="direct"]
svr=svd.groupby("locus")["has_SV"].mean()
ref["sv_rate"]=ref["locus"].map(svr).fillna(0)

ref=ref.merge(frac.groupby("locus")["fraction_lifted"].mean().rename("mean_retention"),
              left_on="locus",right_index=True,how="left")
ref.to_csv(f"{CR}/zamore_locus_annotation.csv",index=False)
print(f"wrote zamore_locus_annotation.csv ({len(ref)} loci)")
print("\nregion_class:\n", ref.region_class.value_counts().to_string())
print("\ndom_te_class:\n", ref.dom_te_class.value_counts().to_string())
print(f"\nmean te_fraction: {ref.te_fraction.mean():.3f} | mean divergence: {ref.divergence.mean():.3f}")
print("\ndivergence by region:\n", ref.groupby("region_class")["divergence"].mean().round(3).to_string())
print("\ndivergence by dom_te_class:\n", ref.groupby("dom_te_class")["divergence"].mean().round(3).to_string())
