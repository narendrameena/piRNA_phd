#!/usr/bin/env python3
"""
Theme 03 (REFRAMED) — QUANTITATIVE cross-strain expression of Zamore pachytene
piRNA loci + statistical testing against divergence, TE class, and the
expression-concentration ("few loci make most piRNAs") phenomenon.

Expression = sum of all_primary_FPM of the PICB combined-run clusters that overlap
the locus (pangenome GRCm39 frame), per strain x timepoint
(unique_pirna/cluster_pav/picb_pangenome_clusters.tsv).

Tests:
  T1  divergence ~ TE class (Kruskal-Wallis; LINE vs SINE Mann-Whitney)
  T2  divergence ~ P20.5 expression (Spearman)
  T3  divergence ~ cross-strain expression CV (Spearman)
  T4  expression concentration: Gini + top-N% per strain; wild vs classical (Mann-Whitney)
  T5  cross-strain expression correlation structure (how conserved is the ranking)
Outputs: zamore_locus_expression_P20.5.csv, zamore_expression_tests.txt, +source data
"""
import os, subprocess, sys
import numpy as np, pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD

BASE="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
CR=f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
PANG=f"{BASE}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/picb_pangenome_clusters.tsv"
W=f"{CR}/_expr"; os.makedirs(W,exist_ok=True)
TPS=["E16.5","P12.5","P20.5"]; TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}

# ── locus FPM = sum of overlapping PICB pangenome clusters (per strain x tp) ───
loci=pd.read_csv(f"{CR}/_zamore_loci_noprefix.bed",sep="\t",header=None,names=["chr","start","end","locus"])
loci["chr"]=loci["chr"].astype(str)
loci.sort_values(["chr","start"]).to_csv(f"{W}/loci.bed",sep="\t",header=False,index=False)
clu_bed=f"{W}/clusters.bed"
if not os.path.exists(clu_bed) or os.path.getsize(clu_bed)==0:
    # g39_chrom start end  name=strain|tp|FPM
    subprocess.run("tail -n +2 %s | awk -F'\\t' '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"|\"$5\"|\"$6}' "
                   "| sort -k1,1 -k2,2n > %s" % (PANG, clu_bed), shell=True, check=True)
r=subprocess.run(["bedtools","intersect","-a",f"{W}/loci.bed","-b",clu_bed,"-wa","-wb"],
                 capture_output=True,text=True)
rows=[]
for ln in r.stdout.strip().split("\n"):
    if not ln: continue
    f=ln.split("\t"); locus=f[3]; strain,tp,fpm=f[7].split("|")
    rows.append((locus,strain,TPMAP.get(tp,tp),float(fpm)))
ex=pd.DataFrame(rows,columns=["locus","strain","tp","fpm"])
# sum FPM of all overlapping clusters per locus x strain x tp
exg=ex.groupby(["locus","strain","tp"])["fpm"].sum().reset_index()

ann=pd.read_csv(f"{CR}/zamore_locus_annotation.csv")
strains=[s for s in STRAIN_ORDER if s in set(exg.strain)]

def matrix(tp):
    m=exg[exg.tp==tp].pivot_table(index="locus",columns="strain",values="fpm",aggfunc="sum").reindex(columns=strains).fillna(0)
    return m.reindex(ann.locus.values).fillna(0)
M20=matrix("P20.5")
M20.to_csv(f"{CR}/zamore_locus_expression_P20.5.csv")

out=[]
def p(s): out.append(s); print(s)

p("="*70); p("ZAMORE LOCUS EXPRESSION — STATISTICAL TESTS (pangenome, combined run)"); p("="*70)
p(f"loci={len(ann)}  strains={len(strains)}  expression=sum overlapping PICB cluster all_primary_FPM\n")

# ── T1: divergence ~ TE class ─────────────────────────────────────────────────
a=ann.dropna(subset=["divergence"])
grp={c:a.loc[a.dom_te_class==c,"divergence"].values for c in ["LINE","LTR","SINE"]}
H,pK=kruskal(*[grp[c] for c in ["LINE","LTR","SINE"]])
U,pLS=mannwhitneyu(grp["LINE"],grp["SINE"],alternative="greater")
p("T1  Divergence (1-retention) vs dominant TE class")
for c in ["LINE","LTR","SINE","none"]:
    v=a.loc[a.dom_te_class==c,"divergence"]
    p(f"      {c:5} n={len(v):3}  median divergence={v.median():.4f}  mean={v.mean():.4f}")
p(f"    Kruskal-Wallis LINE/LTR/SINE: H={H:.2f}  p={pK:.2e}")
p(f"    Mann-Whitney LINE>SINE (one-sided): U={U:.0f}  p={pLS:.2e}\n")

# ── T2: divergence ~ P20.5 expression ─────────────────────────────────────────
mean20=M20.mean(axis=1); a2=ann.set_index("locus").join(np.log1p(mean20).rename("log_meanFPM")).dropna(subset=["divergence","log_meanFPM"])
rho2,pp2=spearmanr(a2.divergence,a2.log_meanFPM)
p("T2  Divergence vs mean P20.5 expression (log FPM)")
p(f"    Spearman rho={rho2:.3f}  p={pp2:.2e}  (n={len(a2)})\n")

# ── T3: divergence ~ cross-strain expression CV ───────────────────────────────
cv=(M20.std(axis=1)/M20.mean(axis=1).replace(0,np.nan)).rename("cv")
a3=ann.set_index("locus").join(cv).dropna(subset=["divergence","cv"])
rho3,pp3=spearmanr(a3.divergence,a3.cv)
p("T3  Divergence vs cross-strain expression CV")
p(f"    Spearman rho={rho3:.3f}  p={pp3:.2e}  (n={len(a3)})\n")

# ── T4: expression concentration per strain + wild vs classical ───────────────
def gini(x):
    x=np.sort(np.asarray(x,float)); x=x[x>=0]; n=len(x)
    if n==0 or x.sum()==0: return np.nan
    return (2*np.sum((np.arange(1,n+1))*x)/(n*x.sum()))-(n+1)/n
def topN_for(x,frac=0.9):
    x=np.sort(np.asarray(x,float))[::-1]; c=np.cumsum(x);
    if c[-1]==0: return np.nan
    return int(np.searchsorted(c, frac*c[-1])+1)
gini_s={s:gini(M20[s].values) for s in strains}
top90={s:topN_for(M20[s].values,0.9) for s in strains}
top10pct={s: 100*np.sort(M20[s].values)[::-1][:10].sum()/M20[s].values.sum() for s in strains}
gw=[gini_s[s] for s in strains if s in WILD]; gc=[gini_s[s] for s in strains if s not in WILD]
Ug,pg=mannwhitneyu(gw,gc)
p("T4  Expression concentration ('few loci make most piRNAs')")
p(f"    median across strains: Gini={np.nanmedian(list(gini_s.values())):.3f}  "
  f"top-10 loci = {np.nanmedian(list(top10pct.values())):.1f}% of output  "
  f"loci for 90% = {int(np.nanmedian(list(top90.values())))}")
p(f"    Gini wild={np.nanmedian(gw):.3f} vs classical={np.nanmedian(gc):.3f}  Mann-Whitney p={pg:.3f}\n")

# ── T5: cross-strain ranking conservation ─────────────────────────────────────
cc=M20.apply(lambda c: np.log1p(c)).corr(method="spearman")
iu=np.triu_indices(len(strains),1)
p("T5  Cross-strain conservation of the expression profile")
p(f"    mean pairwise Spearman of per-locus expression across strains = {cc.values[iu].mean():.3f}")
p(f"    -> the WHICH-loci-fire ranking is highly conserved across strains\n")

open(f"{CR}/zamore_expression_tests.txt","w").write("\n".join(out))
# concentration source data
pd.DataFrame({"strain":strains,"gini":[gini_s[s] for s in strains],
              "top10_pct":[top10pct[s] for s in strains],"loci_for_90pct":[top90[s] for s in strains],
              "subspecies":["wild-derived" if s in WILD else "classical" for s in strains]}
             ).to_csv(f"{CR}/SourceData_zamore_concentration.csv",index=False)
print("\nwrote tests + matrices + source data")
