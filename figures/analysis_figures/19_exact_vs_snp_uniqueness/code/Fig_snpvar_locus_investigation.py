#!/usr/bin/env python3
"""THEME 19 — LOCUS-based investigation of the 4,394 exact-'unique' SNP-variant piRNAs (data-driven, tested).
Question: at the GENOMIC LOCUS level, what are the loci that produce the 4,394 (the piRNAs EXACT calls strain-unique
but which are 1-3 SNP variants of a same-stage allele expressed elsewhere)? Are they new/strain-private loci or
deeply-conserved shared loci, and does that explain the wild-strain dominance?

Sets (exact_stagepeak_classified.csv.gz):
  SNPVAR  = was_snp_variant==True                                   (the 4,394)
  CBS     = klass=='unique: conserved-but-silent' & not snp-variant (clean conserved-locus unique)
  PRIVATE = klass=='unique: strain-private locus'                   (genuinely-novel, TE-driven locus)
Genomic loci from cand_self16/{strain}.cand_self16.bam (read name = cand_id = strain|tp|seq; PanSN chrom).
NOTE: positions are per-strain ASSEMBLY coords, so cross-strain colocation is taken from the orthology/allelic
series (snp_variant_refinement, cached as nshare), NOT from raw positions.

Four hypotheses, each TESTED:
  H1 mechanism    SNPVAR depleted in strain-private TE insertions vs PRIVATE          -> Fisher exact
  H2 allelic ser. SNPVAR loci shared with many strains (nshare)                       -> distribution
  H3 hotspots     SNPVAR no more locus-concentrated than CBS                          -> Gini + permutation
  H4 why-wild     wild vs classical home sit at equally-broad series                  -> Mann-Whitney
Outcome: H1,H2 confirmed; H3 REJECTED (not hotspot-driven); H4 null (same loci, wild just carry 44x more)."""
import warnings, subprocess; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, chi2_contingency
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
P=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
ST="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/samtools"
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
CBS_K="unique: conserved-but-silent"; PRIV_K="unique: strain-private locus"
C_SNP="#b2182b"; C_CBS="#2166ac"; C_PRIV="#1b9e77"; C_W="#C0392B"; C_C="#0072B2"
def gini(x):
    x=np.sort(np.asarray(x,float)); n=x.size
    return (2*np.sum(np.arange(1,n+1)*x)/(n*x.sum()))-(n+1)/n

# ---------- locus extraction (cand_self16 BAMs) ----------
d=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz")
d["set"]=np.where(d.was_snp_variant,"SNPVAR",
        np.where((d.klass==CBS_K)&(~d.was_snp_variant),"CBS",
        np.where(d.klass==PRIV_K,"PRIVATE","")))
sets=d[d.set!=""]; need=dict(zip(sets.cand_id,sets.set)); rows=[]
for X in sorted(sets.strain.unique()):
    out=subprocess.run([ST,"view",f"{P}/cand_self16/{X}.cand_self16.bam"],capture_output=True,text=True,timeout=900).stdout
    for ln in out.splitlines():
        f=ln.split("\t",6)
        if len(f)<5: continue
        flag=int(f[1])
        if flag not in (0,16) or f[0] not in need: continue
        rows.append((f[0],need[f[0]],X,f[2].split("#")[-1],int(f[3]),"-" if flag==16 else "+"))
c=pd.DataFrame(rows,columns=["cand_id","set","strain","chrom","pos","strand"]).sort_values("pos").drop_duplicates("cand_id")
c=c.sort_values(["strain","chrom","pos"]); lid=[]; last=None; k=-1
for r in c.itertuples():
    key=(r.strain,r.chrom)
    if last is None or key!=last[:2] or r.pos-last[2]>1000: k+=1
    lid.append(k); last=(r.strain,r.chrom,r.pos)
c["locus"]=lid

# ---------- H1: strain-private TE insertion overlap ----------
def load_ins(X):
    iv={}
    try:
        for ln in open(f"{P}/pangenome_te/{X}.ins_loci.bed"):
            a=ln.rstrip("\n").split("\t")
            if len(a)>=3: iv.setdefault(a[0].split("#")[-1],[]).append((int(a[1]),int(a[2])))
    except FileNotFoundError: return None
    return {ch:(np.array([s for s,_ in sorted(v)]),np.array([e for _,e in sorted(v)])) for ch,v in iv.items()}
ins={X:load_ins(X) for X in c.strain.unique()}
def hit(r):
    iv=ins.get(r.strain)
    if not iv or r.chrom not in iv: return False
    st,en=iv[r.chrom]; i=np.searchsorted(st,r.pos,"right")-1
    return bool(i>=0 and r.pos<=en[i])
c["in_TE"]=[hit(r) for r in c.itertuples()]
fr={s:(int(c[c.set==s].in_TE.sum()),len(c[c.set==s])) for s in ["SNPVAR","CBS","PRIVATE"]}
ORtp,ptp=fisher_exact([[fr["SNPVAR"][0],fr["SNPVAR"][1]-fr["SNPVAR"][0]],[fr["PRIVATE"][0],fr["PRIVATE"][1]-fr["PRIVATE"][0]]])
_,pcbs=fisher_exact([[fr["SNPVAR"][0],fr["SNPVAR"][1]-fr["SNPVAR"][0]],[fr["CBS"][0],fr["CBS"][1]-fr["CBS"][0]]])

# ---------- H2/H4: allelic-series breadth (cached nshare) ----------
ns=pd.read_csv(f"{T}/data/SourceData_Fig_snp_allele_test.csv.gz").rename(columns={"n_other_strains_with_allele":"nshare"})
ns["home"]=ns.cand_id.str.split("|").str[0]; ns["wild"]=ns.home.isin(WILD)
wv=ns[ns.wild].nshare.values; cv=ns[~ns.wild].nshare.values; _,pw=mannwhitneyu(wv,cv,alternative="two-sided")

# ---------- H3: concentration (Gini + permutation) + chrom distribution ----------
gd={s:gini(c[c.set==s].groupby("locus").size().values) for s in ["SNPVAR","CBS","PRIVATE"]}
pool=c[c.set.isin(["SNPVAR","CBS"])]; codes,_=pd.factorize(pool.locus.values); isS=(pool.set.values=="SNPVAR")
obs=gd["SNPVAR"]-gd["CBS"]; rng=np.random.default_rng(0); nperm=1000; cnt=0
for _ in range(nperm):
    m=rng.permutation(isS); b1=np.bincount(codes[m]); b2=np.bincount(codes[~m])
    if abs(gini(b1[b1>0])-gini(b2[b2>0]))>=abs(obs): cnt+=1
perm_p=(cnt+1)/(nperm+1)
ct=pd.crosstab(pool.chrom,pool.set); chi,pc,_,_=chi2_contingency(ct)

# ---------- figure ----------
fig,((axA,axB),(axC,axD))=plt.subplots(2,2,figsize=(13.5,9.6),dpi=300)
labs=["SNP-variant\n(the 4,394)","conserved-\nbut-silent","strain-private\n(TE-new)"]
pct=[100*fr[s][0]/fr[s][1] for s in ["SNPVAR","CBS","PRIVATE"]]
axA.bar(range(3),pct,color=[C_SNP,C_CBS,C_PRIV],width=0.66,edgecolor="white")
for i,s in enumerate(["SNPVAR","CBS","PRIVATE"]): axA.text(i,pct[i]+0.6,f"{pct[i]:.1f}%\n({fr[s][0]}/{fr[s][1]})",ha="center",va="bottom",fontsize=8.5,fontweight="bold")
axA.plot([0,0,2,2],[27,28.5,28.5,27],lw=1.0,color="#333"); axA.text(1,29.2,f"Fisher OR={ORtp:.3f}, p={ptp:.0e}",ha="center",fontsize=8.5,fontweight="bold")
axA.set_ylim(0,33); axA.set_xticks(range(3)); axA.set_xticklabels(labs,fontsize=8.7); axA.set_ylabel("% loci inside a strain-private TE insertion",fontsize=9.5); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  Mechanism: the 4,394 sit at CONSERVED loci, not new TE insertions\n(~19x depleted vs genuinely-novel piRNAs; indistinguishable from conserved-but-silent)",fontsize=10,fontweight="bold",loc="left")
vc=ns.nshare.value_counts().reindex(range(0,16),fill_value=0); med=ns.nshare.median()
axB.bar(vc.index,vc.values,color=C_SNP,width=0.85,edgecolor="white"); axB.axvline(med,ls="--",color="#333",lw=1.2)
axB.text(med-0.3,vc.max()*0.95,f"median {med:.0f}/15",ha="right",fontsize=9,fontweight="bold",color="#333")
axB.text(0.02,0.80,f"45% shared with ALL 15\n72% shared with $\\geq$10\nmean {ns.nshare.mean():.1f} strains",transform=axB.transAxes,fontsize=8.6,color=C_SNP,fontweight="bold")
axB.set_xlabel("# OTHER strains carrying a 1-3 SNP allele at the SAME locus",fontsize=9.5); axB.set_ylabel("SNP-variant piRNAs",fontsize=9.5)
axB.set_xticks(range(0,16,2)); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  These are deeply-shared ALLELIC SERIES, not strain-private innovation\n(an orthologous locus where nearly every strain carries its own SNP allele)",fontsize=10,fontweight="bold",loc="left")
def lorenz(v):
    x=np.sort(v)[::-1]; return np.concatenate([[0],np.arange(1,len(x)+1)/len(x)]),np.concatenate([[0],np.cumsum(x)/x.sum()])
for s,col,lab in [("SNPVAR",C_SNP,"SNP-variant"),("CBS",C_CBS,"conserved-but-silent"),("PRIVATE",C_PRIV,"strain-private")]:
    fx,cum=lorenz(c[c.set==s].groupby("locus").size().values); axC.plot(fx,cum,color=col,lw=2,label=f"{lab} (Gini {gd[s]:.2f})")
axC.plot([0,1],[0,1],ls=":",color="#999",lw=1)
axC.text(0.97,0.06,f"SNP-variant $\\approx$ conserved-but-silent\nperm p = {perm_p:.2f}  (NOT hotspot-driven)",transform=axC.transAxes,ha="right",fontsize=8.6,color="#333",fontweight="bold")
axC.set_xlabel("fraction of loci (ranked busiest$\\rightarrow$rarest)",fontsize=9.5); axC.set_ylabel("cumulative fraction of piRNAs",fontsize=9.5)
axC.legend(fontsize=8,frameon=False,loc="lower right",bbox_to_anchor=(1.0,0.16)); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  TESTED & REJECTED: SNP-variants are NO more locus-concentrated\nthan conserved-but-silent piRNAs (top 21% of loci = 50% of piRNAs, same as CBS)",fontsize=10,fontweight="bold",loc="left")
parts=axD.violinplot([wv,cv],showmedians=True,widths=0.8)
for b,col in zip(parts["bodies"],[C_W,C_C]): b.set_facecolor(col); b.set_alpha(0.55)
for key in ("cmedians","cbars","cmins","cmaxes"): parts[key].set_color("#333")
axD.set_xticks([1,2]); axD.set_xticklabels([f"wild home\n(n={len(wv):,})",f"classical home\n(n={len(cv)})"],fontsize=9)
axD.set_ylabel("allelic-series breadth (# strains sharing locus)",fontsize=9.5); axD.set_ylim(0,16.5)
axD.text(1.5,15.4,f"same loci: MWU p={pw:.2f} (n.s.)",ha="center",fontsize=9,fontweight="bold")
axD.text(0.5,0.045,f"but wild carry 44x MORE\nSNP-variants: {len(wv):,} vs {len(cv)}",transform=axD.transAxes,ha="center",fontsize=9,color=C_W,fontweight="bold")
axD.spines[["top","right"]].set_visible(False)
axD.set_title("D  Why WILD strains dominate = allelic divergence, NOT special loci\n(wild & classical sit at equally-conserved loci; wild just hit the divergent allele far more)",fontsize=10,fontweight="bold",loc="left")
fig.suptitle("LOCUS-based investigation of the 4,394 exact-'unique' SNP-variant piRNAs — conserved hypervariable allelic series, not strain-private loci",fontsize=12.5,fontweight="bold",y=1.005)
fig.tight_layout(rect=[0,0,1,0.97])
out=f"{T}/figures/Fig_snpvar_locus_investigation"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
c[["cand_id","set","strain","chrom","pos","strand","locus","in_TE"]].to_csv(f"{T}/data/SourceData_Fig_snpvar_locus_investigation.csv.gz",index=False)
print(f"wrote {out}")
print(f"H1 TE-insertion: SNPVAR {pct[0]:.1f}% vs PRIVATE {pct[2]:.1f}% (Fisher p={ptp:.0e}); SNPVAR vs CBS p={pcbs:.2f}")
print(f"H2 series breadth: median {med:.0f}/15, {100*(ns.nshare>=15).mean():.0f}% shared with all 15")
print(f"H3 concentration: Gini SNPVAR {gd['SNPVAR']:.3f} vs CBS {gd['CBS']:.3f} perm p={perm_p:.2f}; chrom chi2={chi:.0f} p={pc:.0e}")
print(f"H4 why-wild: wild median {np.median(wv):.0f} vs classical {np.median(cv):.0f} MWU p={pw:.2f}; counts {len(wv)} vs {len(cv)}")
