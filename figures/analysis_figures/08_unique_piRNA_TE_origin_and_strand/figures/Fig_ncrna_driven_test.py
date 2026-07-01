#!/usr/bin/env python3
"""ncRNA-driven test (the lncRNA analog of the pangenome TE-driven test). For each Step-4 piRNA class, the
fold-enrichment of UNIQUELY-mapping piRNA production loci inside annotated lncRNA genes, over the random-locus
expectation (lncRNA bp / genome bp), per pilot strain. Data-driven: cand_self own-genome loci (Step 4) x
v3.5 GFF3 lncRNA genes (biotype=lncRNA). Uniquely-mapping only (1 own-genome locus) => clean fold-enrichment
(no multimapping inflation). Main chromosomes only. Mirrors Fig_te_driven_coord."""
import warnings; warnings.filterwarnings("ignore")
import pysam, pandas as pd, numpy as np, bisect
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"; ANN="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation"
STR=["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]; COL={"C57BL_6NJ":"#1f77b4","CAST_EiJ":"#2ca02c","SPRET_EiJ":"#d95f02"}
CMAP={"expressed elsewhere (exact)":"expressed elsewhere\n(common)","SNP-variant of expressed (1-3mm)":"SNP-variant",
      "unique: conserved-but-silent":"unique:\nconserved-silent","unique: strain-private locus":"unique:\nstrain-private locus"}
CLASSES=list(dict.fromkeys(CMAP.values())); MAIN={f"chr{i}" for i in list(range(1,20))+["X","Y"]}
rows=[]
for X in STR:
    d=pd.read_csv(f"{U}/step4/{X}.step4_classified.csv.gz"); cls=dict(zip(d.id,d.klass))
    bam=pysam.AlignmentFile(f"{U}/step4/{X}.cand_self.Aligned.sortedByCoord.out.bam","rb")
    gsize=sum(l for r,l in zip(bam.references,bam.lengths) if r.split("#")[-1] in MAIN)
    nloc=Counter(); loc={}
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped: continue
        c=a.reference_name.split("#")[-1]
        if c in MAIN: nloc[a.query_name]+=1; loc[a.query_name]=(c,a.reference_start,a.reference_end)
    bam.close()
    # CONFOUNDING-FREE: lncRNA space with protein_coding SUBTRACTED (drops antisense/intronic lncRNAs that
    # overlap a coding gene, e.g. Gm52992/Dlg2). A locus counts as clean-lncRNA only if it is in a lncRNA AND
    # in no protein_coding gene; the random expectation uses clean-lncRNA bp (lncRNA minus protein_coding).
    lnc=defaultdict(list); pcg=defaultdict(list)
    for ln in open(f"{ANN}/{X}_v3.5.gff3"):
        if ln[0]=="#": continue
        f=ln.split("\t")
        if len(f)<9 or f[2] not in ("gene","ncRNA_gene"): continue
        c=f[0].split("#")[-1]
        if c not in MAIN: continue
        if "biotype=lncRNA" in f[8]: lnc[c].append((int(f[3]),int(f[4])))
        elif "biotype=protein_coding" in f[8]: pcg[c].append((int(f[3]),int(f[4])))
    def merge(dd):
        out={}
        for c,iv in dd.items():
            iv.sort(); m=[]
            for s,e in iv:
                if m and s<=m[-1][1]: m[-1]=[m[-1][0],max(m[-1][1],e)]
                else: m.append([s,e])
            out[c]=m
        return out
    lncm=merge(lnc); pcm=merge(pcg)
    pss={c:[s for s,e in pcm[c]] for c in pcm}; pee={c:[e for s,e in pcm[c]] for c in pcm}
    def in_pc(c,s,e):
        if c not in pss: return False
        i=bisect.bisect_right(pss[c],e)-1
        return i>=0 and pee[c][i]>s
    cleanbp=0; ss={}; ee={}
    for c,m in lncm.items():
        ss[c]=[s for s,e in m]; ee[c]=[e for s,e in m]
        for s,e in m:
            ov=0
            if c in pss:
                for k in range(max(0,bisect.bisect_left(pss[c],s)-1), bisect.bisect_left(pss[c],e)+1):
                    if k<len(pcm[c]): a,b=pcm[c][k]; ov+=max(0,min(b,e)-max(a,s))
            cleanbp+=(e-s)-ov
    lncfrac=cleanbp/gsize
    def inlnc(c,s,e):
        if c not in ss: return False
        i=bisect.bisect_right(ss[c],e)-1
        return i>=0 and ee[c][i]>s and not in_pc(c,s,e)
    print(f"{X}: genome(main)={gsize/1e9:.2f} Gb | clean-lncRNA(pc-excluded)={cleanbp/1e6:.0f} Mb ({lncfrac*100:.1f}%)")
    for orig,lab in CMAP.items():
        ids=[i for i in cls if cls[i]==orig and nloc.get(i)==1]
        n_in=sum(inlnc(*loc[i]) for i in ids); fold=(n_in/len(ids))/lncfrac if ids and lncfrac else np.nan
        rows.append(dict(strain=X,cls=lab,fold=fold,n_uniq=len(ids),n_in_clean_lncRNA=n_in,obs_pct=round(100*n_in/len(ids),2) if ids else np.nan))
        print(f"   {lab.replace(chr(10),' ')}: uniq={len(ids):,} in-lncRNA={n_in:,} ({100*n_in/max(len(ids),1):.1f}%) fold={fold:.2f}x")
df=pd.DataFrame(rows); df.to_csv(f"{PG}/SourceData_ncrna_driven_test.csv",index=False)
# ---- plot (grouped bars, log y) ----
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(11,6.2),dpi=300)
xb=np.arange(len(CLASSES)); bw=0.26
for k,X in enumerate(STR):
    sub=df[df.strain==X].set_index("cls").reindex(CLASSES)
    bars=ax.bar(xb+(k-1)*bw,sub.fold,bw,color=COL[X],label=f"{X.replace('_','/')} (uniquely-mapping: {int(sub.n_uniq.sum()):,})")
    for x,v in zip(xb+(k-1)*bw,sub.fold):
        if pd.notna(v): ax.text(x,v*1.04,f"{v:.1f}×",ha="center",va="bottom",fontsize=7,color=COL[X],fontweight="bold")
ax.axhline(1,color="#555",ls="--",lw=1); ax.text(len(CLASSES)-1,1.05,"random-locus expectation (1×)",ha="right",fontsize=7,color="#555")
ax.set_yscale("log"); ax.set_xticks(xb); ax.set_xticklabels(CLASSES,fontsize=8.5)
ax.set_ylabel("fold-enrichment: piRNA locus inside a lncRNA gene\n(observed ÷ genomic expectation)",fontsize=9)
ax.set_title("ncRNA-driven test (protein-coding-EXCLUDED, confounding-free): modest, largely class-independent lncRNA enrichment (~0.9–2.4×)",fontsize=10,fontweight="bold")
ax.legend(fontsize=8,frameon=False,loc="upper right"); ax.spines[['top','right']].set_visible(False)
fig.text(0.5,0.005,"Confounding-free: clean-lncRNA space = lncRNA MINUS protein_coding (drops antisense/intronic lncRNAs over coding genes, e.g. Gm52992/Dlg2; ~21% of annotated lncRNA bp removed). Random expectation = clean-lncRNA bp ÷ genome bp. "
  "Unlike the TE-driven test (strongly class-specific), lncRNA-overlap does NOT separate the classes — all modest (~0.9–2.4×); in SPRET strain-private is highest (2.38×) while common is ~1× (0.88×). Uniquely-mapping loci only.",ha="center",fontsize=6.0,color="#555")
fig.tight_layout(rect=[0,0.02,1,1])
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_ncrna_driven_test.{e}",bbox_inches="tight")
print("wrote Fig_ncrna_driven_test.{png,pdf,svg} + SourceData_ncrna_driven_test.csv")
