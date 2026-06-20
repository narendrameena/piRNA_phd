#!/usr/bin/env python3
"""16-strain TE families of origin for the TWO genuinely-unique subcategories (klass5 ≥2-read): strain-private
(locus NEW) vs conserved-but-silent (locus SHARED). For each strain, intersect each subcategory's piRNA loci
(cand_self16 BAM, own genome) with the per-strain RepeatMasker BED (col4 = name|class/family), assign each locus
its largest-overlap TE family, count. (A) strain-private TE-family heatmap; (B) conserved-but-silent TE-family
heatmap; (C) TE-derived fraction per strain, grouped (strain-private vs CBS); (D) TE-family small-RNA expression
(sRNA-on-TE). CAVEAT: index = main chr + MT only -> TE fraction is a LOWER BOUND. Cache: SourceData_TE_private_families16_byclass.csv."""
import warnings; warnings.filterwarnings("ignore")
import sys, subprocess, tempfile, os; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np, pysam
from collections import Counter, defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mc
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
SUB = [("strain-private", "unique: strain-private locus"), ("conserved-but-silent", "unique: conserved-but-silent")]
cache = f"{PG}/SourceData_TE_private_families16_byclass.csv"
if os.path.exists(cache):
    T = pd.read_csv(cache); print("loaded cached by-class TE-family table")
else:
    d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")
    rows = []   # (klass, family, strain, count) + special families "__nseen__"/"__nte__" for the fraction
    for X in CANON:
        rm = f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
        if not os.path.exists(rm): continue
        g = d[(d.strain == X) & (d.klass5.isin([k for _, k in SUB]))]
        cls = {X + "|" + r.timepoint + "|" + r.sequence: ("strain-private" if r.klass5.endswith("private locus") else "conserved-but-silent") for r in g.itertuples()}
        if not cls: continue
        bam = pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam", "rb")
        tb = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=PG); seen = set()
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.query_name not in cls or a.query_name in seen: continue
            seen.add(a.query_name); tb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
        tb.close(); bam.close()
        out = subprocess.run(f"sort -k1,1 -k2,2n {tb.name} | {BT} intersect -a - -b {rm} -wa -wb", shell=True, capture_output=True, text=True).stdout
        os.unlink(tb.name)
        best = {}
        for ln in out.splitlines():
            f = ln.split("\t")
            if len(f) < 8: continue
            ov = min(int(f[2]), int(f[6])) - max(int(f[1]), int(f[5])); fam = f[7].split("|")[-1] if "|" in f[7] else f[7]
            if f[3] not in best or ov > best[f[3]][0]: best[f[3]] = (ov, fam)
        cnt = defaultdict(Counter); nseen = Counter(); nte = Counter()
        for cid in seen: nseen[cls[cid]] += 1
        for cid, (ov, fam) in best.items(): cnt[cls[cid]][fam] += 1; nte[cls[cid]] += 1
        for klab in ("strain-private", "conserved-but-silent"):
            for fam, n in cnt[klab].items(): rows.append((klab, fam, X, n))
            rows.append((klab, "__nseen__", X, nseen[klab])); rows.append((klab, "__nte__", X, nte[klab]))
        print(f"{X}: strain-private={nseen['strain-private']:,} conserved-but-silent={nseen['conserved-but-silent']:,}")
    T = pd.DataFrame(rows, columns=["klass", "family", "strain", "count"]); T.to_csv(cache, index=False)
# ---- build matrices ----
def fammat(klab, topn=14):
    sub = T[(T.klass == klab) & (~T.family.isin(["__nseen__", "__nte__"]))]
    M = sub.pivot_table(index="family", columns="strain", values="count", aggfunc="sum").reindex(columns=CANON).fillna(0)
    return M.loc[M.sum(1).sort_values(ascending=False).head(topn).index]
def tefrac(klab):
    a = T[(T.klass == klab) & (T.family == "__nte__")].set_index("strain")["count"]; b = T[(T.klass == klab) & (T.family == "__nseen__")].set_index("strain")["count"]
    return {X: (100 * a.get(X, 0) / b.get(X, np.nan) if b.get(X, 0) else np.nan) for X in CANON}
Msp, Mcbs = fammat("strain-private"), fammat("conserved-but-silent"); fr_sp, fr_cbs = tefrac("strain-private"), tefrac("conserved-but-silent")
# ---- sRNA TE-family expression (reuse cached expression matrix, family x strain) ----
E_csv = f"{PG}/SourceData_TE_private_families16_expression.csv"
E = pd.read_csv(E_csv, index_col=0).reindex(index=Msp.index, columns=CANON).fillna(0.0) if os.path.exists(E_csv) else pd.DataFrame(0.0, index=Msp.index, columns=CANON)
plt.rcParams.update({"font.family": "Liberation Sans"})
CMAP = mc.LinearSegmentedColormap.from_list("vivBlue", ["#eaf3fb", "#9ecae8", "#3a8fd4", "#1565a8"]); CMAP.set_bad("white")
def _tc(n): return "white" if n > 0.55 else "#222"
fig = plt.figure(figsize=(17, 11), dpi=300); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.42, wspace=0.16)
def heat(ax, M, title, fs_lab=4.6):
    L = np.log10(M.values + 1); vmax = L.max() if L.max() > 0 else 1
    im = ax.imshow(np.ma.masked_where(M.values == 0, L), aspect="auto", cmap=CMAP, vmin=0, vmax=vmax)
    ax.set_xticks(range(len(CANON))); ax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=7)
    ax.set_xticks(np.arange(-.5, len(CANON), 1), minor=True); ax.set_yticks(np.arange(-.5, len(M), 1), minor=True)
    ax.grid(which="minor", color="#e6e6e6", linewidth=0.5); ax.tick_params(which="minor", length=0)
    for i in range(len(M)):
        for j in range(len(CANON)):
            v = int(M.values[i, j])
            if v > 0: ax.text(j, i, f"{v:,}" if v < 1000 else f"{v//1000}k", ha="center", va="center", fontsize=fs_lab, color=_tc(np.log10(v+1)/vmax if vmax else 0))
    ax.set_title(title, fontsize=9.4, fontweight="bold", loc="left"); return im
imA = heat(fig.add_subplot(gs[0, 0]), Msp, "A  strain-private (locus NEW) — TE family of origin (counts, log)")
imB = heat(fig.add_subplot(gs[0, 1]), Mcbs, "B  conserved-but-silent (locus SHARED) — TE family of origin (counts, log)")
axC = fig.add_subplot(gs[1, 0]); x = np.arange(len(CANON)); ww = 0.4
axC.bar(x - ww/2, [fr_sp[s] for s in CANON], ww, color="#7a3b9a", label="strain-private")
axC.bar(x + ww/2, [fr_cbs[s] for s in CANON], ww, color="#0072B2", label="conserved-but-silent")
axC.set_xticks(x); axC.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7); axC.set_ylabel("% TE-derived (lower bound)", fontsize=8.5)
axC.legend(fontsize=8, frameon=False); axC.spines[["top", "right"]].set_visible(False)
axC.set_title("C  TE-derived fraction per strain — strain-private vs conserved-but-silent (mappable loci; lower bound)", fontsize=9, fontweight="bold", loc="left")
axE = fig.add_subplot(gs[1, 1]); LE = np.log10(E.values + 1); vmn = LE[E.values > 0].min() if (E.values > 0).any() else 0; vmx = LE.max() if LE.max() > 0 else 1
imE = axE.imshow(np.ma.masked_where(E.values <= 0, LE), aspect="auto", cmap=CMAP, vmin=vmn, vmax=vmx)
axE.set_xticks(range(len(CANON))); axE.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7)
axE.set_yticks(range(len(E))); axE.set_yticklabels(E.index, fontsize=7)
for i in range(len(E)):
    for j in range(len(CANON)):
        v = E.values[i, j]
        if v > 0: axE.text(j, i, f"{v/1e6:.1f}M" if v >= 1e6 else (f"{v/1e3:.0f}k" if v >= 1e3 else f"{v:.0f}"), ha="center", va="center", fontsize=3.8, color=_tc((np.log10(v+1)-vmn)/(vmx-vmn) if vmx > vmn else 0.5))
axE.set_title("D  TE-family small-RNA expression (sRNA-on-TE, summed, log) — strain-private families", fontsize=9.2, fontweight="bold", loc="left")
fig.colorbar(imA, ax=fig.axes[0], fraction=0.025, pad=0.01).ax.tick_params(labelsize=6)
fig.suptitle("TE families driving the two genuinely-unique subcategories across 16 strains — strain-private (A) vs conserved-but-silent (B); LTR/ERVK(IAP)+LINE/L1 dominate strain-private, more so in wild-derived", fontsize=10, fontweight="bold", y=1.0)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_TE_private_families16.{e}", bbox_inches="tight")
print("wrote Fig_TE_private_families16 (by subcategory)")
