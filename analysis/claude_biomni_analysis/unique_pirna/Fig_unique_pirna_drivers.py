#!/usr/bin/env python3
"""What DRIVES each strain's UNIQUE (strain-private-locus) piRNAs — TE vs ncRNA(lncRNA) vs protein-coding gene vs
intergenic? Per strain, each strain-private piRNA locus (cand_self16 BAM, own genome) is intersected with the
per-strain RepeatMasker BED (TE) and per-strain Ensembl gene models (protein_coding gene / lncRNA gene); each locus
gets ONE mutually-exclusive driver by priority TE > lncRNA > protein-coding > intergenic. Keeps timepoint (id =
strain|tp|seq). Panels (canonical order, wild red): (A) driver per strain (pooled tp); (B) driver per timepoint
(pooled strains); (C) STRAIN × TIMEPOINT together — per-strain driver composition faceted by window, so each
strain's prepachytene(TE)->pachytene(lncRNA/intergenic) switch is visible. Cache-aware (per-locus driver table ->
instant re-plot). Biology framing queued for BioMNI. Linked: [[project_te_driven_finding]] [[project_ncrna_driven_finding]]."""
import warnings; warnings.filterwarnings("ignore")
import sys, subprocess, tempfile, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
import pandas as pd, numpy as np, pysam
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; GB = f"{PG}/gene_beds"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPS = ["E16.5", "P12.5", "P20.5"]
DRIVERS = [("TE", "#1B7837"), ("lncRNA", "#762A83"), ("protein-coding", "#4575B4"), ("intergenic", "#BBBBBB")]
DCOLS = [k for k, _ in DRIVERS]
cache = f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/data/source_data/SourceData_unique_pirna_drivers.csv"
if os.path.exists(cache):
    T = pd.read_csv(cache); print("loaded cached driver table")
else:
    d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz")   # mm0-clean strain-private (klass5)
    def ids_overlapping(locibed, bed):
        if not os.path.exists(bed) or os.path.getsize(bed) == 0: return set()
        out = subprocess.run(f"sort -k1,1 -k2,2n {locibed} | {BT} intersect -a - -b {bed} -u", shell=True, capture_output=True, text=True).stdout
        return {ln.split('\t')[3] for ln in out.splitlines() if ln}
    KU = {"unique: conserved-but-silent": "conserved-but-silent", "unique: strain-private locus": "strain-private"}
    rows = []
    for X in CANON:
        g = d[(d.strain == X) & (d.klass5.isin(KU))].copy(); g["id"] = X + "|" + g.timepoint + "|" + g.sequence
        id2cls = dict(zip(g.id, g.klass5.map(KU)))
        sp = set(g.id); rm = f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
        if not sp: continue
        bam = pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam", "rb")
        tb = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=PG); seen = {}
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.query_name not in sp or a.query_name in seen: continue
            seen[a.query_name] = 1; tb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
        tb.close(); bam.close()
        TEset = ids_overlapping(tb.name, rm); LNCset = ids_overlapping(tb.name, f"{GB}/{X}.lnc.bed"); PCset = ids_overlapping(tb.name, f"{GB}/{X}.pc.bed")
        os.unlink(tb.name)
        for locid in seen:
            tp = TPMAP.get(locid.split("|")[1], locid.split("|")[1])
            drv = "TE" if locid in TEset else "lncRNA" if locid in LNCset else "protein-coding" if locid in PCset else "intergenic"
            rows.append((X, tp, id2cls[locid], drv))
        print(f"{X}: n={len(seen):,} TE={len(TEset):,} lnc={len(LNCset):,} pc={len(PCset):,}")
    T = pd.DataFrame(rows, columns=["strain", "tp", "klass", "driver"]); T.to_csv(cache, index=False)
def comp(df, idx):  # fraction matrix
    m = df.groupby([idx, "driver"]).size().unstack(fill_value=0).reindex(columns=DCOLS).fillna(0)
    return m.div(m.sum(1).replace(0, np.nan), axis=0) * 100
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
KLS = [("strain-private", "unique: strain-private — locus NEW", "#7a3b9a"), ("conserved-but-silent", "unique: conserved-but-silent — locus SHARED", "#0072B2")]
fig = plt.figure(figsize=(16.5, 9.2), dpi=300)
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.62, wspace=0.22)
def stack(ax, frac, cats, labels, wild_red=True, w=0.82):
    xx = np.arange(len(cats)); bottom = np.zeros(len(cats))
    for k, col in DRIVERS:
        ax.bar(xx, frac[k].values, bottom=bottom, width=w, color=col, edgecolor="white", linewidth=0.3)
        bottom = bottom + np.nan_to_num(frac[k].values)
    ax.set_xticks(xx); ax.set_xticklabels(labels, rotation=90, fontsize=6.5); ax.set_ylim(0, 100)
    if wild_red:
        for t, X in zip(ax.get_xticklabels(), cats):
            if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
for ri, (kl, ktitle, kc) in enumerate(KLS):
    Tk = T[T.klass == kl]
    axA = fig.add_subplot(gs[ri, 0:2]); axB = fig.add_subplot(gs[ri, 2])
    stack(axA, comp(Tk, "strain").reindex(CANON), CANON, [s.replace("_", "/") for s in CANON])
    axA.set_ylabel("driver composition (%)", fontsize=9)
    axA.set_title(f"{'AB'[ri]}  {ktitle}  (n={len(Tk):,}) — driver per strain (pooled timepoints)", fontsize=9.3, fontweight="bold", loc="left", color=kc)
    stack(axB, comp(Tk, "tp").reindex(TPS), TPS, TPS, wild_red=False, w=0.7); axB.set_ylim(0, 100)
    axB.set_title("by timepoint (pooled strains)", fontsize=8.6, fontweight="bold", loc="left")
fig.legend(handles=[Patch(facecolor=c, label=k) for k, c in DRIVERS], ncol=4, fontsize=9.5, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.005))
fig.suptitle("What drives the two genuinely-unique subcategories — strain-private (locus NEW) vs conserved-but-silent (locus SHARED)\n"
             "TE > lncRNA > protein-coding > intergenic; per strain (canonical order, wild red) + per timepoint", fontsize=11, fontweight="bold", y=1.0)
axA.set_xticklabels([])   # bottom-row strain labels carried by the classical/wild companion below
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
# classical(blue)/wild(orange) companion below the bottom composition row: candidates per strain (subspecies colour scheme)
fig.subplots_adjust(bottom=0.20)
_tot=Tk.groupby("strain").size().reindex(CANON).fillna(0).values   # Tk = bottom-row class subset after loop
_cax=add_classical_wild_companion(fig, axA, CANON, _tot, gap=0.085, height_frac=0.22, ylabel="n per\nstrain (log)")
_cax.set_xticks(np.arange(len(CANON))); _cax.set_xticklabels([s.replace("_","/") for s in CANON], rotation=90, fontsize=6.5)
for lab,s in zip(_cax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title(f"classical (blue) vs wild-derived (orange) — {ktitle} candidates per strain", fontsize=7.5, fontweight="bold", loc="left")
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_unique_pirna_drivers.{e}", bbox_inches="tight")
print("wrote Fig_unique_pirna_drivers (split: strain-private vs conserved-but-silent)")
print(T.groupby(["klass", "driver"]).size().unstack(fill_value=0).to_string())
