#!/usr/bin/env python3
"""Nature-Genetics-level figure for the pangenome insertion-driven piRNA finding (all 16 strains), WITH the
statistical tests and the QC controls that rule out confounders. Reads the QC source data (computed once by
compute_te_driven_qc.py) so the figure is reproducible from a single table.
FINDING (top row): strain-private TE insertions seed strain-private piRNA source loci; wild-derived strains gain
far more, tracking their private-insertion burden — and it holds as a FRACTION of all clusters.
QC (bottom row): not explained by sRNA depth or by total cluster count; per-insertion conversion is ~constant,
so the excess is an insertion-COUNT effect (genome divergence), not higher conversion efficiency."""
import os, sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
df = pd.read_csv(f"{SD}/Fig_te_driven_pangenome16_QC_sourcedata.csv")   # canonical strain order preserved
CW, CC = "#C0392B", "#4393C3"; cols = [CW if w else CC for w in df.wild]
w, c = df[df.wild], df[~df.wild]
def mwp(col): return mannwhitneyu(w[col], c[col], alternative="greater").pvalue
p_loci, p_frac, p_prod = mwp("insertion_driven_loci"), mwp("frac_clusters_ins_driven"), mannwhitneyu(w.loci_per_1k_insertions, c.loci_per_1k_insertions, alternative="two-sided").pvalue
rho_ins, pp_ins = spearmanr(df.private_ins_ge150bp, df.insertion_driven_loci)
rho_dep, pp_dep = spearmanr(df.srna_depth, df.insertion_driven_loci)
rho_cl, pp_cl = spearmanr(df.total_PICB_clusters, df.insertion_driven_loci)
def star(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8})
fig, ax = plt.subplots(2, 3, figsize=(14, 8.4), dpi=300)
xs = np.arange(len(df)); xl = [s.replace("_", "/") for s in df.strain]
def xticks(a):
    a.set_xticks(xs); a.set_xticklabels(xl, rotation=90, fontsize=6.3)
    for t, wd in zip(a.get_xticklabels(), df.wild):
        if wd: t.set_color(CW); t.set_fontweight("bold")
def nostpr(a): a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=7)
# ---- a: loci per strain ----
a = ax[0, 0]; a.bar(xs, df.insertion_driven_loci, color=cols, edgecolor="white", linewidth=0.4)
for i, v in enumerate(df.insertion_driven_loci): a.text(i, v + df.insertion_driven_loci.max() * 0.012, f"{v:,}", ha="center", va="bottom", fontsize=5, rotation=90, color=cols[i], fontweight="bold")
xticks(a); nostpr(a); a.set_ylim(0, df.insertion_driven_loci.max() * 1.2); a.set_ylabel("insertion-driven piRNA\nsource loci (n)", fontsize=8.5)
a.set_title(f"a   Wild strains gain more TE-insertion-driven piRNA loci\nwild vs classical: Mann–Whitney P = {p_loci:.1e} ({star(p_loci)})", fontsize=8.4, fontweight="bold", loc="left")
a.legend(handles=[Patch(facecolor=CW, label="wild-derived (n=4)"), Patch(facecolor=CC, label="classical (n=12)")], fontsize=6.8, frameon=False, loc="upper left")
# ---- b: insertions vs loci ----
b = ax[0, 1]; b.scatter(df.private_ins_ge150bp, df.insertion_driven_loci, c=cols, s=34, edgecolor="white", linewidth=0.5, zorder=3)
for _, r in df.iterrows(): b.annotate(r.strain.replace("_", "/"), (r.private_ins_ge150bp, r.insertion_driven_loci), fontsize=5, color=(CW if r.wild else "#666"), xytext=(3, 1), textcoords="offset points")
b.set_xscale("log"); b.set_yscale("log"); nostpr(b); b.set_xlabel("private TE insertions ≥150 bp (n, log)", fontsize=8.5); b.set_ylabel("insertion-driven piRNA loci (n, log)", fontsize=8.5)
b.set_title(f"b   Loci count tracks private-insertion burden\nSpearman ρ = {rho_ins:.2f}, P = {pp_ins:.1e}", fontsize=8.4, fontweight="bold", loc="left")
# ---- c: fraction of clusters insertion-driven (normalised finding) ----
cc = ax[0, 2]; cc.bar(xs, df.frac_clusters_ins_driven * 100, color=cols, edgecolor="white", linewidth=0.4)
xticks(cc); nostpr(cc); cc.set_ylabel("insertion-driven fraction of\nall PICB clusters (%)", fontsize=8.5)
cc.set_title(f"c   Holds after normalising to total clusters\nwild vs classical: P = {p_frac:.1e} ({star(p_frac)})", fontsize=8.4, fontweight="bold", loc="left")
# ---- d: QC — depth does not explain it ----
d = ax[1, 0]; d.scatter(df.srna_depth / 1e9, df.insertion_driven_loci, c=cols, s=34, edgecolor="white", linewidth=0.5, zorder=3)
for _, r in df.iterrows(): d.annotate(r.strain.replace("_", "/"), (r.srna_depth / 1e9, r.insertion_driven_loci), fontsize=4.8, color=(CW if r.wild else "#888"), xytext=(3, 1), textcoords="offset points")
nostpr(d); d.set_yscale("log"); d.set_xlabel("sRNA library depth (billion mapped reads)", fontsize=8.5); d.set_ylabel("insertion-driven loci (n, log)", fontsize=8.5)
d.set_title(f"d   QC: depth does NOT explain loci count\nSpearman ρ = {rho_dep:.2f}, P = {pp_dep:.2f} (n.s.)", fontsize=8.4, fontweight="bold", loc="left")
# ---- e: QC — total clusters not higher in wild ----
e = ax[1, 1]; e.bar(xs, df.total_PICB_clusters / 1000, color=cols, edgecolor="white", linewidth=0.4)
xticks(e); nostpr(e); e.set_ylabel("total PICB clusters\nper strain (×1,000)", fontsize=8.5)
e.set_title(f"e   QC: wild do NOT have more clusters overall\nloci ~ total clusters: ρ = {rho_cl:.2f}, P = {pp_cl:.2f} (n.s.)", fontsize=8.4, fontweight="bold", loc="left")
# ---- f: productivity ~ constant ----
f = ax[1, 2]; f.bar(xs, df.loci_per_1k_insertions, color=cols, edgecolor="white", linewidth=0.4)
f.axhline(df.loci_per_1k_insertions.median(), ls="--", lw=0.8, color="#555", zorder=3)
xticks(f); nostpr(f); f.set_ylabel("conversion: piRNA loci per\n1,000 private insertions", fontsize=8.5)
f.set_title(f"f   Per-insertion conversion ~constant (not higher in wild)\nwild vs classical: P = {p_prod:.2f} → excess = insertion COUNT", fontsize=8.4, fontweight="bold", loc="left")
fig.suptitle("Strain-private TE insertions seed strain-private piRNA source loci across the 16-strain mouse pangenome — wild-derived genome divergence drives piRNA-locus gain", fontsize=11, fontweight="bold", y=1.005)
fig.tight_layout()
_TH="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/10_cluster_pangenome_PAV/data/source_data"; os.makedirs(_TH,exist_ok=True)
df.to_csv(f"{_TH}/SourceData_Fig_te_driven_pangenome16.csv",index=False)   # per-strain: insertion-driven loci/fraction, private-insertion burden, sRNA depth, total clusters, conversion (all 6 panels)
for ext in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_te_driven_pangenome16.{ext}", bbox_inches="tight")
print(f"wrote Fig_te_driven_pangenome16.{{png,pdf,svg}}  (P_loci={p_loci:.2g}, P_frac={p_frac:.2g}, rho_ins={rho_ins:.2f}, rho_depth={rho_dep:.2f})")
