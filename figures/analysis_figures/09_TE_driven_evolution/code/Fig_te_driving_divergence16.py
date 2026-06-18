#!/usr/bin/env python3
"""Is TE driving the DIVERGENCE class? (a) divergence loci ARE TE-derived, enriched for active retrotransposons
(ERVK/L1) vs baseline clusters [from divergence_te_overlap.py, n=200 each]; (b) TE-family composition; (c,d) cis
vs trans test [halSnps, divergence_cis_vs_trans.csv]: carrier-vs-silent SNP rate AT the locus vs a matched random-
region baseline → the strain-restriction is MOSTLY REGULATORY (trans; sequence as conserved as genome average),
with a MINORITY of cis TE-sequence divergence. TEs are the substrate, not (mostly) the cis cause."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
DIV, BASE = "#6a3d9a", "#9e9e9e"; C_REG, C_CIS = "#1B7837", "#C0392B"
# (a,b) TE enrichment + family composition — from divergence_te_overlap.py (n=200 divergence, 200 baseline)
driver_pct = {"divergence": 98, "baseline": 85}; OR, ORp = 8.65, 2.5e-6
fam = {"LTR/ERVK": (87, 28), "LINE/L1": (45, 28), "LTR/ERVL-MaLR": (27, 33), "LTR/ERV1": (20, 9), "SINE/B2": (14, 68)}   # (divergence, baseline) dominant-family counts /200
# (c,d) cis vs trans
ct = pd.read_csv(f"{PG}/divergence_cis_vs_trans.csv"); ct["ratio"] = ct.locus_rate / ct.base_rate.replace(0, np.nan); ct = ct.dropna(subset=["ratio"])
reg = (ct.locus_rate <= 1.5 * ct.base_rate).mean() * 100; cis = (ct.locus_rate > 2 * ct.base_rate).mean() * 100; inter = 100 - reg - cis
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8})
fig, ax = plt.subplots(2, 2, figsize=(11, 8.6), dpi=300)
def nostpr(a): a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=8)
# a
a = ax[0, 0]; a.bar([0, 1], [driver_pct["divergence"], driver_pct["baseline"]], color=[DIV, BASE], width=0.6, edgecolor="white")
for i, k in enumerate(["divergence", "baseline"]): a.text(i, driver_pct[k] + 1, f"{driver_pct[k]}%", ha="center", fontsize=9, fontweight="bold")
a.set_xticks([0, 1]); a.set_xticklabels(["divergence\nloci", "random\nclusters"], fontsize=8.5); a.set_ylim(0, 110); a.set_ylabel("driver-TE overlap (%)", fontsize=9); nostpr(a)
a.set_title(f"a   Divergence loci ARE TE-derived (active retrotransposons)\nenrichment OR = {OR:.1f}, Fisher P = {ORp:.0e}", fontsize=9, fontweight="bold", loc="left")
# b family composition
b = ax[0, 1]; fams = list(fam); xb = np.arange(len(fams)); w = 0.38
b.bar(xb - w / 2, [fam[f][0] / 2 for f in fams], w, color=DIV, label="divergence", edgecolor="white")
b.bar(xb + w / 2, [fam[f][1] / 2 for f in fams], w, color=BASE, label="baseline", edgecolor="white")
b.set_xticks(xb); b.set_xticklabels(fams, rotation=30, ha="right", fontsize=7.5); b.set_ylabel("loci (%)", fontsize=9); b.legend(fontsize=8, frameon=False); nostpr(b)
b.set_title("b   ERVK/L1 enriched in divergence loci; SINE/B2 in baseline", fontsize=9, fontweight="bold", loc="left")
# c scatter cis vs trans
c = ax[1, 0]; col = [C_CIS if r > 2 else ("#888" if r > 1.5 else C_REG) for r in ct.ratio]
c.scatter(ct.base_rate * 100, ct.locus_rate * 100, s=34, c=col, edgecolor="white", linewidth=0.5, zorder=3)
mx = max(ct.locus_rate.max(), ct.base_rate.max()) * 100 * 1.1
c.plot([0, mx], [0, mx], "--", color="#aaa", lw=0.9, label="locus = baseline (trans)"); c.plot([0, mx / 2], [0, mx], ":", color=C_CIS, lw=0.9, label="locus = 2× baseline (cis)")
c.set_xlabel("genome baseline SNP rate, carrier vs silent (%)", fontsize=9); c.set_ylabel("LOCUS SNP rate, carrier vs silent (%)", fontsize=9); nostpr(c); c.legend(fontsize=7, frameon=False, loc="upper left")
c.set_title(f"c   cis vs trans: locus vs genome-average divergence\nmedian ratio {ct.ratio.median():.2f}× (≈ genome average → trans)", fontsize=9, fontweight="bold", loc="left")
# d summary
d = ax[1, 1]; d.bar([0, 1, 2], [reg, inter, cis], color=[C_REG, "#888", C_CIS], width=0.62, edgecolor="white")
for i, v in enumerate([reg, inter, cis]): d.text(i, v + 1, f"{v:.0f}%", ha="center", fontsize=9, fontweight="bold")
d.set_xticks([0, 1, 2]); d.set_xticklabels(["REGULATORY\n(trans, ≤1.5×)", "intermediate\n(1.5–2×)", "cis-sequence\n(>2×)"], fontsize=8); d.set_ylim(0, max(reg, cis) * 1.25); d.set_ylabel("% of divergence loci", fontsize=9); nostpr(d)
d.set_title(f"d   Divergence is mostly REGULATORY ({reg:.0f}%); minority cis ({cis:.0f}%)", fontsize=9, fontweight="bold", loc="left")
fig.suptitle("Is TE driving the divergence? — TEs are the SUBSTRATE (ERVK/L1, OR=8.6) but the strain-restricted expression is mostly REGULATORY (trans), not cis TE-sequence divergence", fontsize=10.2, fontweight="bold", y=1.005)
fig.tight_layout()
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_te_driving_divergence16.{e}", bbox_inches="tight")
print(f"wrote Fig_te_driving_divergence16.{{png,pdf,svg}} | regulatory={reg:.0f}% intermediate={inter:.0f}% cis={cis:.0f}% | n_loci={len(ct)}")
