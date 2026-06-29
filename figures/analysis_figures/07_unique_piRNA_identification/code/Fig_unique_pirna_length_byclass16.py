#!/usr/bin/env python3
"""piRNA length distribution per klass5 classification class, with a LINE per developmental timepoint (16-strain,
≥2-read). The 5 classes are shown SIDE BY SIDE (single row, shared y-scale) for direct visual comparison. For each
class the length distribution (% of that class's distinct candidate sequences at each length 20-36 nt) is drawn as a
line per timepoint (E16.5 / P12.5 / P20.5). Shaded band = data-driven piRNA window (26-30 nt, FWHM, mode 27);
vertical dashed guides at 27 nt (prepachytene, blue) & 30 nt (pachytene, red). Source: unique16/final_classified_clean_2read.csv.gz."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
KL5 = ["expressed elsewhere (exact)", "SNP-variant (1-3mm)", "low-quality: no mm0 own-genome locus", "unique: conserved-but-silent", "unique: strain-private locus"]
KLAB = ["expressed-elsewhere", "SNP-variant", "low-quality", "conserved-but-silent\n(locus SHARED)", "strain-private\n(locus NEW)"]
KCOL = ["#9e9e9e", "#E69F00", "#cdb892", "#0072B2", "#7a3b9a"]
TP = ["16.5dpc", "12.5dpp", "20.5dpp"]; TPLAB = ["E16.5", "P12.5", "P20.5"]; TPCOL = ["#4393C3", "#E8852B", "#B2182B"]
LEN = list(range(20, 37))
d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz", usecols=["sequence", "timepoint", "length", "klass5"])
d = d[(d.length >= 20) & (d.length <= 36)]
plt.rcParams.update({"font.family": "Liberation Sans"})
fig, axes = plt.subplots(1, 5, figsize=(19, 4.9), dpi=300)   # 5 classes side by side
src = []; gmax = 0.0
for ci, (k, klab, kc) in enumerate(zip(KL5, KLAB, KCOL)):
    ax = axes[ci]; sub = d[d.klass5 == k]
    ax.axvspan(25.5, 30.5, color="#FFF3D6", zorder=0)
    for t, tl, tc in zip(TP, TPLAB, TPCOL):
        st = sub[sub.timepoint == t].drop_duplicates("sequence")
        if len(st) < 5: continue
        h = st.length.value_counts().reindex(LEN, fill_value=0).astype(float); h = 100 * h / h.sum()
        gmax = max(gmax, float(h.max())); mode = int(h.idxmax())
        ax.plot(LEN, h.values, color=tc, lw=1.7, marker="o", ms=2.4, label=f"{tl} (n={len(st):,}, mode {mode})")
        for L, v in zip(LEN, h.values): src.append((klab.replace(chr(10), " "), tl, L, round(v, 3), len(st)))
    for xv, xc in ((27, "#3b6ea5"), (30, "#B2182B")): ax.axvline(xv, ls=(0, (4, 3)), lw=1.1, color=xc, alpha=0.85, zorder=2)
    ax.set_title(klab, fontsize=9, fontweight="bold", color=kc, loc="left")
    ax.set_xlabel("piRNA length (nt)", fontsize=8)
    ax.legend(fontsize=5.7, frameon=False, loc="upper right")
    ax.tick_params(labelsize=7); ax.set_xticks(range(20, 37, 3)); ax.spines[["top", "right"]].set_visible(False)
YT = gmax * 1.12
for ci in range(5):
    ax = axes[ci]; ax.set_ylim(0, YT)
    for xv, xc in ((27, "#3b6ea5"), (30, "#B2182B")): ax.text(xv, gmax * 1.015, str(xv), ha="center", va="bottom", fontsize=7, color=xc, fontweight="bold")
axes[0].set_ylabel("% of distinct sequences", fontsize=8.5)
fig.suptitle("piRNA length distribution per classification class, by developmental timepoint (16-strain, ≥2-read) — five classes side by side for comparison", fontsize=12, fontweight="bold", y=1.0)
fig.text(0.5, -0.02, "One line per timepoint (E16.5 / P12.5 / P20.5); shaded band = 26–30 nt piRNA window (FWHM, mode 27); vertical dashed guides at 27 nt (prepachytene, blue) & 30 nt (pachytene, red). "
  "% = distinct candidate sequences of that class at that timepoint (≥2-read adopted); n and mode per timepoint shown in each panel's legend. Shared y-scale across panels.", ha="center", fontsize=7.2, color="#555")
fig.tight_layout(rect=[0, 0, 1, 0.95])
pd.DataFrame(src, columns=["class", "timepoint", "length_nt", "pct_of_sequences", "n_sequences"]).to_csv(f"{U}/pangenome_te/SourceData_unique_pirna_length_byclass16.csv", index=False)
for e in ("pdf", "svg", "png"): fig.savefig(f"{U}/pangenome_te/Fig_unique_pirna_length_byclass16.{e}", bbox_inches="tight")
print("wrote Fig_unique_pirna_length_byclass16 (5 classes side by side, shared y-scale)")
print(d.groupby("klass5").length.median())
