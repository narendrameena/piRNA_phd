#!/bin/bash
#SBATCH --job-name=PICB_cmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_comparison_results/slurm_%j.log
#SBATCH --error=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_comparison_results/slurm_%j.err

# PICB comparison: chr-by-chr vs whole-BAM
# Sample: C57BL/6 GRCm38 E15.5 rep1
# Both runs use the same BAM (Aligned.sortedByCoord.out.bam) and parameters (25-36 nt).
# The chr-by-chr run should reproduce the existing 328-cluster result.
# The whole-BAM run should produce an identical result.

set -euo pipefail

export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"

# ── Paths ────────────────────────────────────────────────────────────────────
BASE=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA

FASTA=${BASE}/resources/black6/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa
BAM=${BASE}/resources/black6/star_aligned/GRCm38/E15_5_rep1/Aligned.sortedByCoord.out.bam

CHUNKED_R=${BASE}/workflow/scripts/R/picb_script_chunked.R
WHOLEBAM_R=${BASE}/workflow/scripts/R/picb_script_wholebam.R
FIGURE_PY=${BASE}/analysis/claude_biomni_analysis/picb_comparison_figure.py

OUT_DIR=${BASE}/analysis/claude_biomni_analysis/picb_comparison_results
OUT_CHRBYCHR=${OUT_DIR}/E15_5_rep1_chrbychr.xlsx
OUT_WHOLEBAM=${OUT_DIR}/E15_5_rep1_wholebam.xlsx

PIRNA_MIN=25
PIRNA_MAX=36

mkdir -p "$OUT_DIR"

ts() { echo "[$(date '+%Y-%m-%d %H:%M:%S')]"; }

echo "$(ts) ============================================================"
echo "$(ts) PICB comparison: chr-by-chr vs whole-BAM"
echo "$(ts) FASTA:  $FASTA"
echo "$(ts) BAM:    $BAM  ($(du -h "$BAM" | cut -f1))"
echo "$(ts) Output: $OUT_DIR"
echo "$(ts) piRNA:  ${PIRNA_MIN}-${PIRNA_MAX} nt"
echo "$(ts) ============================================================"

# ── Validate inputs ──────────────────────────────────────────────────────────
for f in "$FASTA" "$BAM" "$CHUNKED_R" "$WHOLEBAM_R" "$FIGURE_PY"; do
    if [ ! -f "$f" ]; then
        echo "$(ts) ERROR: required file not found: $f"
        exit 1
    fi
done
echo "$(ts) All input files confirmed."

# ── RUN 1: Chr-by-chr ────────────────────────────────────────────────────────
echo ""
echo "$(ts) ============================================================"
echo "$(ts) RUN 1: Chr-by-chr (picb_script_chunked.R)"
echo "$(ts) ============================================================"

LOG_CHRBYCHR=${OUT_DIR}/E15_5_rep1_chrbychr.log

conda run -n snakemake Rscript "$CHUNKED_R" \
    "$FASTA" "$BAM" "$OUT_CHRBYCHR" "$PIRNA_MIN" "$PIRNA_MAX" \
    2>&1 | tee "$LOG_CHRBYCHR"

if [ ! -f "$OUT_CHRBYCHR" ] || [ ! -s "$OUT_CHRBYCHR" ]; then
    echo "$(ts) ERROR: chr-by-chr output not created or empty."
    exit 1
fi
echo "$(ts) Chr-by-chr complete: $OUT_CHRBYCHR"

# ── RUN 2: Whole-BAM ─────────────────────────────────────────────────────────
echo ""
echo "$(ts) ============================================================"
echo "$(ts) RUN 2: Whole-BAM (picb_script_wholebam.R)"
echo "$(ts) ============================================================"

LOG_WHOLEBAM=${OUT_DIR}/E15_5_rep1_wholebam.log

conda run -n snakemake Rscript "$WHOLEBAM_R" \
    "$FASTA" "$BAM" "$OUT_WHOLEBAM" "$PIRNA_MIN" "$PIRNA_MAX" \
    2>&1 | tee "$LOG_WHOLEBAM"

if [ ! -f "$OUT_WHOLEBAM" ] || [ ! -s "$OUT_WHOLEBAM" ]; then
    echo "$(ts) ERROR: whole-BAM output not created or empty."
    exit 1
fi
echo "$(ts) Whole-BAM complete: $OUT_WHOLEBAM"

# ── RUN 3: Comparison figure ─────────────────────────────────────────────────
echo ""
echo "$(ts) ============================================================"
echo "$(ts) RUN 3: Generating comparison figure"
echo "$(ts) ============================================================"

conda run -n snakemake python3 "$FIGURE_PY" \
    "$OUT_CHRBYCHR" "$OUT_WHOLEBAM" "$OUT_DIR"

echo ""
echo "$(ts) ============================================================"
echo "$(ts) PICB comparison complete."
echo "$(ts) Results in: $OUT_DIR"
echo "$(ts) ============================================================"
