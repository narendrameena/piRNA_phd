#!/bin/bash
# Scale-16 pangenome phase STEP 1: exact cross-strain expression classifier (loads the 16 pools as a
# per-sequence bitmask). Chain after edgeR DA (3302624) + pools (3302627). Big-mem (pool union ~30-50 GB).
#SBATCH --job-name=pirclassify16
#SBATCH --partition=2004,TEST
#SBATCH --cpus-per-task=4
#SBATCH --mem=240G
#SBATCH --time=8:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/classify16.log
set -euo pipefail
PY=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
echo "[start] $(hostname) $(date)"
"$PY" "$U/classify_unique16.py"
echo "[done] $(date)"
