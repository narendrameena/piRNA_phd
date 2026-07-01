#!/bin/bash
#SBATCH --job-name=pirmatrix
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger/matrix_%a.log
set -euo pipefail
TPS=(16.5dpc 12.5dpp 20.5dpp)
TP=${TPS[$SLURM_ARRAY_TASK_ID]}
PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
"$PY" "$U/build_count_matrix.py" "$TP"
echo "[done] $TP $(date)"
