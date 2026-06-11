#!/bin/bash
#SBATCH --job-name=pirpca
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pca/pca_%a.log
set -euo pipefail
TPS=(16.5dpc 12.5dpp 20.5dpp); LABS=(E16.5 P12.5 P20.5)
i=$SLURM_ARRAY_TASK_ID
RS=/mnt/home3/miska/nm667/miniconda3/envs/cosmaxSpatial/bin/Rscript
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
"$RS" "$U/pca_unique.R" "${TPS[$i]}" "${LABS[$i]}"
echo "[done] ${LABS[$i]} $(date)"
