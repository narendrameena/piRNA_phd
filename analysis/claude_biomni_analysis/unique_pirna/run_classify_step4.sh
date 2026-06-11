#!/bin/bash
#SBATCH --job-name=step4cls
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4/classify_%a.log
set -euo pipefail
STRAINS=(C57BL_6NJ CAST_EiJ SPRET_EiJ)
X=${STRAINS[$SLURM_ARRAY_TASK_ID]}
PY=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/python
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
"$PY" "$U/classify_step4.py" "$X"
echo "[done] $X $(date)"
