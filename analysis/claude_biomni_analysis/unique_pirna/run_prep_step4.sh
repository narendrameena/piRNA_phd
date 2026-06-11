#!/bin/bash
#SBATCH --job-name=prep_s4
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4/prep_%a.log
set -euo pipefail
STRAINS=(C57BL_6NJ CAST_EiJ SPRET_EiJ)
S=${STRAINS[$SLURM_ARRAY_TASK_ID]}
PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
"$PY" "$U/prep_step4.py" "$S"
echo "[done] $S $(date)"
