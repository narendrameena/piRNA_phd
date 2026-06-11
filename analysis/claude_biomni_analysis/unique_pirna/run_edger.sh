#!/bin/bash
#SBATCH --job-name=piredger
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=3:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger/edger_%a.log
set -euo pipefail
TPS=(16.5dpc 12.5dpp 20.5dpp)
TP=${TPS[$SLURM_ARRAY_TASK_ID]}
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
"$RS" "$U/edger_da.R" "$TP"
echo "[done] $TP $(date)"
