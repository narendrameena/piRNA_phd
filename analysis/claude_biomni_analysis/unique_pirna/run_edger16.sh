#!/bin/bash
#SBATCH --job-name=piredger16
#SBATCH --partition=2004,TEST
#SBATCH --cpus-per-task=8
#SBATCH --mem=480G
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16/edger16_%a.log
set -euo pipefail
TPS=(16.5dpc 12.5dpp 20.5dpp)
TP=${TPS[$SLURM_ARRAY_TASK_ID]}
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
echo "[start] $TP $(hostname) $(date)"
"$RS" "$U/edger16.R" "$TP"
echo "[done] $TP $(date)"
