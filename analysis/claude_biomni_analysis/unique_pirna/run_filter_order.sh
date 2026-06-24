#!/bin/bash
# Filter-order effect test (27/30 nt BEFORE vs AFTER DESeq2), one task per timepoint. Order B uses the FULL
# 25-32 filterByExpr set (no cap) for a faithful BH/normalization comparison -> heavy; generous resources.
#SBATCH --job-name=fo_test
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=47:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/deseq16_lenfilt/fo_%a.log
set -euo pipefail
TPS=(16.5dpc 12.5dpp 20.5dpp)
tp=${TPS[$SLURM_ARRAY_TASK_ID]}
R=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/Rscript
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
echo "[$tp] filter-order test start $(date)"
"$R" "$U/deseq_filter_order_test.R" "$tp" 0
echo "[$tp] done $(date)"
