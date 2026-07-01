#!/bin/bash
# Full data-based edgeR-vs-DESeq2 benchmark, one task per timepoint (tp-wise, as requested).
#SBATCH --job-name=bench_da
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=23:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/bench_da/bench_%a.log
set -euo pipefail
TPS=(16.5dpc 12.5dpp 20.5dpp)
tp=${TPS[$SLURM_ARRAY_TASK_ID]}
R=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/Rscript
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
NPERM=${NPERM:-5}; MAXF=${MAXF:-300000}   # 300k random features = large representative sample; metrics are rates
echo "[$tp] benchmark start $(date) (nperm=$NPERM maxf=$MAXF)"
"$R" "$U/benchmark_da_methods.R" "$tp" "$NPERM" 1 "$MAXF"
echo "[$tp] done $(date)"
