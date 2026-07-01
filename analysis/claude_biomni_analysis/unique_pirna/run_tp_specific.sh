#!/bin/bash
#SBATCH --job-name=tpspec
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger_tp/tpspec_%a.log
set -euo pipefail
STR=(C57BL_6NJ CAST_EiJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
"$PY" "$U/build_count_matrix_perstrain.py" "$X"
"$RS" "$U/edger_da_timepoint.R" "$X"
echo "[done] $X $(date)"
