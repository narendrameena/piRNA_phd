#!/bin/bash
#SBATCH --job-name=pirpool16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/pool_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
echo "[start] $X $(hostname) $(date)"
"$PY" "$U/build_pools16.py" "$X"
echo "[done] $X $(date)"
