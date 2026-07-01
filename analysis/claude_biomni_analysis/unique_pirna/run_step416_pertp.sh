#!/bin/bash
# WITHIN-TIMEPOINT D4 (SNP-variant refinement) for all 16 strains. Reuses the existing step4_16 cand_to_Y
# BAMs (no STAR) + the new per-tp pools (pools_pertp/) -> lightweight (pysam read + pool streaming).
#SBATCH --job-name=snp_withintp
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4_16/snp_withintp_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
PY=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
echo "[$X] start $(date)"
"$PY" "$U/classify_step416_pertp.py" "$X"
echo "[done] $X $(date)"
