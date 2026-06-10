#!/bin/bash
#SBATCH --job-name=phasing_C57
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_C57BL_6NJ_exact/sbatch_%j.log
bash /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_C57BL_6NJ_exact/run_exact.sh
