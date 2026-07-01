#!/bin/bash
#SBATCH --job-name=te_sa
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --array=0-15
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/te_sa_%a.log
/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/build_te_sense_antisense.py
