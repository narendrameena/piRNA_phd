#!/bin/bash
#SBATCH --job-name=uniqpirna
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=2
#SBATCH --mem=220G
#SBATCH --time=4:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/caller.log
set -euo pipefail
PY=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/python
# presence = reproducibility (detected in >=2/3 reps), 24-32 nt window; no magic RPM threshold (user-locked 2026-06-10)
"$PY" /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/build_unique_pirna_pilot.py 2
echo "[done] $(date)"
