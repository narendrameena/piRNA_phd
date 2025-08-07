#!/bin/bash
#SBATCH --job-name=bamCoverageJob
#SBATCH --output=bamCoverageJob_%A_%a.out
#SBATCH --error=bamCoverageJob_%A_%a.err
#SBATCH --array=1-3%3  # Adjust the array range according to your needs

# Load necessary modules if required
# module load <module_name>

# Change directory to the location of your bash script
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/sh

# Run the bash script for each array task
bash run_bam_covrage.sh
