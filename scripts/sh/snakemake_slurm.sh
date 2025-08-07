#!/bin/bash

#SBATCH --job-name=Snakemake_Job
#SBATCH --output=Snakemake_Job.out
#SBATCH --error=Snakemake_Job.err
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --qos=group1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

# load modules
module load snakemake

# activate the virtual environment (if you are using one)
source activate myenv

# run Snakemake with the desired number of cores
snakemake --cores $SLURM_CPUS_PER_TASK --snakefile my_snakefile --jobs 10 --latency-wait 30

