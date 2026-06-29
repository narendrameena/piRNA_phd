#!/bin/bash
export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"
mkdir -p log/cluster

snakemake --printshellcmds --use-conda --conda-frontend conda -s workflow/Strain_Snakefile_black6_E16_5 --latency-wait 100 --output-wait 100 --jobs 999 $1 --max-status-checks-per-second 500 --profile /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/config/
#snakemake --printshellcmds --use-conda -s workflow/Strain_Snakefile --latency-wait 500 --output-wait 500 --jobs 5 $1 -T 5 --profile  config/cluster_config_slurm.json --cluster 'sbatch --output=log/cluster/{rule}.{wildcards}.out --error=log/cluster/{rule}.{wildcards}.err --partition={cluster.queue} --cpus-per-task={cluster.slots} --mem={cluster.mem}'
# --nodelist={cluster.nodelist}'