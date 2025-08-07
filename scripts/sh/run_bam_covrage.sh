#!/bin/bash

# Define arrays for tp, rep, and sid
tp=("16.5dpc" "12.5dpp" "20.5dpp")
rep=(1 2 3)
sid=(
    "129S1_SvImJ"
    "AKR_J"
    "A_J"
    "BALB_cJ"
    "C3H_HeJ"
    "CAST_EiJ"
    "CBA_J"
    "DBA_2J"
    "FVB_NJ"
    "LP_J"
    "NOD_ShiLtJ"
    "NZO_HlLtJ"
    "PWK_PhJ"
    "SPRET_EiJ"
    "WSB_EiJ"
    "C57BL_6NJ"
)

# Function to submit each bamCoverage command as a separate Slurm job
submit_job() {
    local cmd=$1
    local log_file=$2
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=bamCoverageJob
#SBATCH --output=${log_file}_%j.out
#SBATCH --error=${log_file}_%j.err
#SBATCH --partition=NXFL  # Adjust with your desired Slurm partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=400G  # Adjust memory as needed

# Load necessary modules if required
# module load <module_name>

${cmd}
EOT
}

# Iterate through all combinations and submit bamCoverage commands as separate Slurm jobs
for t in "${tp[@]}"; do
    for r in "${rep[@]}"; do
        for s in "${sid[@]}"; do
            bam_cmd_plus="bamCoverage --normalizeUsing CPM -b \"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/$s/$s-$t.$r/Aligned.sortedByCoord.out.bam\" -o \"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/bamCoverageSrna/pacBio/$s/$s-$t.${r}_plusStrand.bw\" --filterRNAstrand forward"
            bam_cmd_minus="bamCoverage --normalizeUsing CPM -b \"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/STAR_srna_strain_wise/$s/$s-$t.$r/Aligned.sortedByCoord.out.bam\" -o \"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/bamCoverageSrna/pacBio/$s/$s-$t.${r}_minusStrand.bw\" --filterRNAstrand reverse"

            submit_job "$bam_cmd_plus" "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/log/bamCoverageSrnaBW/pacBio/$s/$s-$t.${r}_forward"
            submit_job "$bam_cmd_minus" "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/log/bamCoverageSrnaBW/pacBio/$s/$s-$t.${r}_reverse"
        done
    done
done
