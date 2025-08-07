#!/bin/bash

# Define the arrays
timepoints=("16.5dpc" "12.5dpp" "20.5dpp")
replicates=(1 2 3)
strains=("129S1_SvImJ" "AKR_J" "A_J" "BALB_cJ" "C3H_HeJ" "CAST_EiJ" "CBA_J" "DBA_2J" "FVB_NJ" "LP_J" "NOD_ShiLtJ" "NZO_HlLtJ" "PWK_PhJ" "SPRET_EiJ" "WSB_EiJ" "C57BL_6NJ")

# Base directory
base_dir="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results"

# Iterate over each combination
for strain in "${strains[@]}"
do
    for timepoint in "${timepoints[@]}"
    do
        for replicate in "${replicates[@]}"
        do
            # Construct the directories and file paths
            result_dir="${base_dir}/TPMCount_rna/strains/${strain}-${timepoint}.${replicate}"
            bam_file="${base_dir}/STAR_rna_strain_wise/${strain}-${timepoint}.${replicate}/Aligned.sortedByCoord.out.bam"
            gtf_file="${base_dir}/annotationGTF/${strain}_v3.2.gtf"

            # Create directory and change to it
            mkdir -p "$result_dir"
            cd "$result_dir"

            # Execute the TPMCalculator command
            TPMCalculator -g "$gtf_file" -k gene_id -e -p -b "$bam_file"
        done
    done
done
