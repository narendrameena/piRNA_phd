#!/bin/bash

# Define the list of strains, timepoints, and replicates
strains=(
    'C57BL_6NJ' 'BALB_cJ' 'A_J' 'FVB_NJ' 'C3H_HeJ' '129S1_SvImJ' 'DBA_2J' 'AKR_J' 'CBA_J' 'NOD_ShiLtJ' 'WSB_EiJ' 'CAST_EiJ' 'PWK_PhJ' 'SPRET_EiJ'
)
timepoints=('16.5dpc' '12.5dpp' '20.5dpp')
replicates=('1' '2' '3')

# Loop over each strain
for strain in "${strains[@]}"; do
    # Loop over each timepoint
    for timepoint in "${timepoints[@]}"; do
        # Loop over each replicate
        for replicate in "${replicates[@]}"; do
            # Define the output folder path
            output_folder="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/filter_precursors_bed/${strain}"
            
            # Create the output folder if it doesn't exist
            mkdir -p "${output_folder}"

            # Define the output file path
            output_path="${output_folder}/${strain}-${timepoint}.${replicate}.500.bed"

            # Run the script with the current combination of strain, timepoint, and replicate
            bash "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/sh/piRNA_filter_precursors.sh" \
                "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/trinity/filter_precursors/${strain}/${strain}-${timepoint}.${replicate}.500.csv" \
                "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/minimap2/strains/${strain}/${strain}-${timepoint}.500.bed12.bed" \
                "${output_path}"
        done
    done
done
