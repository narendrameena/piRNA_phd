#!/bin/bash

# Define arrays for tps, reps, and strains
tps=("16.5dpc" "12.5dpp" "20.5dpp")
reps=(1 2 3)
strains=(
    "129S1_SvImJ" "AKR_J" "A_J" "BALB_cJ" "C3H_HeJ" "CAST_EiJ" "CBA_J" "DBA_2J" "FVB_NJ" "LP_J"
    "NOD_ShiLtJ" "NZO_HlLtJ" "PWK_PhJ" "SPRET_EiJ" "WSB_EiJ" "C57BL_6NJ" "C57BL_6"
)

# Loop through combinations of tps, reps, and strains
for tp in "${tps[@]}"; do
    for rep in "${reps[@]}"; do
        for strain in "${strains[@]}"; do
            # Run the command for each combination
            bash /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/sh/piRNA_filter_precursors.sh \
                "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/trinity/filter_precursors/$strain/$strain-$tp.$rep.500.csv" \
                "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/minimap2/strains/$strain/$strain-$tp.500.bed12.bed" \
                "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/filter_precursors_bed/$strain/$strain-$tp.$rep.bed"
        done
    done
done
