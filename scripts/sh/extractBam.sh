#!/bin/bash

# Check for required arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input.bam input.aligned_bam output.bam"
    exit 1
fi

# Assign arguments to variables
input_bam=$1
input_aligned_bam=$2
output_bam=$3

# Create a temporary directory
TMP_DIR=$(mktemp -d -t extract_reads_XXXXXX)

# Define temporary file path for unique read sequences
TMP_SEQ="$TMP_DIR/read_sequences.txt"

# Extract unique read sequences from input.bam
samtools view "$input_bam" | cut -f 10 | sort | uniq > "$TMP_SEQ"

# Filter reads and create the output BAM file with header
samtools view -H "$input_aligned_bam" | cat - <(awk 'NR==FNR {seq[$1]; next} $10 in seq' "$TMP_SEQ" <(samtools view "$input_aligned_bam")) | samtools view -b -o "$output_bam"

# Sort and index the output BAM
samtools sort "$output_bam" -o "$output_bam"
samtools index "$output_bam"

# Cleanup
rm -r "$TMP_DIR"
