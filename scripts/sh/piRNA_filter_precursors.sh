#!/bin/bash

# Check if three arguments are given (for the two input files and the output file)
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 first_file.csv second_file.txt output_file.txt"
    exit 1
fi

first_file=$1
second_file=$2
output_file=$3

# Extract unique values from the first file and use grep to filter the second file
awk -F, '($3 > 100 && $4 > 100) { print $1 }' "$first_file" | sort -u | grep -Ff - "$second_file" > "$output_file"

echo "Processing complete. Output saved to $output_file"
