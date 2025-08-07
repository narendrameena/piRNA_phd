#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <first_file> <second_file> <output_file>"
    exit 1
fi

# Store input and output filenames
first_file=$1
second_file=$2
output_file=$3

# Extract transcript IDs with RPM >= 100 from the first file
rpm_filtered_ids=$(awk -F'\t' '$8 >= 100 { print $1 }' "$first_file")

# Filter the second file based on an exact match using grep
grep -Fwf <(echo "$rpm_filtered_ids") "$second_file" > "$output_file"
