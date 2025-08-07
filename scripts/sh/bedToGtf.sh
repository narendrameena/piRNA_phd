#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

# Assign input and output filenames to variables
input_file=$1
output_file=$2

# Execute the awk command with provided arguments
awk -v min_difference=500 'BEGIN { FS="\t"; OFS="\t"; id=1 } {
    if ($2 == 0) {
        $2 = 1
        $3 += 1
    }
    if (($3 - $2) > min_difference) {
        print $1, ".", "piRNA", $2, $3, ".", $4, ".", "ID \"piRNA" id "\";"
        id++
    }
}' "$input_file" > "$output_file"

echo "Conversion complete. Output written to $output_file."
