#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0  <regions.bed> <genes.gff3> <output_folder>"
    exit 1
fi

GFF3_FILE="$2"
BED_FILE="$1"
OUTPUT_FOLDER="$3"  # Assigning the third argument as the output folder

# Create the output folder if it doesn't exist
if [ ! -d "$OUTPUT_FOLDER" ]; then
    mkdir -p "$OUTPUT_FOLDER"
fi

# Check if the input files are non-empty
if [ ! -s "$GFF3_FILE" ] || [ ! -s "$BED_FILE" ]; then
    echo "Input files are empty or do not exist."
    exit 1
fi

# Extract chromosomes from GFF3 and BED files
awk '$3 == "gene" {print $1}' "$GFF3_FILE" | sort -u > "${OUTPUT_FOLDER}/gff3_chromosomes.txt"
if [ ! -s "${OUTPUT_FOLDER}/gff3_chromosomes.txt" ]; then
    echo "No gene chromosomes found in GFF3 file."
    exit 1
fi

cut -f1 "$BED_FILE" | sort -u > "${OUTPUT_FOLDER}/bed_chromosomes.txt"
if [ ! -s "${OUTPUT_FOLDER}/bed_chromosomes.txt" ]; then
    echo "No chromosomes found in BED file."
    exit 1
fi

# Find common chromosomes, exit if none found
comm -12 "${OUTPUT_FOLDER}/gff3_chromosomes.txt" "${OUTPUT_FOLDER}/bed_chromosomes.txt" > "${OUTPUT_FOLDER}/common_chromosomes.txt"
if [ ! -s "${OUTPUT_FOLDER}/common_chromosomes.txt" ]; then
    echo "No common chromosomes found."
    exit 1
fi

# Filter GFF3 and BED files for common chromosomes only
awk 'NR==FNR {chr[$1]; next} $1 in chr' "${OUTPUT_FOLDER}/common_chromosomes.txt" "$GFF3_FILE" > "${OUTPUT_FOLDER}/filtered.gff3"
awk 'NR==FNR {chr[$1]; next} $1 in chr' "${OUTPUT_FOLDER}/common_chromosomes.txt" "$BED_FILE" > "${OUTPUT_FOLDER}/filtered.bed"

# Convert GFF3 to BED and store it in the output location
awk '$3 == "gene" {print $1"\t"$4-1"\t"$5"\t"$7"\t"$9}' "${OUTPUT_FOLDER}/filtered.gff3" | tr -d '' | tr -s ' ' | cut -d ' ' -f1,2,3,4,6 > "${OUTPUT_FOLDER}/genes_unsorted.bed"

# Create a fake genome file assuming the maximum chromosome length is 1e9
cut -f1 "${OUTPUT_FOLDER}/genes_unsorted.bed" | sort -u | awk 'BEGIN{OFS="\t"}{print $1, "1000000000"}' > "${OUTPUT_FOLDER}/genome.file"

# Sort the genes BED file
bedtools sort -i "${OUTPUT_FOLDER}/genes_unsorted.bed" -g "${OUTPUT_FOLDER}/genome.file" > "${OUTPUT_FOLDER}/genes.bed"

# Sort the original BED file
bedtools sort -i "${OUTPUT_FOLDER}/filtered.bed" -g "${OUTPUT_FOLDER}/genome.file" > "${OUTPUT_FOLDER}/sorted_regions.bed"

# Find overlaps with bedtools intersect, considering same strand
bedtools intersect -a "${OUTPUT_FOLDER}/sorted_regions.bed" -b "${OUTPUT_FOLDER}/genes.bed" -wa -wb  > "${OUTPUT_FOLDER}/overlaps.bed"

# Rename the BED entries with overlaps, preserving original name
awk -F$'\t' 'BEGIN{OFS=FS}{
  if(NF > 12){
    $4=$4":overlap-"$17
  }
  print $0
}' "${OUTPUT_FOLDER}/overlaps.bed" > "${OUTPUT_FOLDER}/final_output.bed"

# Uncomment the line below if you want to remove intermediate files
# rm "${OUTPUT_FOLDER}/gff3_chromosomes.txt" "${OUTPUT_FOLDER}/bed_chromosomes.txt" "${OUTPUT_FOLDER}/common_chromosomes.txt" "${OUTPUT_FOLDER}/filtered.gff3" "${OUTPUT_FOLDER}/filtered.bed" "${OUTPUT_FOLDER}/genes_unsorted.bed" "${OUTPUT_FOLDER}/genome.file" "${OUTPUT_FOLDER}/genes.bed" "${OUTPUT_FOLDER}/sorted_regions.bed" "${OUTPUT_FOLDER}/overlaps.bed"

# Final echo statement
echo "Renamed BED file with overlaps created in $OUTPUT_FOLDER"
