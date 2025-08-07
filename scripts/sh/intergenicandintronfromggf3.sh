#!/bin/bash

# Check if input GFF3 is provided
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input.gff3>"
    exit 1
fi

input_gff3=$1

# Extract chromosome sizes from the GFF3 file
awk '{if($1 !~ /^#/ && $5+0>0) arr[$1]=($5>arr[$1]?$5:arr[$1])} END{for (i in arr) print i "\t" arr[i]}' $input_gff3 > chromSizes.txt

# Convert chromSizes.txt to BED format
awk 'OFS="\t" {print $1, "0", $2}' chromSizes.txt | sort -k1,1 -k2,2n > chromSizes.bed

# Sort the GFF3 file
cat $input_gff3 | awk '$1 ~ /^##sequence-region/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > sorted_input.gff3

# Get intergenic regions
bedtools complement -i sorted_input.gff3 -g chromSizes.bed > intergenic_sorted.bed

# Extract exonic coordinates in BED format
awk 'BEGIN{OFS="\t"} $1 !~ /^#/ {if ($3 == "exon") print $1, $4-1, $5}' sorted_input.gff3 > exon_sorted.bed

# Get introns
bedtools complement -i <(cat exon_sorted.bed intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chromSizes.bed > intron_sorted.bed
#check for overlap and take the union or max of exons(use bedtools) eg: use gene or exon  

# Convert intergenic_sorted.bed to intergenic.gff3
awk 'BEGIN{OFS="\t"} {print $1, "BEDtools", "intergenic", $2+1, $3, ".", ".", ".", "ID=intergenic_" NR}' intergenic_sorted.bed > intergenic.gff3

# Convert intron_sorted.bed to intron.gff3
awk 'BEGIN{OFS="\t"} {print $1, "BEDtools", "intron", $2+1, $3, ".", ".", ".", "ID=intron_" NR}' intron_sorted.bed > intron.gff3

echo "Processing complete. Output files: intergenic.gff3, intron.gff3"

