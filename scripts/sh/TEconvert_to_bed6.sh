#!/bin/bash
input=$1
awk -v filename=$(basename "$input" .shared.tab) 'BEGIN {OFS="\t"} {print $1, $2, $3, filename, 0, $4}' "$input" | grep -v  'chrm' > "${input%}.bed"
