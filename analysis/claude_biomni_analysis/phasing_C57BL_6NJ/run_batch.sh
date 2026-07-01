#!/bin/bash
set -u
SAM=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
OUT=$ROOT/analysis/claude_biomni_analysis/phasing_C57BL_6NJ
PRE=$OUT/_prefilt
SCRIPT=$ROOT/workflow/scripts/R/phasing_analysis.R
PICB=$ROOT/results/STAR_srna_strain_wise_picb/C57BL_6NJ
NP=$ROOT/results/STAR_srna_strain_wise/C57BL_6NJ

# label : bam   (16.5dpc from picb; 12.5dpp & 20.5dpp from non-picb; identical STAR params)
declare -A BAMS=(
 [16.5dpc.1]=$PICB/C57BL_6NJ-16.5dpc.1/Aligned.sortedByCoord.out.bam
 [16.5dpc.3]=$PICB/C57BL_6NJ-16.5dpc.3/Aligned.sortedByCoord.out.bam
 [12.5dpp.1]=$NP/C57BL_6NJ-12.5dpp.1/Aligned.sortedByCoord.out.bam
 [12.5dpp.2]=$NP/C57BL_6NJ-12.5dpp.2/Aligned.sortedByCoord.out.bam
 [12.5dpp.3]=$NP/C57BL_6NJ-12.5dpp.3/Aligned.sortedByCoord.out.bam
 [20.5dpp.1]=$NP/C57BL_6NJ-20.5dpp.1/Aligned.sortedByCoord.out.bam
 [20.5dpp.2]=$NP/C57BL_6NJ-20.5dpp.2/Aligned.sortedByCoord.out.bam
 [20.5dpp.3]=$NP/C57BL_6NJ-20.5dpp.3/Aligned.sortedByCoord.out.bam
)
ORDER="16.5dpc.1 16.5dpc.3 12.5dpp.1 12.5dpp.2 12.5dpp.3 20.5dpp.1 20.5dpp.2 20.5dpp.3"
for L in $ORDER; do
  B=${BAMS[$L]}
  echo "[$(date +%T)] === $L ==="
  if [ ! -f "$B" ]; then echo "  MISSING BAM $B — skip"; continue; fi
  P=$PRE/$L.bam
  $SAM view -h -F 0x104 "$B" | awk '/^@/{print;next}{l=length($10);if(l>=25&&l<=33)print}' | $SAM view -b -o "$P" -
  $SAM index "$P"
  echo "  [$(date +%T)] prefiltered $(du -h $P|cut -f1), reads $($SAM view -c $P)"
  $RS $SCRIPT "$P" "$OUT/C57BL_6NJ-$L" 25 33 50 0 all follow 2>&1 | grep -E "\+1 nt|usable|ERROR|error" 
  rm -f "$P" "$P.bai"
  echo "  [$(date +%T)] done $L"
done
echo "[$(date +%T)] ALL DONE. Combined summary:"
head -1 $(ls $OUT/C57BL_6NJ-*.phasing_summary.csv 2>/dev/null|head -1) > $OUT/ALL_summary.csv
for f in $OUT/C57BL_6NJ-*.phasing_summary.csv; do tail -n +2 "$f" >> $OUT/ALL_summary.csv; done
cat $OUT/ALL_summary.csv
