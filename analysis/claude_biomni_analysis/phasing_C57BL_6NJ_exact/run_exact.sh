#!/bin/bash
set -u
SAM=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
OUT=$ROOT/analysis/claude_biomni_analysis/phasing_C57BL_6NJ_exact
PRE=$OUT/_prefilt
SCRIPT=$ROOT/workflow/scripts/R/phasing_analysis.R
PICB=$ROOT/results/STAR_srna_strain_wise_picb/C57BL_6NJ
NP=$ROOT/results/STAR_srna_strain_wise/C57BL_6NJ
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
# 16.5dpc first (picb temp BAMs at risk of SLURM deletion); ALL multimappers, 24-32 nt, follow()
for L in 16.5dpc.1 16.5dpc.3 12.5dpp.1 12.5dpp.2 12.5dpp.3 20.5dpp.1 20.5dpp.2 20.5dpp.3; do
  B=${BAMS[$L]}
  echo "[$(date +%T)] === $L ==="
  if [ ! -f "$B" ]; then echo "  MISSING BAM $B — skip"; continue; fi
  P=$PRE/$L.bam
  $SAM view -b -F 0x4 -e 'qlen>=24 && qlen<=32' -@ 4 "$B" -o "$P"
  $SAM index "$P"
  echo "  [$(date +%T)] prefiltered ALL-mappers $(du -h $P|cut -f1), alns $($SAM view -c $P)"
  $RS $SCRIPT "$P" "$OUT/C57BL_6NJ-$L" 24 32 50 0 all follow 2>&1 | grep -E "\+1 nt|usable|loaded|Error|error"
  rm -f "$P" "$P.bai"
  echo "  [$(date +%T)] done $L"
done
echo "[$(date +%T)] ALL DONE. Combined:"
head -1 $(ls $OUT/C57BL_6NJ-*.phasing_summary.csv 2>/dev/null|head -1) > $OUT/ALL_summary.csv
for f in $OUT/C57BL_6NJ-*.phasing_summary.csv; do tail -n +2 "$f" >> $OUT/ALL_summary.csv; done
cat $OUT/ALL_summary.csv
