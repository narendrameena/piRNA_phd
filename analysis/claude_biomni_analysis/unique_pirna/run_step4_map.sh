#!/bin/bash
#SBATCH --job-name=step4map
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4/map_%a.log
set -euo pipefail
STRAINS=(C57BL_6NJ CAST_EiJ SPRET_EiJ)
X=${STRAINS[$SLURM_ARRAY_TASK_ID]}
STAR=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/STAR
SAMTOOLS=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/samtools
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
S4=$ROOT/analysis/claude_biomni_analysis/unique_pirna/step4
# sRNA STAR params (identical to config params.srna.STAR); ONLY --outFilterMismatchNmax varies:
# 0 for X's own candidates (origin loci), 3 for mapping candidates to OTHER genomes (data-driven
# Poisson SNP cutoff k*=3). --limitBAMsortRAM reduced to 30G (resource setting, not biological).
SRNA="--outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 \
 --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All \
 --scoreDelOpen -10000 --scoreInsOpen -10000 --outSAMtype BAM SortedByCoordinate \
 --outBAMsortingThreadN 1 --limitBAMsortRAM 30000000000 --runThreadN 8"
case $X in
  C57BL_6NJ) OTHERS="CAST_EiJ SPRET_EiJ";;
  CAST_EiJ)  OTHERS="C57BL_6NJ SPRET_EiJ";;
  SPRET_EiJ) OTHERS="C57BL_6NJ CAST_EiJ";;
esac
FA="$S4/$X.candidates.fasta"
echo "[$X] candidates=$(grep -c '^>' "$FA"); others=$OTHERS"

# (1) candidates -> X's OWN genome, exact (mm0): origin loci + confirms X-expressed mapping
$STAR --runMode alignReads --genomeDir "$ROOT/results/indexs/$X" --readFilesIn "$FA" \
  --outFilterMismatchNmax 0 $SRNA --outFileNamePrefix "$S4/$X.cand_self."
$SAMTOOLS index "$S4/$X.cand_self.Aligned.sortedByCoord.out.bam"

# (2) candidates -> each OTHER genome, <=3 mm: homologous locus + NM (SNP distance)
for Y in $OTHERS; do
  $STAR --runMode alignReads --genomeDir "$ROOT/results/indexs/$Y" --readFilesIn "$FA" \
    --outFilterMismatchNmax 3 $SRNA --outFileNamePrefix "$S4/$X.cand_to_$Y."
  $SAMTOOLS index "$S4/$X.cand_to_$Y.Aligned.sortedByCoord.out.bam"
done
echo "[done] $X $(date)"
