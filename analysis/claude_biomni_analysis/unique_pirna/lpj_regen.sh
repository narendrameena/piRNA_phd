#!/bin/bash
# Regenerate the CORRUPT LP_J-12.5dpp.2 sRNA BAM (truncated/Invalid BGZF) with the EXACT validated sRNA params
# (config_para.json params.srna): cutadapt SE + STAR to the LP_J PanSN index. Verifies before overwriting.
#SBATCH --job-name=lpj_regen
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=8
#SBATCH --mem=72G
#SBATCH --time=5:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/lpj_regen.log
set -e
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
RAW=$ROOT/resources/mice16_data/srna/LP_J-12.5dpp.2.1s.fastq.gz
IDX=$ROOT/results/indexs/LP_J
OUT=$ROOT/results/STAR_srna_strain_wise/LP_J/LP_J-12.5dpp.2
TMP=$ROOT/analysis/claude_biomni_analysis/unique_pirna/.lpj_regen_tmp
CUTADAPT=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/cutadapt
STAR=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/STAR
SAMTOOLS=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
mkdir -p $TMP
echo "[cutadapt] $(date)"
$CUTADAPT -j 8 --minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -o $TMP/trimmed.fastq.gz "$RAW" > $TMP/cutadapt.log 2>&1
echo "[STAR] $(date)"
$STAR --runThreadN 8 --genomeDir $IDX --readFilesIn $TMP/trimmed.fastq.gz --readFilesCommand zcat --outFileNamePrefix $TMP/ \
  --outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 \
  --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 \
  --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 --limitBAMsortRAM 571734623688
echo "[verify] $(date)"
if $SAMTOOLS quickcheck $TMP/Aligned.sortedByCoord.out.bam; then
  cp $TMP/Aligned.sortedByCoord.out.bam $OUT/Aligned.sortedByCoord.out.bam
  $SAMTOOLS index -c $OUT/Aligned.sortedByCoord.out.bam
  echo "REGEN OK: $($SAMTOOLS view -c $OUT/Aligned.sortedByCoord.out.bam) records; $(date)"
  rm -rf $TMP
else
  echo "REGEN FAILED - corrupt STAR output, original NOT overwritten"
fi
