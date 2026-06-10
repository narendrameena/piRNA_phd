#!/bin/bash
#SBATCH --job-name=phas1rand
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_C57BL_6NJ_1random/sbatch_%j.log
set -u
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
CUT=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/cutadapt
STAR=/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA/.snakemake/conda/3ed2c25f2a7c3a9c29f4adb59cbb09c0_/bin/STAR
SAM=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
IDX=$ROOT/results/indexs/C57BL_6NJ
RAW=$ROOT/resources/mice16_data/srna
OUT=$ROOT/analysis/claude_biomni_analysis/phasing_C57BL_6NJ_1random
TMP=$OUT/_tmp
SCRIPT=$ROOT/workflow/scripts/R/phasing_analysis.R
# STAR: identical to original BAMs + 1-random-coordinate (multNmax 1 + Random order, fixed seed)
SP="--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --limitBAMsortRAM 150000000000 --outSAMmultNmax 1 --outMultimapperOrder Random --runRNGseed 777"

# raw-file suffixes per sample (prefix C57BL_6NJ-<L>. , suffix .fastq.gz)
declare -A SUF=(
 [16.5dpc.1]="1g 1s 2s" [16.5dpc.2]="3s 1g 1s 2s" [16.5dpc.3]="3s 1g 1s 4s 2s"
 [12.5dpp.1]="1g 1s 2s" [12.5dpp.2]="1g 1s 2s" [12.5dpp.3]="1g 1s 2s"
 [20.5dpp.1]="1g 1s 2s" [20.5dpp.2]="3s 1g 1s 2s" [20.5dpp.3]="3s 1g 1s 2s"
)
for L in 16.5dpc.1 16.5dpc.2 16.5dpc.3 12.5dpp.1 12.5dpp.2 12.5dpp.3 20.5dpp.1 20.5dpp.2 20.5dpp.3; do
  echo "[$(date +%T)] === $L : cutadapt ==="
  files=""; for s in ${SUF[$L]}; do files="$files $RAW/C57BL_6NJ-$L.$s.fastq.gz"; done
  zcat $files | $CUT --minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -o $TMP/$L.trim.fastq - > $TMP/$L.cutadapt.log 2>&1
  echo "[$(date +%T)] $L : STAR (1-random) ; trimmed reads=$(( $(wc -l < $TMP/$L.trim.fastq)/4 ))"
  $STAR --runThreadN 16 --genomeDir $IDX --readFilesIn $TMP/$L.trim.fastq $SP --outFileNamePrefix $TMP/$L. > $TMP/$L.star.log 2>&1
  B=$TMP/$L.Aligned.sortedByCoord.out.bam
  $SAM index $B
  echo "[$(date +%T)] $L : aligned=$($SAM view -c $B) ; phasing"
  $RS $SCRIPT $B $OUT/C57BL_6NJ-$L 24 32 50 0 all follow 2>&1 | grep -E "\+1 nt|total|Error"
  rm -f $TMP/$L.trim.fastq $B ${B}.bai; rm -rf $TMP/$L._STARtmp $TMP/$L.Unmapped* $TMP/$L.SJ* $TMP/$L.Log.progress.out
  echo "[$(date +%T)] $L DONE"
done
echo "[$(date +%T)] ALL DONE"
head -1 $(ls $OUT/C57BL_6NJ-*.phasing_summary.csv 2>/dev/null|head -1) > $OUT/ALL_summary.csv
for f in $OUT/C57BL_6NJ-*.phasing_summary.csv; do tail -n +2 "$f" >> $OUT/ALL_summary.csv; done
cat $OUT/ALL_summary.csv
