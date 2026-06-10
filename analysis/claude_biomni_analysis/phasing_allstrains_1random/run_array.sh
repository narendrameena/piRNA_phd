#!/bin/bash
#SBATCH --job-name=phas16
#SBATCH --partition=1804
#SBATCH --array=1-144%24
#SBATCH --cpus-per-task=6
#SBATCH --mem=72G
#SBATCH --time=6:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_allstrains_1random/logs/task_%a.log
set -u
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
CUT=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/cutadapt
STAR=/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA/.snakemake/conda/3ed2c25f2a7c3a9c29f4adb59cbb09c0_/bin/STAR
SAM=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
RAW=$ROOT/resources/mice16_data/srna
OUT=$ROOT/analysis/claude_biomni_analysis/phasing_allstrains_1random
TMP=$OUT/_tmp
SCRIPT=$ROOT/workflow/scripts/R/phasing_analysis.R
MAN=$ROOT/analysis/claude_biomni_analysis/phasing_allstrains_manifest.tsv
SP="--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 30000000000 --outSAMmultNmax 1 --outMultimapperOrder Random --runRNGseed 777"

LN=$((SLURM_ARRAY_TASK_ID+1))
ROW=$(sed -n "${LN}p" $MAN)
STRAIN=$(echo "$ROW"|cut -f1); TP=$(echo "$ROW"|cut -f2); REP=$(echo "$ROW"|cut -f3)
IDX=$(echo "$ROW"|cut -f4); RAWF=$(echo "$ROW"|cut -f5)
L=$STRAIN-$TP.$REP
OUTPREF=$OUT/$L
if [ -f "$OUTPREF.phasing_summary.csv" ]; then echo "$L already done"; exit 0; fi
files=""; IFS=',' read -ra FA <<< "$RAWF"; for f in "${FA[@]}"; do files="$files $RAW/$f"; done
echo "[$(date +%T)] $L cutadapt ($RAWF)"
zcat $files | $CUT --minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -o $TMP/$L.trim.fastq - > $TMP/$L.cut.log 2>&1
echo "[$(date +%T)] $L STAR (1-random) on $IDX"
$STAR --runThreadN 6 --genomeDir $ROOT/$IDX --readFilesIn $TMP/$L.trim.fastq $SP --outFileNamePrefix $TMP/$L. > $TMP/$L.star.log 2>&1
B=$TMP/$L.Aligned.sortedByCoord.out.bam
$SAM index $B
echo "[$(date +%T)] $L phasing"
$RS $SCRIPT $B $OUTPREF 24 32 50 0 all follow > $TMP/$L.phas.log 2>&1
grep -E "\+1 nt|Error" $TMP/$L.phas.log
rm -f $TMP/$L.trim.fastq $B ${B}.bai; rm -rf $TMP/$L._STARtmp $TMP/$L.SJ.out.tab $TMP/$L.Log.progress.out $TMP/$L.Log.out
echo "[$(date +%T)] $L DONE"
