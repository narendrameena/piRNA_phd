#!/bin/bash
#SBATCH --job-name=picbcomb
#SBATCH --partition=1804,2204,2004
#SBATCH --array=1-48%16
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=14:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/logs/task_%a.log
set -u
R=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
CUT=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/cutadapt
STAR=/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA/.snakemake/conda/3ed2c25f2a7c3a9c29f4adb59cbb09c0_/bin/STAR
SAM=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
RAW=$R/resources/mice16_data/srna
PICBIDX=$R/resources/PICB/index
FASTAD=$R/resources/PICB/refFasta
CHUNK=$R/workflow/scripts/run_picb_analysis_chunked.sh
RS=$R/workflow/scripts/R/picb_script_chunked.R
OUT=$R/results/picb_result_combined
TMP=$R/analysis/claude_biomni_analysis/picb_combined_array/_tmp
MAN=$R/analysis/claude_biomni_analysis/picb_combined_array/manifest.tsv
SP="--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --limitBAMsortRAM 150000000000"
LN=$((SLURM_ARRAY_TASK_ID+1))
ROW=$(sed -n "${LN}p" $MAN); S=$(echo "$ROW"|cut -f1); TP=$(echo "$ROW"|cut -f2); RAWF=$(echo "$ROW"|cut -f3)
L=$S-$TP.combined; XLSX=$OUT/$S/$S-$TP.combined.xlsx
if [ -f "$XLSX" ]; then echo "$L already done"; exit 0; fi
mkdir -p $OUT/$S $TMP
files=""; IFS=',' read -ra FA <<< "$RAWF"; for f in "${FA[@]}"; do files="$files $RAW/$f"; done
echo "[$(date +%T)] $L cutadapt (reps combined)"
zcat $files | $CUT --minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -o $TMP/$L.trim.fastq - > $TMP/$L.cut.log 2>&1
echo "[$(date +%T)] $L STAR (PICB idx, all multimappers)"
$STAR --runThreadN 8 --genomeDir $PICBIDX/$S --readFilesIn $TMP/$L.trim.fastq $SP --outFileNamePrefix $TMP/$L. > $TMP/$L.star.log 2>&1
B=$TMP/$L.Aligned.sortedByCoord.out.bam; $SAM index $B
echo "[$(date +%T)] $L PICB (fixed env script)"
bash $CHUNK $FASTAD/${S}_chromosomes_MT.fasta $B $B.bai $XLSX $TMP/$L.picb.log 8 $RS 999 5000 18 50 $L > $TMP/$L.picbrun.log 2>&1
echo "[$(date +%T)] $L done xlsx=$([ -f $XLSX ] && echo OK || echo FAIL)"
rm -f $TMP/$L.trim.fastq $B $B.bai; rm -rf $TMP/$L._STARtmp $TMP/$L.SJ.out.tab $TMP/$L.Log.progress.out $TMP/$L.Log.out
