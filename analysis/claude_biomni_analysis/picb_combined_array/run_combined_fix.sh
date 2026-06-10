#!/bin/bash
#SBATCH --job-name=picbcombfix
#SBATCH --partition=1804,2204,2004
#SBATCH --array=1-29%16
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_combined_array/logs_fix/task_%a.log
set -u
# RERUN of the 29 incomplete/failed combined-PICB samples, with the race-condition FIX:
# each array task gets its OWN work dir ($TMPBASE/$L) so the chunked script's
# picb_chrom_temp (= dirname(BAM)/picb_chrom_temp) is unique per task and tasks can
# no longer rm -rf each other's chromosome BAMs. All cutadapt/STAR/PICB params unchanged.
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
TMPBASE=$R/analysis/claude_biomni_analysis/picb_combined_array/_tmp_fix
MAN=$R/analysis/claude_biomni_analysis/picb_combined_array/manifest_rerun.tsv
SP="--outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --limitBAMsortRAM 150000000000"
LN=$((SLURM_ARRAY_TASK_ID+1))
ROW=$(sed -n "${LN}p" $MAN); S=$(echo "$ROW"|cut -f1); TP=$(echo "$ROW"|cut -f2); RAWF=$(echo "$ROW"|cut -f3)
L=$S-$TP.combined; XLSX=$OUT/$S/$S-$TP.combined.xlsx
# NOTE: no "already done" skip -> this rerun overwrites any incomplete/corrupt xlsx.
mkdir -p $OUT/$S
WD=$TMPBASE/$L                      # <-- per-task isolated work dir (the fix)
rm -rf "$WD"; mkdir -p "$WD"
files=""; IFS=',' read -ra FA <<< "$RAWF"; for f in "${FA[@]}"; do files="$files $RAW/$f"; done
echo "[$(date +%T)] $L cutadapt (reps combined) -> $WD"
zcat $files | $CUT --minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -o $WD/$L.trim.fastq - > $WD/$L.cut.log 2>&1
echo "[$(date +%T)] $L STAR (PICB idx, all multimappers)"
$STAR --runThreadN 8 --genomeDir $PICBIDX/$S --readFilesIn $WD/$L.trim.fastq $SP --outFileNamePrefix $WD/$L. > $WD/$L.star.log 2>&1
B=$WD/$L.Aligned.sortedByCoord.out.bam; $SAM index $B
echo "[$(date +%T)] $L PICB (fixed env script, isolated temp)"
bash $CHUNK $FASTAD/${S}_chromosomes_MT.fasta $B $B.bai $XLSX $WD/$L.picb.log 8 $RS 999 5000 18 50 $L > $WD/$L.picbrun.log 2>&1
echo "[$(date +%T)] $L done xlsx=$([ -f $XLSX ] && echo OK || echo FAIL)"
# keep the small logs, drop the big intermediates; whole WD is self-contained
cp -f $WD/$L.picb.log $WD/$L.star.log $WD/$L.cut.log $TMPBASE/ 2>/dev/null
rm -rf "$WD"
