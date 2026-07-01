#!/bin/bash
#SBATCH --job-name=pcap16
#SBATCH --partition=1804
#SBATCH --array=13,19,34
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=8:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity/logs/cap_%a.log
# re-run of the 3 timed-out E16.5 tasks (56-65 GB BAMs) with threaded decompression (-@ 4) + 8 h limit.
ST=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
DIR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity
BEDS=$DIR/beds16; OUTD=$DIR/cap16_results; mkdir -p $OUTD
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DIR/samplesheet16.tsv)
strain=$(echo "$LINE" | cut -f1); tpn=$(echo "$LINE" | cut -f2); tp=$(echo "$LINE" | cut -f3); bam=$(echo "$LINE" | cut -f4)
picb=$BEDS/$strain-$tp.picb.bed; trin=$BEDS/$strain-$tp.trin.bed
echo "[$(date +%T)] task $SLURM_ARRAY_TASK_ID: $strain $tpn ($tp) -@4"
total=$($ST view -@ 4 "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
in_picb=$($ST view -@ 4 -L "$picb" "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
in_trin=$($ST view -@ 4 -L "$trin" "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
echo "$strain,$tpn,$total,$in_picb,$in_trin" > "$OUTD/$strain-$tp.csv"
echo "[$(date +%T)] done $strain $tpn: total=$total in_picb=$in_picb in_trin=$in_trin"
