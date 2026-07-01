#!/bin/bash
#SBATCH --job-name=pcap16
#SBATCH --partition=1804
#SBATCH --array=1-48%16
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity/logs/cap_%a.log
# all-16-strain piRNA read-capture: one sample per task; counts 24-32 nt alignments total / in PICB / in Trinity exon blocks.
ST=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
DIR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity
BEDS=$DIR/beds16; OUTD=$DIR/cap16_results; mkdir -p $OUTD $DIR/logs
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DIR/samplesheet16.tsv)
strain=$(echo "$LINE" | cut -f1); tpn=$(echo "$LINE" | cut -f2); tp=$(echo "$LINE" | cut -f3); bam=$(echo "$LINE" | cut -f4)
picb=$BEDS/$strain-$tp.picb.bed; trin=$BEDS/$strain-$tp.trin.bed
echo "[$(date +%T)] task $SLURM_ARRAY_TASK_ID: $strain $tpn ($tp)"
total=$($ST view "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
in_picb=$($ST view -L "$picb" "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
in_trin=$($ST view -L "$trin" "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
echo "$strain,$tpn,$total,$in_picb,$in_trin" > "$OUTD/$strain-$tp.csv"
echo "[$(date +%T)] done $strain $tpn: total=$total in_picb=$in_picb in_trin=$in_trin"
