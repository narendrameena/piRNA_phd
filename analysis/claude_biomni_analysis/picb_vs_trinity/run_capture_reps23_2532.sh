#!/bin/bash
#SBATCH --job-name=pcaprep
#SBATCH --partition=1804
#SBATCH --array=1-96%16
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity/logs/caprep_%a.log
# reps 2,3 piRNA read-capture at 25-32 (per-rep) for replicate dots/error bars.
ST=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
DIR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity
BEDS=$DIR/beds16; OUTD=$DIR/cap16_reps
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DIR/samplesheet16_reps23.tsv)
strain=$(echo "$LINE"|cut -f1); tpn=$(echo "$LINE"|cut -f2); tp=$(echo "$LINE"|cut -f3); rep=$(echo "$LINE"|cut -f4); bam=$(echo "$LINE"|cut -f5)
picb=$BEDS/$strain-$tp.picb.bed; trin=$BEDS/$strain-$tp.trin.bed
AW='length($10)>=25 && length($10)<=32'
out=$OUTD/$strain-$tp.$rep.csv
[ -f "$out" ] && { echo "$strain-$tp.$rep done"; exit 0; }
echo "[$(date +%T)] $strain $tpn rep$rep 25-32"
total=$($ST view -@ 4 "$bam" | awk "$AW" | wc -l)
in_picb=$($ST view -@ 4 -L "$picb" "$bam" | awk "$AW" | wc -l)
in_trin=$($ST view -@ 4 -L "$trin" "$bam" | awk "$AW" | wc -l)
echo "$strain,$tpn,$rep,$total,$in_picb,$in_trin" > "$out"
echo "[$(date +%T)] done $strain $tpn rep$rep"
