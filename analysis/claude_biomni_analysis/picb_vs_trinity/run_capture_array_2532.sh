#!/bin/bash
#SBATCH --job-name=pcap2532
#SBATCH --partition=1804
#SBATCH --array=1-48%16
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity/logs/cap2532_%a.log
# piRNA read-capture at the 25-32 window (purer; 24 nt is 1U-impure). total + in_picb + in_trin all recomputed at 25-32.
ST=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
DIR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity
BEDS=$DIR/beds16; OUTD=$DIR/cap16_results; mkdir -p $OUTD
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DIR/samplesheet16.tsv)
strain=$(echo "$LINE"|cut -f1); tpn=$(echo "$LINE"|cut -f2); tp=$(echo "$LINE"|cut -f3); bam=$(echo "$LINE"|cut -f4)
picb=$BEDS/$strain-$tp.picb.bed; trin=$BEDS/$strain-$tp.trin.bed
AW='length($10)>=25 && length($10)<=32'
echo "[$(date +%T)] $strain $tpn 25-32"
total=$($ST view -@ 4 "$bam" | awk "$AW" | wc -l)
in_picb=$($ST view -@ 4 -L "$picb" "$bam" | awk "$AW" | wc -l)
in_trin=$($ST view -@ 4 -L "$trin" "$bam" | awk "$AW" | wc -l)
echo "$strain,$tpn,$total,$in_picb,$in_trin" > "$OUTD/$strain-$tp.csv"
echo "[$(date +%T)] done $strain $tpn: total=$total in_picb=$in_picb in_trin=$in_trin"
