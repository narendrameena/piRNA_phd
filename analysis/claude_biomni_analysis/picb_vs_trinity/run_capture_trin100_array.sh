#!/bin/bash
#SBATCH --job-name=pcap100
#SBATCH --partition=1804
#SBATCH --array=1-48%16
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=6:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity/logs/cap100_%a.log
# LEAN re-run at thesis 100/100: only in_trin changes (more Trinity precursors); total + in_picb are identical
# (same BAM, same PICB bed) -> reuse them from cap16_results_200bak. Threaded view -L (-@4).
ST=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
DIR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity
BEDS=$DIR/beds16; OUTD=$DIR/cap16_results; OLD=$DIR/cap16_results_200bak; mkdir -p $OUTD
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DIR/samplesheet16.tsv)
strain=$(echo "$LINE"|cut -f1); tpn=$(echo "$LINE"|cut -f2); tp=$(echo "$LINE"|cut -f3); bam=$(echo "$LINE"|cut -f4)
oldcsv=$OLD/$strain-$tp.csv
if [ ! -f "$oldcsv" ]; then echo "no old csv for $strain-$tp (skipping in_picb reuse)"; total=NA; in_picb=NA
else total=$(cut -d, -f3 "$oldcsv"); in_picb=$(cut -d, -f4 "$oldcsv"); fi
trin=$BEDS/$strain-$tp.trin.bed
echo "[$(date +%T)] task $SLURM_ARRAY_TASK_ID: $strain $tpn ($tp) in_trin @100/100"
in_trin=$($ST view -@ 4 -L "$trin" "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
echo "$strain,$tpn,$total,$in_picb,$in_trin" > "$OUTD/$strain-$tp.csv"
echo "[$(date +%T)] done $strain $tpn: total=$total in_picb=$in_picb in_trin=$in_trin"
