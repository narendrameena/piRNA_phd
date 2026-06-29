#!/bin/bash
#SBATCH --job-name=pcapC3H
#SBATCH --partition=1804
#SBATCH --array=13
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity/logs/capC3H_%a.log
# recover total + in_picb for C3H_HeJ-16.5dpc (65 GB BAM; lost in the re-run chain). in_trin (100/100) reused from cap16_results.
ST=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
DIR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity
BEDS=$DIR/beds16; OUTD=$DIR/cap16_results
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DIR/samplesheet16.tsv)
strain=$(echo "$LINE"|cut -f1); tpn=$(echo "$LINE"|cut -f2); tp=$(echo "$LINE"|cut -f3); bam=$(echo "$LINE"|cut -f4)
in_trin=$(cut -d, -f5 "$OUTD/$strain-$tp.csv")   # existing 100/100 in_trin
picb=$BEDS/$strain-$tp.picb.bed
echo "[$(date +%T)] $strain $tpn: total+in_picb (reuse in_trin=$in_trin)"
total=$($ST view -@ 8 "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
in_picb=$($ST view -@ 8 -L "$picb" "$bam" | awk 'length($10)>=24 && length($10)<=32' | wc -l)
echo "$strain,$tpn,$total,$in_picb,$in_trin" > "$OUTD/$strain-$tp.csv"
echo "[$(date +%T)] done: total=$total in_picb=$in_picb in_trin=$in_trin"
