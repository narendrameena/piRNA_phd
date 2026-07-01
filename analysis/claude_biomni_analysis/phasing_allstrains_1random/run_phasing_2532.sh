#!/bin/bash
#SBATCH --job-name=phas2532
#SBATCH --partition=1804
#SBATCH --array=1-144%24
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_allstrains_25_32/logs/task_%a.log
# phasing at 25-32 on the EXISTING sRNA STAR BAMs (multimap=primary -> 1 primary alignment/read) — NO re-alignment.
set -u
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
RS=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/Rscript
SCRIPT=$ROOT/workflow/scripts/R/phasing_analysis.R
MAN=$ROOT/analysis/claude_biomni_analysis/phasing_allstrains_manifest.tsv
OUT=$ROOT/analysis/claude_biomni_analysis/phasing_allstrains_25_32; mkdir -p $OUT/logs
LN=$((SLURM_ARRAY_TASK_ID+1)); ROW=$(sed -n "${LN}p" $MAN)
STRAIN=$(echo "$ROW"|cut -f1); TP=$(echo "$ROW"|cut -f2); REP=$(echo "$ROW"|cut -f3)
L=$STRAIN-$TP.$REP
BAM=$ROOT/results/STAR_srna_strain_wise/$STRAIN/$STRAIN-$TP.$REP/Aligned.sortedByCoord.out.bam
OUTPREF=$OUT/$L
[ -f "$OUTPREF.phasing_summary.csv" ] && { echo "$L done"; exit 0; }
[ -f "$BAM" ] || { echo "NO BAM $BAM"; exit 1; }
[ -f "$BAM.bai" ] || /mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools index "$BAM"
echo "[$(date +%T)] $L phasing 25-32 primary on existing BAM"
$RS $SCRIPT "$BAM" "$OUTPREF" 25 32 50 0 primary follow > $OUT/logs/$L.phas.log 2>&1
grep -E "\+1 nt|frac|Error" $OUT/logs/$L.phas.log | head -2
echo "[$(date +%T)] $L DONE"
