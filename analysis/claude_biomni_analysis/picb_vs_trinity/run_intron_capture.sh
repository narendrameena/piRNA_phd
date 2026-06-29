#!/bin/bash
# Read-budget extras for the capture figure: per (strain,tp,rep) compute (1) TOTAL small-RNA alignments (any length)
# and (2) 25-32 nt piRNA reads inside Trinity FULL-SPAN precursors (intron-INCLUDED). Merged with existing
# cap16_reps (which has 25-32 total, in_picb, in_trin exon-blocks). Reuses results/STAR_srna_strain_wise BAMs.
#SBATCH --job-name=introncap
#SBATCH --partition=1804,2204,2004,NXFL
#SBATCH --array=1-144%24
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity/logs/introncap_%a.log
set -uo pipefail
ST=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools
DIR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity
BEDS=$DIR/beds16; OUTD=$DIR/intron_cap; mkdir -p $OUTD
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DIR/allreps_samplesheet.tsv)
strain=$(echo "$LINE"|cut -f1); tpn=$(echo "$LINE"|cut -f2); tp=$(echo "$LINE"|cut -f3); rep=$(echo "$LINE"|cut -f4); bam=$(echo "$LINE"|cut -f5)
span=$BEDS/$strain-$tp.trin_span.bed
AW='length($10)>=25 && length($10)<=32'
[ -f "$bam" ] || { echo "missing bam $bam"; exit 0; }
echo "[$(date +%T)] $strain $tpn rep$rep"
total_all=$($ST view -@4 "$bam" | wc -l)                                # ALL small-RNA alignments (any length)
in_trin_span=$($ST view -@4 -L "$span" "$bam" | awk "$AW" | wc -l)      # 25-32 in Trinity FULL SPAN (intron-incl)
echo "$strain,$tpn,$rep,$total_all,$in_trin_span" > "$OUTD/$strain-$tp.$rep.csv"
echo "[done] $strain $tpn rep$rep: total_all=$total_all in_trin_span=$in_trin_span"
