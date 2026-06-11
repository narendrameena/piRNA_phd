#!/bin/bash
# Lift per-timepoint expr+strand PICB clusters to GRCm39 (FPM in col4, strand col6 preserved/flipped). Array=48.
#SBATCH --job-name=halflift_expr
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --array=0-47
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/bytp/lifte_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
TPS=(16.5dpc 12.5dpp 20.5dpp)
i=$SLURM_ARRAY_TASK_ID; X=${STR[$((i/3))]}; TP=${TPS[$((i%3))]}
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
SIF=$ROOT/cactus_v2.9.3.sif; HAL=$ROOT/results/pangenome/output/mouse_17strain_pangenome.full.hal
PAV=$ROOT/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/bytp
BED=$PAV/$X.$TP.expr.clusters.bed
[ -s "$BED" ] || { echo "[skip] $X $TP"; exit 0; }
singularity exec --bind /mnt "$SIF" halLiftover "$HAL" "$X" "$BED" GRCm39 "$PAV/$X.$TP.expr.in_GRCm39.bed"
echo "[$X $TP] src=$(wc -l < "$BED") lifted=$(wc -l < "$PAV/$X.$TP.expr.in_GRCm39.bed") $(date)"
