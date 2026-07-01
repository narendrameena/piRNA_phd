#!/bin/bash
#SBATCH --job-name=halflift
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/lift_%a.log
set -euo pipefail
STR=(WSB_EiJ CAST_EiJ BALB_cJ C57BL_6NJ A_J DBA_2J NOD_ShiLtJ SPRET_EiJ AKR_J C3H_HeJ CBA_J PWK_PhJ NZO_HlLtJ FVB_NJ 129S1_SvImJ LP_J)
X=${STR[$SLURM_ARRAY_TASK_ID]}
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
SIF=$ROOT/cactus_v2.9.3.sif
HAL=$ROOT/results/pangenome/output/mouse_17strain_pangenome.full.hal
PAV=$ROOT/analysis/claude_biomni_analysis/unique_pirna/cluster_pav
# project this strain's PICB clusters into the GRCm39 common frame through the cactus graph/HAL
singularity exec --bind /mnt "$SIF" halLiftover "$HAL" "$X" "$PAV/$X.clusters.bed" GRCm39 "$PAV/$X.in_GRCm39.bed"
src=$(wc -l < "$PAV/$X.clusters.bed"); lifted=$(cut -f4 "$PAV/$X.in_GRCm39.bed" 2>/dev/null | sort -u | wc -l)
echo "[$X] source clusters=$src | lifted-interval rows=$(wc -l < "$PAV/$X.in_GRCm39.bed")"
echo "[done] $X $(date)"
