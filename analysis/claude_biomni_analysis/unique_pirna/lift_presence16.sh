#!/bin/bash
# Scale-16 pangenome phase — STEP 2c: for cross-strain LOCUS PRESENCE, halLiftover the UNION of all novel
# candidates' GRCm39 loci OUT to each strain Y. A candidate whose GRCm39 locus lifts to Y => the orthologous
# locus EXISTS in Y (homolog present); if it lifts to no other strain => strain-private locus.
# PREREQ (run once after lift_cand16.sh): cat unique16/loci/*.cand_GRCm39.bed > unique16/loci/all.cand_GRCm39.bed
# Array = 16 target strains.
#SBATCH --job-name=pirpresence16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/presence_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
Y=${1:-${STR[${SLURM_ARRAY_TASK_ID:-0}]}}   # $1 = strain (Snakemake) else SLURM array index
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
SIF=$ROOT/cactus_v2.9.3.sif
HAL=$ROOT/results/pangenome/output/mouse_17strain_pangenome.full.hal
LOCI=$ROOT/analysis/claude_biomni_analysis/unique_pirna/unique16/loci
[ -s "$LOCI/all.cand_GRCm39.bed" ] || { echo "[err] build $LOCI/all.cand_GRCm39.bed first (cat *.cand_GRCm39.bed)"; exit 1; }
singularity exec --bind /mnt "$SIF" halLiftover "$HAL" GRCm39 "$LOCI/all.cand_GRCm39.bed" "$Y" "$LOCI/present_in_$Y.bed"
echo "[done] GRCm39 -> $Y  present rows=$(wc -l < "$LOCI/present_in_$Y.bed")  $(date)"
