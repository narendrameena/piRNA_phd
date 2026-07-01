#!/bin/bash
# Scale-16 pangenome phase — STEP 2b: halLiftover each strain's candidate production loci from the strain's
# own genome to the GRCm39 common frame, THROUGH the minigraph-cactus pangenome HAL (same tool/graph that
# built the cluster PAV). Run AFTER cand_loci16.sh. Array = 16 source strains.
#SBATCH --job-name=pirliftcand16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/liftcand_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${1:-${STR[${SLURM_ARRAY_TASK_ID:-0}]}}   # $1 = strain (Snakemake) else SLURM array index
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
SIF=$ROOT/cactus_v2.9.3.sif
HAL=$ROOT/results/pangenome/output/mouse_17strain_pangenome.full.hal
LOCI=$ROOT/analysis/claude_biomni_analysis/unique_pirna/unique16/loci
[ -s "$LOCI/$X.cand_loci.ens.bed" ] || { echo "[skip] no loci bed for $X"; exit 0; }
singularity exec --bind /mnt "$SIF" halLiftover "$HAL" "$X" "$LOCI/$X.cand_loci.ens.bed" GRCm39 "$LOCI/$X.cand_GRCm39.bed"
echo "[done] $X  src=$(wc -l < "$LOCI/$X.cand_loci.ens.bed") -> GRCm39 rows=$(wc -l < "$LOCI/$X.cand_GRCm39.bed")  $(date)"
