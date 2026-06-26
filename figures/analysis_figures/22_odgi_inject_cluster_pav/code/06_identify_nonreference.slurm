#!/bin/bash
#SBATCH --job-name=t22_nonref
#SBATCH --array=0-15
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav/logs/nonref_%a.log
# THEME 22 step 6 — CATCH non-reference clusters. Re-lift each strain's piCB clusters (tagged with a per-cluster ID)
# to GRCm39 via the HAL; a native cluster that produces NO lifted interval = NON-REFERENCE (its sequence is absent
# from GRCm39, e.g. a strain-specific TE-insertion locus) -> exactly the clusters the GRCm39-anchored PAV misses.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
CP=$B/analysis/claude_biomni_analysis/unique_pirna/cluster_pav
T22=$B/figures/analysis_figures/22_odgi_inject_cluster_pav; D=$T22/data
SIF=$B/cactus_v2.9.3.sif; HAL=$B/results/pangenome/output/mouse_17strain_pangenome.full.hal
STR=(129S1_SvImJ A_J AKR_J BALB_cJ C3H_HeJ C57BL_6NJ CAST_EiJ CBA_J DBA_2J FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ PWK_PhJ SPRET_EiJ WSB_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}; mkdir -p $D/nonref
awk -v X=$X 'BEGIN{OFS="\t"}{print $1,$2,$3,X"|"NR}' $CP/$X.clusters.bed > $D/nonref/$X.named.bed
singularity exec --bind /mnt "$SIF" halLiftover "$HAL" "$X" "$D/nonref/$X.named.bed" GRCm39 "$D/nonref/$X.lifted.bed" 2>/dev/null
cut -f4 $D/nonref/$X.named.bed | sort -u > $D/nonref/$X.all_ids
cut -f4 $D/nonref/$X.lifted.bed 2>/dev/null | sort -u > $D/nonref/$X.lifted_ids
comm -23 $D/nonref/$X.all_ids $D/nonref/$X.lifted_ids > $D/nonref/$X.nonref_ids
# keep the native BED of the non-reference clusters (for TE/dev-class characterisation)
awk 'NR==FNR{nr[$1];next} ($4 in nr)' $D/nonref/$X.nonref_ids $D/nonref/$X.named.bed > $D/nonref/$X.nonref.bed
NAT=$(wc -l < $D/nonref/$X.all_ids); LIF=$(wc -l < $D/nonref/$X.lifted_ids); NR=$(wc -l < $D/nonref/$X.nonref_ids)
printf "%s\t%d\t%d\t%d\t%.1f\n" "$X" "$NAT" "$LIF" "$NR" "$(awk -v a=$NR -v b=$NAT 'BEGIN{print (b>0)?100*a/b:0}')" > $D/nonref/$X.summary
echo "$X: native=$NAT lifted=$LIF NON-REFERENCE=$NR ($(awk -v a=$NR -v b=$NAT 'BEGIN{printf "%.0f",(b>0)?100*a/b:0}')%)"
rm -f $D/nonref/$X.all_ids $D/nonref/$X.lifted_ids $D/nonref/$X.named.bed
