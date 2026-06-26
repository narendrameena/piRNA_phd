#!/bin/bash
#SBATCH --job-name=t22_pav_cluster
#SBATCH --partition=TEST,NXFL,2204,2004
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --time=6:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav/logs/pav_cluster_%j.log
# THEME 22 step 4 (2nd pivot — overlap too slow, untangle crashes on the injected graph's step index). Use odgi pav
# (fast + worked in theme 21, no step index) with -p path-groups: group the injected cluster-paths by STRAIN, and at
# each GRCm39 master locus report the fraction of the locus's nodes covered by each strain's CLUSTER-paths = graph
# cluster-presence. Reuses theme-21 master_loci_pav_fixed.bed. Compared (step 5) to the HAL-liftover cluster-PAV.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T21=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav
T22=$B/figures/analysis_figures/22_odgi_inject_cluster_pav; D=$T22/data
S="singularity exec --bind /mnt $B/cactus_v2.9.3.sif"; INJ=$D/graph_inj.og
awk -F'|' 'BEGIN{OFS="\t"}{print $0,$2}' "$D/cluster_paths.txt" > "$D/cluster_to_strain.tsv"
echo "[$(date)] path-groups: $(wc -l < $D/cluster_to_strain.tsv) cluster-paths -> $(cut -f2 $D/cluster_to_strain.tsv|sort -u|wc -l) strains | loci: $(wc -l < $T21/data/master_loci_pav_fixed.bed)"
echo "[$(date)] odgi pav -p (cluster-paths grouped by strain) ..."
$S odgi pav -i "$INJ" -b "$T21/data/master_loci_pav_fixed.bed" -p "$D/cluster_to_strain.tsv" -M -t 32 > "$D/graph_cluster_pav_matrix.tsv" 2>"$T22/logs/pav_cluster.err"
echo "[$(date)] DONE rows=$(wc -l < $D/graph_cluster_pav_matrix.tsv) cols=$(head -1 $D/graph_cluster_pav_matrix.tsv|awk '{print NF}')"
head -2 "$D/graph_cluster_pav_matrix.tsv" | cut -c1-220
