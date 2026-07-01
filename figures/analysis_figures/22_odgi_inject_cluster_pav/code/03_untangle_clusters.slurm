#!/bin/bash
#SBATCH --job-name=t22_untangle
#SBATCH --partition=TEST,NXFL,2204,2004
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav/logs/untangle_%j.log
# THEME 22 step 3 (pivot — odgi overlap was O(n) per path, ~300h for 365k paths). Project the injected cluster-paths
# onto GRCm39 in BATCH with odgi untangle -> each cluster gets graph-derived GRCm39 coords (BEDPE). Cluster-PAV =
# group cluster-paths by GRCm39 overlap -> strains/locus; clusters that DON'T project = non-reference (graph-private).
# This is the GRAPH projection (vs the HAL-liftover projection) of the same clusters.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T21=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav
T22=$B/figures/analysis_figures/22_odgi_inject_cluster_pav; D=$T22/data
S="singularity exec --bind /mnt $B/cactus_v2.9.3.sif"; INJ=$D/graph_inj.og
grep "^GRCm39" $T21/data/og_paths.txt > $D/grcm39_paths.txt
echo "[$(date)] target GRCm39 paths: $(wc -l < $D/grcm39_paths.txt) | query cluster-paths: $(wc -l < $D/cluster_paths.txt)"
echo "[$(date)] odgi untangle ..."
$S odgi untangle -i "$INJ" -Q "$D/cluster_paths.txt" -R "$D/grcm39_paths.txt" -t 32 -P > "$D/untangle.bedpe" 2>"$T22/logs/untangle.err"
echo "[$(date)] DONE — projections: $(wc -l < $D/untangle.bedpe)"; head -3 "$D/untangle.bedpe"
