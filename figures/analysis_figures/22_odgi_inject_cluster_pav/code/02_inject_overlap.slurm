#!/bin/bash
#SBATCH --job-name=t22_inject
#SBATCH --partition=NXFL,2204,2004
#SBATCH --cpus-per-task=32
#SBATCH --mem=380G
#SBATCH --time=16:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav/logs/inject_overlap_%j.log
# THEME 22 step 2 — inject the 364,924 native piRNA clusters as paths onto the pangenome graph (reuses theme-21
# graph.og, which has the strain paths), then odgi overlap to find which strains' cluster-paths co-locate in graph
# space (= reference-free cross-strain cluster sharing).
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T21=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav
T22=$B/figures/analysis_figures/22_odgi_inject_cluster_pav; D=$T22/data
S="singularity exec --bind /mnt $B/cactus_v2.9.3.sif"
OG=$T21/data/graph.og; INJ=$D/graph_inj.og
echo "[$(date)] inject $(wc -l < $D/inject.bed) cluster records onto $(ls -lh $OG|awk '{print $5}') graph"
$S odgi inject -i "$OG" -b "$D/inject.bed" -o "$INJ" -t 32 -P 2>&1 | tail -3
echo "[$(date)] injected: $($S odgi paths -i "$INJ" -L 2>/dev/null | grep -c '^CLU') CLU paths (expect 364924)"
echo "[$(date)] odgi overlap — cross-cluster co-location"
$S odgi overlap -i "$INJ" -R "$D/cluster_paths.txt" -s "$D/cluster_paths.txt" -t 32 > "$D/cluster_overlap.tsv" 2>"$T22/logs/overlap.err"
echo "[$(date)] DONE — overlap rows: $(wc -l < $D/cluster_overlap.tsv)"; head -3 "$D/cluster_overlap.tsv"
