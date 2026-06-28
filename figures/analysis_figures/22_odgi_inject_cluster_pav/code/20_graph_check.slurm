#!/bin/bash
#SBATCH --job-name=t22_gcheck
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --partition=TEST,NXFL,2204,2004
#SBATCH --time=4:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav/logs/graph_check.log
# THEME 22 step 20 — GRAPH-NATIVE cross-check of step-18's deletion/present classification of the shared-subset loci.
# odgi pav on graph_inj.og: each shared-subset cluster-path (CLU|strain|id, in cluster-path frame [0,len]) x path-groups
# {GRCm39 + 16 strains}. GRCm39's coverage of a cluster's nodes = does GRCm39 carry that sequence in the GRAPH.
# Expectation: the 15 'PRESENT_in_GRCm39' (minimap2) loci -> high GRCm39 coverage; the 6 'C57BL6J_DELETION' -> ~0.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T22=$B/figures/analysis_figures/22_odgi_inject_cluster_pav; D=$T22/data
SIF=$B/cactus_v2.9.3.sif; INJ=$D/graph_inj.og
S="singularity exec --bind /mnt $SIF"
echo "[$(date)] odgi pav: GRCm39 + strain coverage of $(wc -l < $D/graph_check_loci.bed) shared-subset cluster-loci"
$S odgi pav -i "$INJ" -b "$D/graph_check_loci.bed" -p "$D/graph_check_groups.tsv" -M -t 32 > "$D/graph_check_pav.tsv" 2>"$T22/logs/graph_check.err"
echo "[$(date)] done: $(wc -l < $D/graph_check_pav.tsv) rows; header:"; head -1 "$D/graph_check_pav.tsv"
