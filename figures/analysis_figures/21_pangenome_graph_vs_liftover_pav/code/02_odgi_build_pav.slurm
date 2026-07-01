#!/bin/bash
#SBATCH --job-name=t21_odgi_pav
#SBATCH --partition=TEST
#SBATCH --cpus-per-task=32
#SBATCH --mem=950G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav/logs/odgi_build_pav_%j.log
# THEME 21 step 2 — build the odgi graph from the FINAL minigraph-cactus GFA, then compute the GRAPH sequence-PAV
# (odgi pav, grouped by strain -S) at the master piRNA-cluster loci (GRCm39-anchored). Output: per-locus x 17-sample
# PAV-ratio matrix = fraction of each GRCm39 cluster-locus's graph nodes traversed by each strain (genetic presence).
# This is compared (step 3, python) against the HAL-liftover cluster-PAV (regulatory presence) -> loss vs silencing.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; W=$T/data
SIF=$B/cactus_v2.9.3.sif; S="singularity exec --bind /mnt $SIF"
GFAGZ=$B/results/pangenome/output/mouse_17strain_pangenome.gfa.gz
GFA=$W/graph.gfa; OG=$W/graph.og
echo "[$(date)] START on $(hostname)  mem=${SLURM_MEM_PER_NODE:-?}"
# 1. decompress GFA (odgi build needs plain GFAv1)
if [ ! -s "$OG" ]; then
  [ -s "$GFA" ] || { echo "[$(date)] decompress $GFAGZ"; zcat "$GFAGZ" > "$GFA"; }
  echo "[$(date)] GFA: $(ls -lh $GFA | awk '{print $5}')"
  # 2. build odgi graph (-O compact ids, -P progress)
  echo "[$(date)] odgi build ..."
  /usr/bin/time -v $S odgi build -g "$GFA" -o "$OG" -t 32 -P -O 2>&1 | tail -25
fi
echo "[$(date)] odgi build done: $(ls -lh $OG | awk '{print $5}')"
$S odgi stats -i "$OG" -S 2>/dev/null > "$W/graph_stats.txt" || true; cat "$W/graph_stats.txt" 2>/dev/null | head
# 3. graph sequence-PAV at master loci, grouped by strain (PanSN sample = 1st field)
echo "[$(date)] odgi pav at $(wc -l < $W/master_loci_pav.bed) loci ..."
$S odgi pav -i "$OG" -b "$W/master_loci_pav.bed" -S -M -t 32 > "$W/graph_pav_matrix.tsv" 2>>"$T/logs/odgi_pav.err"
echo "[$(date)] graph_pav_matrix rows: $(wc -l < $W/graph_pav_matrix.tsv)  cols: $(head -1 $W/graph_pav_matrix.tsv | awk '{print NF}')"
head -2 "$W/graph_pav_matrix.tsv" | cut -c1-200
# 4. free the 30 GB plain GFA (keep the .og for downstream untangle/snarl steps)
rm -f "$GFA"
echo "[$(date)] DONE"
