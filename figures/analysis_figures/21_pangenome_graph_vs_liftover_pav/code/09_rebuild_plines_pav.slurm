#!/bin/bash
#SBATCH --job-name=t21_rebuild_pav
#SBATCH --partition=TEST
#SBATCH --cpus-per-task=32
#SBATCH --mem=950G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav/logs/rebuild_pav_%j.log
# FIX: odgi build reads GFAv1 P-lines only; the cactus GFA uses W-lines -> 0 paths imported. Convert the gbz to a
# P-line GFA (vg convert -W), rebuild graph.og WITH paths, then odgi pav (graph sequence-PAV) at the master loci.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; D=$T/data
S="singularity exec --bind /mnt $B/cactus_v2.9.3.sif"
GBZ=$B/results/pangenome/output/mouse_17strain_pangenome.gbz
GFA=$D/graph_plines.gfa; OG=$D/graph.og
echo "[$(date)] vg convert gbz -> P-line GFA (-W)"
[ -s "$GFA" ] || $S vg convert "$GBZ" -f -W -t 32 > "$GFA" 2>$T/logs/vg_convert.err
echo "[$(date)] GFA $(ls -lh $GFA|awk '{print $5}')  P-lines=$(grep -c '^P' $GFA 2>/dev/null||echo ?)  W-lines=$(grep -c '^W' $GFA 2>/dev/null||echo ?)"
echo "[$(date)] odgi build (with paths)"; rm -f "$OG"
$S odgi build -g "$GFA" -o "$OG" -t 32 -P -O 2>&1 | tail -4
echo "[$(date)] stats (length nodes edges paths steps):"; $S odgi stats -i "$OG" -S 2>/dev/null | tee $D/graph_stats.txt
$S odgi paths -i "$OG" -L > $D/og_paths.txt 2>/dev/null
echo "[$(date)] total paths=$(wc -l < $D/og_paths.txt)  GRCm39 sample:"; grep "^GRCm39" $D/og_paths.txt | head -3
# build pav BED with the ACTUAL GRCm39 path names (one path per chrom expected)
grep "^GRCm39" $D/og_paths.txt | awk -F'#' '{print $3"\t"$0}' | sort -k1,1 -u > $D/.g39map
awk 'NR==FNR{p[$1]=$2;next} ($1 in p){print p[$1]"\t"$2"\t"$3"\t"$4}' $D/.g39map $D/master_loci_GRCm39.bed > $D/master_loci_pav_fixed.bed
echo "[$(date)] pav BED rows=$(wc -l < $D/master_loci_pav_fixed.bed)/$(wc -l < $D/master_loci_GRCm39.bed)  e.g.: $(head -1 $D/master_loci_pav_fixed.bed)"
echo "[$(date)] odgi pav"
$S odgi pav -i "$OG" -b $D/master_loci_pav_fixed.bed -S -M -t 32 > $D/graph_pav_matrix.tsv 2>$T/logs/odgi_pav2.err
echo "[$(date)] DONE matrix rows=$(wc -l < $D/graph_pav_matrix.tsv) cols=$(head -1 $D/graph_pav_matrix.tsv|awk '{print NF}')"
head -2 $D/graph_pav_matrix.tsv | cut -c1-200
rm -f "$GFA"
