#!/bin/bash
#SBATCH --job-name=t21_pav_fix
#SBATCH --partition=2204,2004,NXFL
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=6:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav/logs/odgi_pav_fix_%j.log
# Re-run odgi pav with the ACTUAL GRCm39 path names as named in graph.og (build kept GRCm39#0#<chr>#<off>, not
# GRCm39#0#<chr>). Map each master-locus chrom -> the real GRCm39 path, then odgi pav -S -M (graph sequence-PAV).
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; D=$T/data
S="singularity exec --bind /mnt $B/cactus_v2.9.3.sif"; OG=$D/graph.og
echo "[$(date)] list paths"
$S odgi paths -i "$OG" -L > $D/og_paths.txt
echo "  GRCm39 paths (sample):"; grep "^GRCm39" $D/og_paths.txt | head
echo "  n GRCm39 paths total: $(grep -c '^GRCm39' $D/og_paths.txt)  per-chrom:"; grep "^GRCm39" $D/og_paths.txt | awk -F'#' '{print $3}' | sort | uniq -c | head -25
# build chrom -> GRCm39 path map (expect ONE path per chrom for the reference)
grep "^GRCm39" $D/og_paths.txt | awk -F'#' '{print $3"\t"$0}' | sort -k1,1 -u > $D/.g39map
# rewrite pav BED with the real path name for each locus's chrom
awk 'NR==FNR{p[$1]=$2; next} ($1 in p){print p[$1]"\t"$2"\t"$3"\t"$4}' $D/.g39map $D/master_loci_GRCm39.bed > $D/master_loci_pav_fixed.bed
echo "[$(date)] pav BED rows: $(wc -l < $D/master_loci_pav_fixed.bed) / $(wc -l < $D/master_loci_GRCm39.bed) master loci"
head -1 $D/master_loci_pav_fixed.bed
echo "[$(date)] odgi pav ..."
$S odgi pav -i "$OG" -b $D/master_loci_pav_fixed.bed -S -M -t 16 > $D/graph_pav_matrix.tsv 2>$T/logs/odgi_pav_fix.err
echo "[$(date)] DONE rows=$(wc -l < $D/graph_pav_matrix.tsv) cols=$(head -1 $D/graph_pav_matrix.tsv|awk '{print NF}')"
head -2 $D/graph_pav_matrix.tsv | cut -c1-220
