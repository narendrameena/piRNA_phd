#!/bin/bash
# THEME 21 step 5 (independent) — what STRUCTURAL VARIANTS drive the liftover-variable (dispensable/private) piRNA
# cluster loci? Intersect the precomputed deconstructed pangenome VCF (vg deconstruct vs GRCm39, 16 genotypes) with the
# variable loci, keep SV-sized alleles (>=300 bp ref/alt length diff), and type them by blastn vs the mouse TE consensus
# library (teref.mouse.fa). TE-insertion SVs at these loci = the genetic origin of strain-variable piRNA clusters.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; D=$T/data
VCF=$B/results/pangenome/output/mouse_17strain_pangenome.vcf.gz
TEREF=$B/resources/tldr/teref.mouse.fa; ML=$D/master_loci_GRCm39.bed
# variable loci (liftover n_strains < 16) with their class
awk -F'\t' 'NR>1{s=0;for(i=5;i<=NF;i++)s+=$i; if(s<16) print $2"\t"$3"\t"$4"\t"$1"\t"((s==1)?"private":"dispensable")}' $D/liftover_pav_matrix.tsv | sort -k1,1 -k2,2n > $D/variable_loci.bed
echo "variable loci (dispensable+private): $(wc -l < $D/variable_loci.bed)"
# SVs overlapping variable loci
bedtools intersect -a "$VCF" -b $D/variable_loci.bed -u > $D/sv_at_variable.vcf 2>/dev/null
echo "VCF records at variable loci: $(grep -vc '^#' $D/sv_at_variable.vcf || echo 0)"
# SV-sized (>=300 bp) -> larger allele to FASTA
awk '/^#/{next}{n=split($5,a,","); rl=length($4); for(i=1;i<=n;i++){al=length(a[i]); d=(al>rl)?al-rl:rl-al; if(d>=300){seq=(al>rl)?a[i]:$4; print ">"$1"_"$2"_a"i"_"d"bp\n"seq}}}' $D/sv_at_variable.vcf > $D/sv_alleles.fa
echo "SV-sized alleles (>=300 bp): $(grep -c '^>' $D/sv_alleles.fa || echo 0)"
# type by blastn vs TE consensus
makeblastdb -in "$TEREF" -dbtype nucl -out $D/.teref_db 2>/dev/null
blastn -query $D/sv_alleles.fa -db $D/.teref_db -outfmt '6 qseqid sseqid pident length evalue bitscore' -evalue 1e-5 -perc_identity 80 -max_target_seqs 1 -num_threads 8 2>/dev/null \
  | sort -k1,1 -k6,6gr | awk '!seen[$1]++' > $D/sv_te_hits.tsv
NT=$(wc -l < $D/sv_te_hits.tsv); NA=$(grep -c '^>' $D/sv_alleles.fa || echo 0)
echo "SV alleles typed as TE: ${NT} / ${NA}"
echo "=== TE family at variable piRNA-cluster SVs ==="; cut -f2 $D/sv_te_hits.tsv | sort | uniq -c | sort -rn
rm -f $D/.teref_db.*
