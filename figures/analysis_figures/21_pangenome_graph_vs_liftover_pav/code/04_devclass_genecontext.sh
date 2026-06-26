#!/bin/bash
# THEME 21 step 4 (independent of graph build) — annotate the 42,384 master piRNA-cluster loci with
#   (a) DEVELOPMENTAL CLASS from the per-tp lifted clusters (E16.5 / P12.5 / P20.5 presence, pooled over strains)
#       pachytene = P20.5-only; pre_pachytene = P12.5+P20.5 (no E16.5); hybrid = E16.5+P20.5; fetal = E16.5 w/o P20.5.
#       (full 3-timepoint scheme — the navin concept lacked E16.5 and could not resolve pre-pachytene vs hybrid.)
#   (b) GENE CONTEXT vs GRCm39 Ensembl gene models (priority 3UTR > CDS > gene_body > downstream-50kb > intergenic)
#       pre-pachytene clusters are mRNA/3'UTR/readthrough-derived; pachytene clusters are intergenic (A-MYB driven).
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
CP=$B/analysis/claude_biomni_analysis/unique_pirna/cluster_pav
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; D=$T/data
GFF=$B/results/strains_genome_annotation/Mus_musculus.GRCm39.113.chr.gff3
ML=$D/master_loci_GRCm39.bed
# ---------- (a) developmental classification ----------
for tp in 16.5dpc 12.5dpp 20.5dpp; do
  cat $CP/bytp/*.${tp}.in_GRCm39.bed | awk 'BEGIN{OFS="\t"} $1!~/_/ && NF>=3{print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -d 1000 > $D/tp_${tp}.bed
  bedtools intersect -a $ML -b $D/tp_${tp}.bed -c | awk '{print ($5>0)?1:0}' > $D/.p_${tp}
done
paste $ML $D/.p_16.5dpc $D/.p_12.5dpp $D/.p_20.5dpp | awk 'BEGIN{OFS="\t"} {
  e=$5;p1=$6;p2=$7
  if(p2&&!p1&&!e)        c="pachytene"
  else if(p1&&p2&&!e)    c="pre_pachytene"
  else if(e&&p2)         c="hybrid"
  else if(e&&!p2)        c="fetal"
  else if(p1&&!p2&&!e)   c="P12.5_only"
  else                   c="other"
  print $4,c,e p1 p2 }' > $D/master_devclass.tsv
echo "=== developmental class (master loci) ==="; cut -f2 $D/master_devclass.tsv | sort | uniq -c | sort -rn
# ---------- (b) gene context vs GRCm39 ----------
awk -F'\t' 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="gene"{c=$1;sub(/^chr/,"",c);print c,$4-1,$5,".",".",$7}' $GFF | sort -k1,1 -k2,2n > $D/g39_genes.bed
awk -F'\t' 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="three_prime_UTR"{c=$1;sub(/^chr/,"",c);print c,$4-1,$5}' $GFF | sort -k1,1 -k2,2n | bedtools merge > $D/g39_3utr.bed
awk -F'\t' 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="CDS"{c=$1;sub(/^chr/,"",c);print c,$4-1,$5}' $GFF | sort -k1,1 -k2,2n | bedtools merge > $D/g39_cds.bed
awk -F'\t' 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="gene"{c=$1;sub(/^chr/,"",c); if($7=="+"){s=$5;e=$5+50000}else{s=($4-1-50000<0?0:$4-1-50000);e=$4-1} if(e>s)print c,s,e}' $GFF | sort -k1,1 -k2,2n | bedtools merge > $D/g39_down50k.bed
( bedtools intersect -a $ML -b $D/g39_3utr.bed   -u | awk '{print $4"\t1\t3UTR"}'
  bedtools intersect -a $ML -b $D/g39_cds.bed    -u | awk '{print $4"\t2\tCDS"}'
  bedtools intersect -a $ML -b $D/g39_genes.bed  -u | awk '{print $4"\t3\tgene_body"}'
  bedtools intersect -a $ML -b $D/g39_down50k.bed -u | awk '{print $4"\t4\tdownstream"}' ) \
  | sort -k1,1 -k2,2n | awk '!seen[$1]++{print $1"\t"$3}' > $D/.ctx_best
awk 'NR==FNR{ctx[$1]=$2; next} {print $4"\t"(($4 in ctx)?ctx[$4]:"intergenic")}' $D/.ctx_best $ML > $D/master_genecontext.tsv
echo "=== gene context (master loci) ==="; cut -f2 $D/master_genecontext.tsv | sort | uniq -c | sort -rn
# ---------- join annotations ----------
paste <(sort -k1,1 $D/master_devclass.tsv) <(sort -k1,1 $D/master_genecontext.tsv | cut -f2) | sort -V > $D/master_annot.tsv 2>/dev/null || \
  join -t$'\t' <(sort -k1,1 $D/master_devclass.tsv) <(sort -k1,1 $D/master_genecontext.tsv) > $D/master_annot.tsv
rm -f $D/.p_* $D/.ctx_best
echo "wrote master_devclass.tsv master_genecontext.tsv master_annot.tsv"
echo "=== dev-class x gene-context cross-tab ==="
join -t$'\t' <(sort -k1,1 $D/master_devclass.tsv|cut -f1,2) <(sort -k1,1 $D/master_genecontext.tsv) | awk -F'\t' '{print $2"\t"$3}' | sort | uniq -c | sort -rn | head -20
