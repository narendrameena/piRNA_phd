#!/bin/bash
# THEME 21 step 7 (independent) — pseudogene fragments (PGFs) in pachytene clusters. Pachytene piRNA clusters can
# regulate mRNAs in trans via embedded ANTISENSE pseudogene fragments (Loubalova et al. 2025). Take pachytene/hybrid
# INTERGENIC master loci (GRCm39 frame), extract their GRCm39 sequence, and blastn vs the GRCm39 CDS set; antisense
# hits = candidate mRNA-targeting piRNA sources. All-GRCm39 frame (no cross-coordinate liftover).
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; D=$T/data
GFF=$B/results/strains_genome_annotation/Mus_musculus.GRCm39.113.chr.gff3
GENOME=$B/results/pangenome/prepared/GRCm39.fa; ML=$D/master_loci_GRCm39.bed
join -t$'\t' <(awk -F'\t' '$2=="pachytene"||$2=="hybrid"{print $1}' $D/master_devclass.tsv|sort) \
             <(awk -F'\t' '$2=="intergenic"{print $1}' $D/master_genecontext.tsv|sort) > $D/.pgf_ids
awk 'NR==FNR{ids[$1];next} ($4 in ids)' $D/.pgf_ids $ML > $D/pachytene_intergenic.bed
echo "pachytene/hybrid intergenic loci: $(wc -l < $D/pachytene_intergenic.bed)"
bedtools getfasta -fi "$GENOME" -bed $D/pachytene_intergenic.bed -nameOnly -fo $D/pachytene_intergenic.fa
# GRCm39 CDS (spliced) + transcript->gene map
[ -s $D/g39_cds_seq.fa ] || /mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/gffread -x $D/g39_cds_seq.fa -g "$GENOME" "$GFF" 2>/dev/null
awk -F'\t' '$3=="mRNA"{match($9,/ID=transcript:([^;]+)/,a); match($9,/Name=([^;]+)/,b); match($9,/Parent=gene:([^;]+)/,p); if(a[1])print a[1]"\t"(b[1]?b[1]:p[1])}' "$GFF" > $D/.tx2gene 2>/dev/null || true
/usr/bin/makeblastdb -in $D/g39_cds_seq.fa -dbtype nucl -out $D/.cds_db 2>/dev/null
/usr/bin/blastn -query $D/pachytene_intergenic.fa -db $D/.cds_db -outfmt '6 qseqid sseqid pident length sstart send evalue bitscore' \
  -evalue 1e-10 -perc_identity 70 -max_target_seqs 5 -num_threads 8 2>/dev/null \
  | awk 'BEGIN{OFS="\t"} $4>=100{o=($5<$6)?"sense":"antisense"; print $0,o}' > $D/pgf_hits.tsv
NA=$(awk -F'\t' '$NF=="antisense"' $D/pgf_hits.tsv | wc -l)
echo "PGF hits (>=70% id, >=100 bp): $(wc -l < $D/pgf_hits.tsv) | antisense (mRNA-targeting candidates): ${NA}"
echo "=== top antisense PGF target transcripts ==="
awk -F'\t' '$NF=="antisense"{print $2}' $D/pgf_hits.tsv | sort | uniq -c | sort -rn | head -12 | \
  awk 'NR==FNR{g[$1]=$2;next}{t=$2; print "    "$1"  "t"  "(t in g?g[t]:"")}' $D/.tx2gene - 2>/dev/null || \
  awk -F'\t' '$NF=="antisense"{print $2}' $D/pgf_hits.tsv | sort | uniq -c | sort -rn | head -12
rm -f $D/.cds_db.*
