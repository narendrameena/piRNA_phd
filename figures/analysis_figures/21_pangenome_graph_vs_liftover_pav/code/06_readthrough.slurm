#!/bin/bash
#SBATCH --job-name=t21_readthrough
#SBATCH --array=0-15
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav/logs/readthrough_%a.log
# THEME 21 step 6 (independent) — piC-DoG transcriptional READTHROUGH (Konstantinidou et al. 2024). Per strain (native
# frame), a piCB cluster within 50 kb downstream of a gene is a candidate piC-DoG if WHOLE-TESTIS bulk RNA shows
# continuous signal from the gene's last exon through the gap into the cluster (gap coverage >= 10% of last-exon).
# One samtools bedcov pass per strain over all (last-exon, gap, cluster) windows pooled across tp/rep BAMs.
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; D=$T/data
ST=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/samtools
STR=(129S1_SvImJ A_J AKR_J BALB_cJ C3H_HeJ C57BL_6NJ CAST_EiJ CBA_J DBA_2J FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ PWK_PhJ SPRET_EiJ WSB_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
GFF=""; for v in 3.3 3.5 3.2; do f=$B/resources/annotation/${X}_v${v}.gff3; [ -f "$f" ] && GFF=$f && break; done
[ -z "$GFF" ] && { echo "no GFF for $X"; exit 0; }
CL=$B/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/${X}.clusters.bed
BAMS=$(ls $B/results/STAR_rna_strain_wise/${X}-*/Aligned.sortedByCoord.out.bam 2>/dev/null)
[ -z "$BAMS" ] && { echo "no bulk RNA BAMs for $X"; exit 0; }
for b in $BAMS; do [ -f "${b}.bai" ] || $ST index -@ 2 "$b"; done   # bedcov needs indexed BAMs
# gene downstream-50kb windows (GFF chrom STRAIN#1#chrN -> bare N to match clusters); carry gene name + 3' end + strand
awk -F'\t' -v OFS='\t' '$3=="gene"{c=$1; sub(/^.*#chr/,"",c); n=$9; if(n~/Name=/){sub(/.*Name=/,"",n); sub(/;.*/,"",n)} else n="."; if($7=="+"){s=$5;e=$5+50000;g3=$5}else{s=($4-1-50000<0?0:$4-1-50000);e=$4-1;g3=$4-1} if(e>s)print c,s,e,n,g3,$7}' "$GFF" | sort -k1,1 -k2,2n > $D/.down_${X}.bed
# clusters in a downstream window = candidate piC-DoG (attach the gene)
bedtools intersect -a <(awk -v OFS='\t' 'NF>=3{print $1,$2,$3}' "$CL"|sort -k1,1 -k2,2n) -b $D/.down_${X}.bed -wa -wb \
  | awk -F'\t' -v OFS='\t' -v X="$X" '$1~/^([1-9]|1[0-9]|X)$/{cc=$1;cs=$2;ce=$3; gn=$7;g3=$8;gst=$9; id=cc":"cs"-"ce"|"gn; bc=X"#1#chr"cc;
      if(gst=="+"){les=(g3-500<0?0:g3-500);lee=g3;gs=g3;ge=cs} else {les=g3;lee=g3+500;gs=ce;ge=g3}
      if(lee>les) print bc,les,lee,id"|LE"; if(ge>gs) print bc,gs,ge,id"|GAP"; if(ce>cs) print bc,cs,ce,id"|CLU"}' \
  | sort -k1,1 -k2,2n > $D/.win_${X}.bed
NW=$(wc -l < $D/.win_${X}.bed); echo "$X candidate windows: $NW"
[ "$NW" -eq 0 ] && { echo "$X no candidates"; rm -f $D/.down_${X}.bed $D/.win_${X}.bed; exit 0; }
$ST bedcov $D/.win_${X}.bed $BAMS 2>/dev/null \
  | awk -F'\t' 'BEGIN{OFS="\t"}{ns=split($4,a,"|"); id=a[1]"|"a[2]; ty=a[3]; L=$3-$2; cov=(L>0)?$5/L:0; cv[id"|"ty]=cov; seen[id]=1}
      END{print "candidate\tlast_exon_cov\tgap_cov\tcluster_cov\tratio\treadthrough";
          for(id in seen){le=cv[id"|LE"]+0; gp=cv[id"|GAP"]+0; cl=cv[id"|CLU"]+0; r=(le>0)?gp/le:0; rt=(r>=0.1&&cl>0)?"YES":"NO"; print id,le,gp,cl,r,rt}}' > $D/readthrough_${X}.tsv
echo "$X readthrough YES: $(awk -F'\t' 'NR>1&&$NF=="YES"' $D/readthrough_${X}.tsv|wc -l) / $(($(wc -l<$D/readthrough_${X}.tsv)-1))"
rm -f $D/.down_${X}.bed $D/.win_${X}.bed
