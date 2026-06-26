#!/bin/bash
# THEME 21 step 1 — define MASTER piRNA-cluster loci (GRCm39 frame) and the HAL-LIFTOVER cluster-PAV.
# Master loci = union of all 16 strains' lifted+merged piCB clusters ({strain}.GRCm39.merged.bed), merged within 1 kb.
# Liftover cluster-PAV = per master locus, which strains have a lifted piCB cluster overlapping it (>=50% of the
# strain's lifted cluster) -> n_strains -> core(16)/dispensable(2-15)/private(1). This is the EXISTING-method baseline
# that the graph (odgi pav, step 3) is compared against. Also count NON-lifting clusters per strain (clusters that
# halLiftover could not place on GRCm39 = candidate graph-only / non-reference loci).
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
CP=$B/analysis/claude_biomni_analysis/unique_pirna/cluster_pav
T=$B/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav; D=$T/data
STRAINS=(129S1_SvImJ A_J AKR_J BALB_cJ C3H_HeJ C57BL_6NJ CAST_EiJ CBA_J DBA_2J FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ PWK_PhJ SPRET_EiJ WSB_EiJ)
# --- master loci (GRCm39 frame) ---
cat ${CP}/*.GRCm39.merged.bed | awk 'BEGIN{OFS="\t"} $1!~/_/ && NF>=3 {print $1,$2,$3}' | sort -k1,1 -k2,2n \
  | bedtools merge -d 1000 | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"L"NR}' > ${D}/master_loci_GRCm39.bed
NM=$(wc -l < ${D}/master_loci_GRCm39.bed); echo "master loci: ${NM}"
# odgi pav BED (path name = GRCm39#0#chrom)
awk 'BEGIN{OFS="\t"}{print "GRCm39#0#"$1,$2,$3,$4}' ${D}/master_loci_GRCm39.bed > ${D}/master_loci_pav.bed
# --- liftover cluster-PAV matrix (locus x strain, 1 if a strain lifted-cluster covers >=50% of itself in the locus) ---
HDR="locus_id\tchrom\tstart\tend"; for X in "${STRAINS[@]}"; do HDR="${HDR}\t${X}"; done
echo -e "${HDR}" > ${D}/liftover_pav_matrix.tsv
PRES=${D}/.pres; mkdir -p ${PRES}
for X in "${STRAINS[@]}"; do
  bedtools intersect -a ${D}/master_loci_GRCm39.bed -b ${CP}/${X}.GRCm39.merged.bed -c | awk '{print ($5>0)?1:0}' > ${PRES}/${X}.col
done
paste ${D}/master_loci_GRCm39.bed $(for X in "${STRAINS[@]}"; do echo ${PRES}/${X}.col; done) \
  | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s",$4,$1,$2,$3; for(i=5;i<=NF;i++) printf "\t%s",$i; print ""}' >> ${D}/liftover_pav_matrix.tsv
rm -rf ${PRES}
# --- summary: core/dispensable/private by liftover ---
tail -n+2 ${D}/liftover_pav_matrix.tsv | awk 'BEGIN{c=0;d=0;p=0} {s=0; for(i=5;i<=NF;i++) s+=$i; if(s>=16)c++; else if(s==1)p++; else if(s>=2)d++} END{print "  liftover core(16): "c"  dispensable(2-15): "d"  private(1): "p}'
# --- non-lifting clusters per strain (candidate non-reference / graph-only) ---
echo -e "strain\tn_clusters_native\tn_clusters_lifted\tn_nonlifting\tpct_nonlifting" > ${D}/nonlifting_per_strain.tsv
for X in "${STRAINS[@]}"; do
  nc=$(wc -l < ${CP}/${X}.clusters.bed); nl=$(wc -l < ${CP}/${X}.GRCm39.merged.bed); nn=$((nc-nl>0?nc-nl:0))
  printf "%s\t%d\t%d\t%d\t%.1f\n" "$X" "$nc" "$nl" "$nn" "$(awk -v a=$nn -v b=$nc 'BEGIN{print b>0?100*a/b:0}')" >> ${D}/nonlifting_per_strain.tsv
done
echo "wrote: master_loci_GRCm39.bed master_loci_pav.bed liftover_pav_matrix.tsv nonlifting_per_strain.tsv"
column -t ${D}/nonlifting_per_strain.tsv | head -20
