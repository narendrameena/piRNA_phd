#!/bin/bash
#SBATCH --job-name=fc_list1_v33
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --array=0-15
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/fc_list1_v33_%a.log
set -e
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
PY=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
FC=/mnt/home3/miska/nm667/miniconda3/pkgs/subread-2.0.8-h577a1d6_0/bin/featureCounts
GFF=$ROOT/results/genicRegionGff3/intergenic/${X}_v3.3.gff3
[ -f "$GFF" ] || $PY $ROOT/workflow/scripts/python/intergenic_list1_gff.py $ROOT/resources/annotation/${X}_v3.3.gff3 "$GFF"
OUT=$ROOT/results/IntronExonfeatureCount/pacBio/list1_v3.3/$X; mkdir -p $OUT
for tp in 16.5dpc 12.5dpp 20.5dpp; do for r in 1 2 3; do
  BAM=$ROOT/results/STAR_srna_strain_wise/$X/$X-$tp.$r/Aligned.sortedByCoord.out.bam
  [ -f "$BAM" ] && $FC -M --fraction -O -t gene,lnc_RNA,intergenic -g ID -F GTF -T 4 -a "$GFF" -o $OUT/${X}-${tp}.${r}.featureCounts "$BAM" 2>>$OUT/fc.log
done; done
echo "[done] $X list1 v3.3 $(ls $OUT/*.featureCounts 2>/dev/null|wc -l) libs"
