#!/bin/bash
# Step-4 (genome-anchored expression test) for ALL 16 strains: map each strain's candidates to its own genome
# (mm0) + each of the 15 other genomes (mm3), then classify into the 4 classes incl. SNP-variant. The actual
# 16x15 pairwise mapping (the pangenome approach was the shortcut; this is the full version). Array = 16.
#SBATCH --job-name=step4map16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4_16/map16_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
STAR=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/STAR
ST=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
PYCLASS=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
U=$ROOT/analysis/claude_biomni_analysis/unique_pirna; S4=$U/step4_16
SRNA="--outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 --alignIntronMax 1 \
 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 \
 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 --limitBAMsortRAM 30000000000 --runThreadN 8"
FA=$S4/$X.candidates16.fasta
[ -s "$FA" ] || { echo "[skip] $X no candidates"; exit 0; }
echo "[$X] candidates=$(grep -c '^>' "$FA") $(date)"
# (1) own genome, mm0
$STAR --genomeDir "$ROOT/results/indexs/$X" --readFilesIn "$FA" --outFilterMismatchNmax 0 $SRNA --outFileNamePrefix "$S4/$X.cand_self16."
$ST index "$S4/$X.cand_self16.Aligned.sortedByCoord.out.bam"
# (2) each other genome, mm<=3
for Y in "${STR[@]}"; do
  [ "$Y" = "$X" ] && continue
  $STAR --genomeDir "$ROOT/results/indexs/$Y" --readFilesIn "$FA" --outFilterMismatchNmax 3 $SRNA --outFileNamePrefix "$S4/$X.cand_to_$Y."
  $ST index "$S4/$X.cand_to_$Y.Aligned.sortedByCoord.out.bam"
done
"$PYCLASS" "$U/classify_step416.py" "$X"
echo "[done] $X $(date)"
