#!/bin/bash
#SBATCH --job-name=senseanti
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/sense_antisense/sa_%a.log
set -euo pipefail
STR=(C57BL_6NJ CAST_EiJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
SA=$ROOT/analysis/claude_biomni_analysis/unique_pirna/sense_antisense
PY=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/python
F=$ROOT/resources/repeatMasker/${X}_chromosomes_MT_unplaced.fasta.out
# RM .out -> stranded TE BED (PanSN coords): col5=chrom col6=begin col7=end col9=strand(+/C) col10=name col11=class/family
awk 'NR>3 && NF>=11 {st=($9=="C")?"-":"+"; print $5"\t"($6-1)"\t"$7"\t"$10"|"$11"\t0\t"st}' "$F" \
  | sort -k1,1 -k2,2n > "$SA/$X.TE_stranded.bed"
echo "[$X] stranded TE annotations: $(wc -l < "$SA/$X.TE_stranded.bed")"
"$PY" "$ROOT/analysis/claude_biomni_analysis/unique_pirna/sense_antisense.py" "$X"
echo "[done] $X $(date)"
