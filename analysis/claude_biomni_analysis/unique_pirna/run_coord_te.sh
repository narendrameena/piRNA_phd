#!/bin/bash
#SBATCH --job-name=coordte
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te/coord_%a.log
set -euo pipefail
STR=(C57BL_6NJ CAST_EiJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
MM=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/minimap2
PY=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/python
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
PG=$ROOT/analysis/claude_biomni_analysis/unique_pirna/pangenome_te
GENOME=$ROOT/resources/REL-2205-Assembly/${X}_chromosomes_MT.fasta
# locate each X-private insertion sequence back in X's own genome (near-identical -> asm5), primary only
"$MM" -x asm5 --secondary=no -t 8 "$GENOME" "$PG/$X.private_insertions.fasta" > "$PG/$X.ins.paf"
# PAF -> BED of confident insertion loci (mapq>=20, aln block >=40 bp); cols: 6=target 8=tstart 9=tend 1=qname 11=alnlen 12=mapq
awk -F'\t' '$12>=20 && $11>=40 {print $6"\t"$8"\t"$9"\t"$1}' "$PG/$X.ins.paf" | sort -k1,1 -k2,2n > "$PG/$X.ins_loci.bed"
echo "[$X] insertion loci localized: $(wc -l < "$PG/$X.ins_loci.bed")"
"$PY" "$ROOT/analysis/claude_biomni_analysis/unique_pirna/coord_classify.py" "$X"
echo "[done] $X $(date)"
