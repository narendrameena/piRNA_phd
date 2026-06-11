#!/bin/bash
# 16-strain coordinate TE-driven: minimap2 each strain's private insertions back to its own genome -> insertion
# loci (PanSN) -> coord_classify16 (locus-in-private-insertion by class). Array = 16 strains. Run AFTER
# run_pangenome_insertions16 (fastas) and cand_self16 (BAMs).
#SBATCH --job-name=coordte16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te/coordte16_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
MM=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/minimap2
PY=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/python
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
U=$ROOT/analysis/claude_biomni_analysis/unique_pirna; PG=$U/pangenome_te
GENOME=$ROOT/resources/REL-2205-Assembly/${X}_chromosomes_MT.fasta
INS=$PG/ins16/$X.private_insertions.fasta
[ -s "$INS" ] || { echo "[skip] $X no private insertions"; exit 0; }
"$MM" -x asm5 --secondary=no -t 8 "$GENOME" "$INS" > "$PG/ins16/$X.ins.paf"
awk -F'\t' '$12>=20 && $11>=40 {print $6"\t"$8"\t"$9"\t"$1}' "$PG/ins16/$X.ins.paf" | sort -k1,1 -k2,2n > "$PG/ins16/$X.ins_loci.bed"
echo "[$X] insertion loci: $(wc -l < "$PG/ins16/$X.ins_loci.bed")"
"$PY" "$U/coord_classify16.py" "$X"
echo "[done] $X $(date)"
