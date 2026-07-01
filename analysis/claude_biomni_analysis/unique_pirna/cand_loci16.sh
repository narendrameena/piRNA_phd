#!/bin/bash
# Scale-16 pangenome phase — STEP 2a: map each strain's NOVEL-SEQUENCE candidates (from classify_unique16.py
# -> {X}.{tp}.novel.fasta) to its OWN unmasked genome with the VALIDATED sRNA piRNA STAR params (0 mismatch,
# 800 multimappers, EndToEnd) -> production-loci BED. Chrom stripped PanSN ({X}#1#chrN) -> HAL Ensembl (N)
# for halLiftover. Run AFTER classify_unique16.py. Array = 16 strains.
#SBATCH --job-name=pircandloci16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/candloci_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${1:-${STR[${SLURM_ARRAY_TASK_ID:-0}]}}   # $1 = strain (Snakemake) else SLURM array index
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
U=$ROOT/analysis/claude_biomni_analysis/unique_pirna
STAR=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/STAR
ST=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
BT=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools
IDX=$ROOT/results/indexs/$X
OUT=$U/unique16/loci; mkdir -p "$OUT"
TMP=$(mktemp -d -p "${TMPDIR:-/tmp}" candloci_${X}_XXXX); trap 'rm -rf "$TMP"' EXIT
cat $U/unique16/$X.*.novel.fasta > "$TMP/in.fasta" 2>/dev/null || { echo "[skip] no novel fasta for $X"; exit 0; }
[ -s "$TMP/in.fasta" ] || { echo "[skip] empty $X"; exit 0; }
$STAR --runThreadN 8 --genomeDir "$IDX" --readFilesIn "$TMP/in.fasta" --outFileNamePrefix "$TMP/$X." \
  --outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 \
  --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 \
  --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 --limitBAMsortRAM 50000000000
# BAM -> BED; read name = candidate id (X|tp|seq); strip PanSN {X}#1#chr -> Ensembl chrom (1..19,X,Y,MT)
$BT bamtobed -i "$TMP/$X.Aligned.sortedByCoord.out.bam" | sed "s/^${X}#1#//; s/^chr//" > "$OUT/$X.cand_loci.ens.bed"
echo "[done] $X loci=$(wc -l < "$OUT/$X.cand_loci.ens.bed") $(date)"
