#!/bin/bash
#SBATCH --job-name=collapse13
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-12
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/collapse13/c13_%a.log
set -euo pipefail
STR=(BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ PWK_PhJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
RAW=$ROOT/resources/mice16_data/srna
COL=$ROOT/results/collapse
CUT=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/cutadapt
COLLAPSE=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/collapse
# EXACT pipeline sRNA params (user-approved: reuse the validated adapter for all strains)
CUTPARAMS="--minimum-length 20 --maximum-length 36 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG"
TMP=$(mktemp -d -p "${TMPDIR:-/tmp}" c13_${X}_XXXX)
trap 'rm -rf "$TMP"' EXIT
for tp in 16.5dpc 12.5dpp 20.5dpp; do
  for r in 1 2 3; do
    lanes=$(ls $RAW/${X}-${tp}.${r}.*.fastq.gz 2>/dev/null || true)
    [ -z "$lanes" ] && { echo "[skip] no fastq for ${X}-${tp}.${r}"; continue; }
    out=$COL/${X}-${tp}.${r}.raw.fasta.gz
    [ -s "$out" ] && { echo "[exists] $out"; continue; }
    zcat $lanes > $TMP/in.fastq
    $CUT $CUTPARAMS -o $TMP/trim.fastq $TMP/in.fastq > $TMP/${X}-${tp}.${r}.cutadapt.log 2>&1
    tr=$(grep -E "Reads written|passing filters" $TMP/${X}-${tp}.${r}.cutadapt.log | head -1)
    $COLLAPSE $TMP/trim.fastq "$out"
    echo "[done] ${X}-${tp}.${r}: $(zcat "$out" | grep -c '^>') uniq seqs | $tr"
    cp $TMP/${X}-${tp}.${r}.cutadapt.log $COL/../../analysis/claude_biomni_analysis/unique_pirna/collapse13/ 2>/dev/null || true
    rm -f $TMP/in.fastq $TMP/trim.fastq
  done
done
echo "[ALLDONE] $X $(date)"
