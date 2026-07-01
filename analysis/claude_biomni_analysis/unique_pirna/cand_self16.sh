#!/bin/bash
# Map ALL candidates (every class) of each strain to its OWN genome (PanSN frame) -> cand_self16 BAM, so the
# coordinate TE-driven test can score locus-in-private-insertion per class (incl. the expressed-elsewhere null).
# Same validated piRNA STAR params as cand_loci16. Array = 16 strains.
#SBATCH --job-name=candself16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-15
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te/candself16_%a.log
set -euo pipefail
STR=(C57BL_6NJ BALB_cJ A_J FVB_NJ C3H_HeJ LP_J 129S1_SvImJ DBA_2J AKR_J CBA_J NZO_HlLtJ NOD_ShiLtJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
U=$ROOT/analysis/claude_biomni_analysis/unique_pirna
STAR=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/STAR
ST=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/samtools
PYBIO=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
IDX=$ROOT/results/indexs/$X
OUT=$U/cand_self16; mkdir -p "$OUT"
TMP=$(mktemp -d -p "${TMPDIR:-/tmp}" candself_${X}_XXXX); trap 'rm -rf "$TMP"' EXIT
"$PYBIO" -c "import pandas as pd; d=pd.read_csv('$U/unique16/final_classified.csv.gz'); d=d[d.strain=='$X']; f=open('$TMP/in.fa','w'); [f.write('>$X|%s|%s\n%s\n'%(t,s,s)) for t,s in zip(d.timepoint,d.sequence)]; f.close()"
[ -s "$TMP/in.fa" ] || { echo "[skip] $X no candidates"; exit 0; }
"$STAR" --runThreadN 8 --genomeDir "$IDX" --readFilesIn "$TMP/in.fa" --outFileNamePrefix "$TMP/$X." \
  --outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 \
  --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes NH --scoreDelOpen -10000 --scoreInsOpen -10000 \
  --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 --limitBAMsortRAM 50000000000
mv "$TMP/$X.Aligned.sortedByCoord.out.bam" "$OUT/$X.cand_self16.bam"
"$ST" index "$OUT/$X.cand_self16.bam"
echo "[done] $X loci-mapped $(date)"
