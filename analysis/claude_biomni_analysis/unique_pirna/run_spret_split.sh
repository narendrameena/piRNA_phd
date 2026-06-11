#!/bin/bash
#SBATCH --job-name=spretsplit
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/spret_split.log
set -euo pipefail
BT=/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin/bowtie
RD=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
# Stream C57BL_6NJ + CAST_EiJ EXPRESSED piRNAs (collapse, all reps) against the SPRET-candidate index.
# -v 3 (<=3 mismatches, end-to-end), --best -k 1 (best target per read). Track min mismatch per candidate.
zcat $RD/collapse/C57BL_6NJ-*.raw.fasta.gz $RD/collapse/CAST_EiJ-*.raw.fasta.gz \
  | $BT -v 3 -k 1 --best -p 8 -f $U/SPRET_candidates_idx - 2> $U/bowtie_spret.log \
  | awk -F'\t' '{mm=($8==""?0:split($8,a,",")); if(!(($3) in M)||mm<M[$3])M[$3]=mm} END{for(r in M)print r"\t"M[r]}' \
  > $U/SPRET_candidate_minmm_to_others_expressed.tsv
echo "[done] $(date)"
echo "  SPRET candidates with a <=3mm EXPRESSED match in C57/CAST: $(wc -l < $U/SPRET_candidate_minmm_to_others_expressed.tsv) / 92181"
echo "  min-mm breakdown:"; cut -f2 $U/SPRET_candidate_minmm_to_others_expressed.tsv | sort | uniq -c
grep -E "reads processed|alignment" $U/bowtie_spret.log | tail -3
