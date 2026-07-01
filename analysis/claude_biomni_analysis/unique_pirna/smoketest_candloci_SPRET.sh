#!/bin/bash
# Smoke test of the cand_loci16 path NOW (does not need DA): take 50 real SPRET 24-32 nt sequences,
# format them as candidate reads (>SPRET_EiJ|tp|seq), run the exact cand_loci16 STAR + bamtobed + chrom-strip,
# and verify (a) STAR accepts the FASTA, (b) candidate read names (with '|') survive into the BED, (c) chroms
# strip PanSN->Ensembl cleanly. Validates the new "candidate-FASTA-as-STAR-reads" path before the real run.
set -euo pipefail
ROOT=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
U=$ROOT/analysis/claude_biomni_analysis/unique_pirna
STAR=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/STAR
BT=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools
IDX=$ROOT/results/indexs/SPRET_EiJ
OUT=$U/unique16/smoketest; mkdir -p "$OUT"
TMP=$(mktemp -d -p "${TMPDIR:-/tmp}" smoke_XXXX); trap 'rm -rf "$TMP"' EXIT
set +o pipefail   # head -50 closes the pipe early -> SIGPIPE in zcat; harmless here
zcat "$ROOT/results/collapse/SPRET_EiJ-20.5dpp.1.raw.fasta.gz" | awk 'NR%2==0 && length($0)>=24 && length($0)<=32' \
  | head -50 | awk '{print ">SPRET_EiJ|20.5dpp|"$0"\n"$0}' > "$TMP/in.fasta"
set -o pipefail
echo "input candidate reads: $(grep -c '^>' "$TMP/in.fasta")"
$STAR --runThreadN 4 --genomeDir "$IDX" --readFilesIn "$TMP/in.fasta" --outFileNamePrefix "$TMP/sp." \
  --outFilterMismatchNmax 0 --outFilterMultimapNmax 800 --winAnchorMultimapNmax 1600 --seedSearchStartLmax 10 \
  --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All --scoreDelOpen -10000 --scoreInsOpen -10000 \
  --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 --limitBAMsortRAM 30000000000
$BT bamtobed -i "$TMP/sp.Aligned.sortedByCoord.out.bam" | sed "s/^SPRET_EiJ#1#chr//" > "$OUT/smoketest.bed"
nrows=$(wc -l < "$OUT/smoketest.bed")
names_ok=$(awk -F'\t' '$4 ~ /\|/' "$OUT/smoketest.bed" | wc -l)
chrom_bad=$(awk -F'\t' '$1 ~ /#/' "$OUT/smoketest.bed" | wc -l)
echo "=== BED rows=$nrows | names_with_pipe=$names_ok | rows_with_unstripped_chrom=$chrom_bad ==="
echo "--- sample BED lines ---"; head -4 "$OUT/smoketest.bed"
if [ "$nrows" -gt 0 ] && [ "$names_ok" -gt 0 ] && [ "$chrom_bad" -eq 0 ]; then echo "SMOKETEST: PASS"; else echo "SMOKETEST: FAIL"; fi
