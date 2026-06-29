#!/bin/bash
#SBATCH --job-name=GRCm39_unmasked_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --partition=TEST
#SBATCH --output=log/GRCm39_star_index_%j.log
#SBATCH --error=log/GRCm39_star_index_%j.log

set -euo pipefail

###############################################################################
# Rebuild GRCm39 STAR index from the correct unmasked genome.
#
# Background: results/ref_genome/GRCm39.106.fasta was previously a hard-masked
# (dna_rm) file. It has been replaced with a symlink to the unmasked
# primary assembly. This script rebuilds the STAR index from the unmasked file.
#
# Input:  results/ref_genome/GRCm39.106.fasta  (symlink → unmasked dna)
# Output: results/indexs/GRCm39.106/           (rebuilt from unmasked genome)
###############################################################################

BASE=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
GENOME=${BASE}/results/ref_genome/GRCm39.106.fasta
INDEX_DIR=${BASE}/results/indexs/GRCm39.106
STAR=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/.snakemake/conda/3ed2c25f2a7c3a9c29f4adb59cbb09c0_/bin/STAR
THREADS=16

timestamp() { date '+%Y-%m-%dT%H:%M:%S'; }

echo "$(timestamp) === GRCm39 STAR index rebuild (unmasked) ==="
echo "Genome:    $GENOME"
echo "Index dir: $INDEX_DIR"
echo "STAR:      $($STAR --version 2>&1 | head -1)"
echo ""

# Verify genome is unmasked before doing anything
HEADER=$(grep '^>' "$GENOME" | head -1)
echo "$(timestamp) Genome header: $HEADER"
if echo "$HEADER" | grep -qE 'dna_rm|dna_sm'; then
    echo "FATAL: Genome appears masked (dna_rm/dna_sm in header). Aborting."
    exit 1
fi
echo "$(timestamp) Genome masking check: PASS (unmasked)"
echo ""

# Rename the old masked index to backup (preserves data)
if [ -d "${INDEX_DIR}" ]; then
    BACKUP="${INDEX_DIR}.dna_rm_backup"
    if [ ! -d "${BACKUP}" ]; then
        echo "$(timestamp) Renaming old masked index to ${BACKUP}"
        mv "${INDEX_DIR}" "${BACKUP}"
    else
        echo "$(timestamp) Backup ${BACKUP} already exists — removing old index dir"
        rm -rf "${INDEX_DIR}"
    fi
fi

mkdir -p "${INDEX_DIR}"

echo "$(timestamp) Running STAR genomeGenerate ..."
echo ""

$STAR \
    --runMode genomeGenerate \
    --runThreadN ${THREADS} \
    --genomeFastaFiles "${GENOME}" \
    --genomeDir "${INDEX_DIR}" \
    --limitGenomeGenerateRAM 96000000000 \
    --outTmpDir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/tmp/GRCm39_star_index_tmp

echo ""
echo "$(timestamp) STAR done."

# Verify output
echo ""
echo "$(timestamp) === Index verification ==="
echo "Files:"
ls -lh "${INDEX_DIR}/" | head -10
echo ""
echo "genomeParameters.txt (key lines):"
grep -E "genomeFastaFiles|genomeSAindexNbases|versionGenome" "${INDEX_DIR}/genomeParameters.txt"

echo ""
echo "$(timestamp) === ALL DONE ==="
echo "New index: ${INDEX_DIR}"
echo "Old masked index backed up at: ${INDEX_DIR}.dna_rm_backup"
