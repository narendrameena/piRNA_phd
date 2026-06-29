#!/bin/bash
set -euo pipefail

# --- Configuration ---
# List of sample groups (timepoints) 
samples=("C57BL_6-20.5dpp" "C57BL_6-16.5dpc" "C57BL_6-12.5dpp")
# Replicates per sample group (assumes replicates 1,2,3 exist)
reps=("1" "2" "3")
# Directories
RAW_DIR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/mice16_data/rna"               # Directory where copied FASTQ files are located (update as needed)
DOWNSAMPLE_DIR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/downsample"
TRINITY_DIR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/trinity"
LOG_DIR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/log"
THREADS_DOWNSAMPLE=2
THREADS_TRINITY=50
# Other parameters
DOWNSAMPLE_COUNT=50000000
SEQTK_SEED=123
MAX_MEMORY="90G"          # Trinity max memory parameter

# --- Create output directories ---
mkdir -p "$DOWNSAMPLE_DIR"
mkdir -p "$TRINITY_DIR"
mkdir -p "$LOG_DIR/downsample"
mkdir -p "$LOG_DIR/trinity"


# --- Trinity Assembly ---
# For each sample group, run Trinity using the downsampled files.
for sample in "${samples[@]}"; do
    echo "Running Trinity for sample group: ${sample}"
    
    # Input files from downsample step
    in_r1="${DOWNSAMPLE_DIR}/${sample}_r1.fastq.gz"
    in_r2="${DOWNSAMPLE_DIR}/${sample}_r2.fastq.gz"
    
    # Trinity output directory and final FASTA file
    tri_out_dir="${TRINITY_DIR}/trinity-${sample}"
    final_fasta="${TRINITY_DIR}/trinity-${sample}.Trinity.fasta"
    
    # Log files for Trinity
    log_tri="${LOG_DIR}/trinity/${sample}.out"
    err_tri="${LOG_DIR}/trinity/${sample}.err"
    
    # Create Trinity output directory if not exists
    mkdir -p "$tri_out_dir"
    
    echo "  Starting Trinity assembly for ${sample} at $(date)" | tee -a "$log_tri"
    
    Trinity --normalize_reads --min_kmer_cov 2 --seqType fq \
            --left "$in_r1" --right "$in_r2" \
            --output "$tri_out_dir" --SS_lib_type RF \
            --CPU ${THREADS_TRINITY} --max_memory ${MAX_MEMORY} \
            >> "$log_tri" 2>> "$err_tri"
    
    echo "  Trinity assembly completed for ${sample}" | tee -a "$log_tri"
    
    # Move the Trinity.fasta file to the final location.
    if [ -f "${tri_out_dir}/Trinity.fasta" ]; then
        mv "${tri_out_dir}/Trinity.fasta" "$final_fasta"
        echo "  Trinity.fasta moved to ${final_fasta}" | tee -a "$log_tri"
    else
        echo "Error: Trinity.fasta not found in ${tri_out_dir}" | tee -a "$err_tri"
        exit 1
    fi
    
    # Optionally clean up Trinity intermediate files:
    # Create an empty directory to sync (this is similar to the Snakefile approach)
    empty_dir="${tri_out_dir}_empty"
    mkdir -p "$empty_dir"
    rsync -a --delete "$empty_dir"/ "$tri_out_dir"/ >> "$log_tri" 2>> "$err_tri"
    rmdir "$tri_out_dir" "$empty_dir"
    
    echo "Trinity cleanup completed for ${sample}" | tee -a "$log_tri"
done

echo "Workflow completed successfully."