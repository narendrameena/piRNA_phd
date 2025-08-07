###############################################################################
# FILE: Snakefile
# 
# A simple pipeline that:
#   1) Reads sample info from config.yaml
#   2) Separates single-end sRNA vs. paired-end RNA samples
#   3) Runs Cutadapt (single or paired)
#   4) Aligns reads with STAR (with sRNA- or default-specific params)
###############################################################################

################################################################################
# 1) Load Configuration
################################################################################
configfile: "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/config/black_6.yaml"

# Extract all sample names from config
SAMPLES = list(config["samples"].keys())

# Separate them by type (RNA vs. sRNA)
SAMPLES_SRNA = [s for s in SAMPLES if config["samples"][s]["type"] == "sRNA"]
SAMPLES_RNA  = [s for s in SAMPLES if config["samples"][s]["type"] == "RNA"]

################################################################################
# 2) Global Parameters
################################################################################
RESULTS_DIR = config['resultDir']   # Where trimmed FASTQs/BAMs will go
LOGS_DIR    = config['logDir']     # Where logs will go
STAR_INDEX  = config["star_index"]  # From config.yaml

# Extract adapter(s) from config, if needed
ADAPTER_SRNA_3P = config['params']['srna']['cutadapt']

# STAR alignment parameters
STAR_PARAMS_SRNA = " ".join(config['params']['srna']['STAR'])
STAR_PARAMS_RNA  = " ".join(config['params']['rna']['STAR'])

################################################################################
# 3) rule all
################################################################################

################################################################################
# 4) Cutadapt (Single-end) for sRNA
################################################################################
rule cutadapt_srna:
    """
    Adapter trim single-end sRNA reads using Cutadapt.
    Pulls 3' adapter from config["adapters"]["sRNA_3prime"].
    """
    input:
        lambda wc: config["samples"][wc.sample]["R1"]
    output:
        fastq = temp(f"{RESULTS_DIR}/cutadapt_srna/{{sample}}.fastq")
    log:
        f"{LOGS_DIR}/cutadapt_srna/{{sample}}.log"
    threads: 1
    shell:
        """
        cutadapt \
            --minimum-length 20 \
            --maximum-length 36 \
            --discard-untrimmed \
            -a {ADAPTER_SRNA_3P} \
            -o {output.fastq} {input} \
        2> {log}
        """

################################################################################
# 5) Cutadapt (Paired-end) for RNA
################################################################################
rule cutadapt_rna:
    """
    Adapter trim paired-end RNA reads using Cutadapt.
    (If you have specific adapters for mRNA, define them in config.yaml
     and reference here.)
    """
    input:
        R1=lambda wc: config["samples"][wc.sample]["R1"],
        R2=lambda wc: config["samples"][wc.sample]["R2"]
    output:
        R1_trim = temp(f"{RESULTS_DIR}/cutadapt_rna/{{sample}}_R1.fastq"),
        R2_trim = temp(f"{RESULTS_DIR}/cutadapt_rna/{{sample}}_R2.fastq")
    log:
        f"{LOGS_DIR}/cutadapt_rna/{{sample}}.log"
    threads: 1
    shell:
        """
        # Example: if you have known 5'/3' adapters, specify them with -g/-a
        # For now, we leave them blank or you can add from config if needed.
        cutadapt \
            -o {output.R1_trim} \
            -p {output.R2_trim} \
            {input.R1} {input.R2} \
        2> {log}
        """

################################################################################
# 6) STAR Alignment for sRNA (Single-end)
################################################################################
rule star_align_srna:
    """
    Align single-end sRNA reads with STAR using sRNA-specific parameters.
    """
    input:
        f"{RESULTS_DIR}/cutadapt_srna/{{sample}}.fastq"
    output:
        bam       = f"{RESULTS_DIR}/STAR_srna_black6/{{sample}}/Aligned.sortedByCoord.out.bam",
        log_final = f"{RESULTS_DIR}/STAR_srna_black6/{{sample}}/Log.final.out"
    log:
        f"{LOGS_DIR}/STAR_srna_black6/{{sample}}.log"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS_DIR}/STAR_srna_black6/{wildcards.sample}
        STAR \
            --genomeDir {STAR_INDEX} \
            --readFilesIn {input} \
            {STAR_PARAMS_SRNA} \
            --outFileNamePrefix {RESULTS_DIR}/STAR_srna_black6/{wildcards.sample}/ \
        1> {log} 2>&1
        """

################################################################################
# 7) STAR Alignment for RNA (Paired-end)
################################################################################
rule star_align_rna:
    """
    Align paired-end RNA reads with STAR using default parameters.
    """
    input:
        R1 = f"{RESULTS_DIR}/cutadapt_rna/{{sample}}_R1.fastq",
        R2 = f"{RESULTS_DIR}/cutadapt_rna/{{sample}}_R2.fastq"
    output:
        bam       = f"{RESULTS_DIR}/STAR_rna_black6/{{sample}}/Aligned.sortedByCoord.out.bam",
        log_final = f"{RESULTS_DIR}/STAR_rna_black6/{{sample}}/Log.final.out"
    log:
        f"{LOGS_DIR}/STAR_rna_black6/{{sample}}.log"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS_DIR}/STAR_rna_black6/{wildcards.sample}
        STAR \
            --genomeDir {STAR_INDEX} \
            --readFilesIn {input.R1} {input.R2} \
            {STAR_PARAMS_RNA} \
            --outFileNamePrefix {RESULTS_DIR}/STAR_rna_black6/{wildcards.sample}/ \
        1> {log} 2>&1
        """



###############################################################################
# (Optional) RULE: Trinity assembly
#
# You might only do Trinity for RNA samples, possibly merging replicates.
# Below is a simplistic example for a single sample. 
# Often you will combine all samples or a subset to build a single assembly.
###############################################################################

rule trinity_assembly:
    input:
        r1 = "results/C57BL_6/trimmed/{sample}_R1_val.fq.gz",
        r2 = "results/C57BL_6/trimmed/{sample}_R2_val.fq.gz"
    output:
        f = "results/C57BL_6/trinity_black6/trinity-{sample}.Trinity.fasta"
    threads: 12
    log:
        stdout = "log/trinity_black6/{sample}.out",
        stderr = "log/trinity_black6/{sample}.err"
    resources:
        mem = 128000
    params:
        # Working directory for Trinity
        d = "results/C57BL_6/trinity_black6/trinity-{sample}"
    conda: "envs/trinity.yaml" 
    shell: r"""
        echo "********************************" | tee {log.stdout} {log.stderr}
        echo "NEW RUN: $(date)" | tee -a {log.stdout} {log.stderr}
        echo "********************************" | tee -a {log.stdout} {log.stderr}

        # Run Trinity
        Trinity \
            --normalize_reads \
            --min_kmer_cov 2 \
            --seqType fq \
            --left {input.r1} \
            --right {input.r2} \
            --output {params.d} \
            --SS_lib_type RF \
            --CPU {threads} \
            --max_memory 90G \
            >> {log.stdout} 2>> {log.stderr}

        echo "*******MOVE FASTA FILE**********" | tee -a {log.stdout} {log.stderr}
        mv {params.d}/Trinity.fasta {output.f}

        echo "*******CREATE EMPTY DIR*********" | tee -a {log.stdout} {log.stderr}
        mkdir {params.d}_empty >> {log.stdout} 2>> {log.stderr}

        echo "*******DELETE TRINITY FILES*********" | tee -a {log.stdout} {log.stderr}
        rsync -a --delete {params.d}_empty/ {params.d}/ >> {log.stdout} 2>> {log.stderr}

        echo "*******DELETE EMPTY DIRS*********" | tee -a {log.stdout} {log.stderr}
        rmdir {params.d} >> {log.stdout} 2>> {log.stderr}
        rmdir {params.d}_empty >> {log.stdout} 2>> {log.stderr}
    """