configfile: "config/config.json"

#tps=["16.5dpc","12.5dpp","20.5dpp"]
#reps=[1,2,3]
#strains=[
#    "129S1_SvImJ",
#    "AKR_J",
#    "A_J",
#    "BALB_cJ",
#    "C3H_HeJ",
#    "C57BL_6NJ",
#    "CAST_EiJ",
#    "CBA_J",
#    "DBA_2J",
#    "FVB_NJ",
#    "LP_J",
#    "NOD_ShiLtJ",
#    "NZO_HlLtJ",
#    "PWK_PhJ",
#    "SPRET_EiJ",
#    "WSB_EiJ"
#]
tps=["16.5dpc","20.5dpp"]
reps=[1]
strains=[
    "SPRET_EiJ",
]
rule all:
    input:
        #config["genomes"]["contaminants"]["indexes"]["STAR"] + "/genomeParameters.txt",
        config["genomes"][config["ref_genome"]]["indexes"]["STAR"]+ "/genomeParameters.txt",
        [f"results/filter_precursors/{sid}-{tp}.{rep}.{config['min_prec_len']}.csv" for sid in strains for tp in tps for rep in reps],
        #f"results/filter_precursors/C3H_HeJ-12.5dpp.1.500.csv",
        #f"results/filter_precursors/CBA_J-12.5dpp.1.500.csv",
        #[f"results/stringtie_tp/{sid}-{tp}/{sid}-{tp}.gtf" for sid in strains for tp in tps],
        [f"results/STAR_srna/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" for sid in strains for tp in tps for rep in reps],
        [f"results/peterplots/{sid}-{tp}.{rep}.{rawali}.svg" for sid in strains for tp in tps for rep in reps for rawali in ["aligned","raw"]]
        
rule get_genome:
    output:
       "results/genomes/ref_genome.fasta"
    resources:
        mem = 4000
    threads: 1
    params:
        url = config["genomes"][config["ref_genome"]]["url"]
    shell:"""
        wget -O {output}.gz {params.url}
        gunzip {output}.gz
    """

rule get_ann:
    output:
       "results/genomes/ref_genome.gtf"
    resources:
        mem = 4000
    threads: 1
    params:
        url = config["genomes"][config["ref_genome"]]["url_ann"]
    shell:"""
        wget -O {output}.gz {params.url}
        gunzip {output}.gz
    """

idx_cmd = "STAR --limitGenomeGenerateRAM  96000000000 --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.outdir} --genomeFastaFiles {input} > {log.stdout} 2> {log.stderr}"

rule STAR_idx_genome:
    shadow:"shallow"
    input:
        "results/genomes/ref_genome.fasta"
    output:
        f=config["genomes"][config["ref_genome"]]["indexes"]["STAR"] + "/genomeParameters.txt"
    threads: 8
    params:
        outdir=config["genomes"][config["ref_genome"]]["indexes"]["STAR"]
    resources:
        mem = 64000
    conda: "envs/star.yaml"
    log:
        stdout="log/STAR_idx/genome.out",
        stderr="log/STAR_idx/genome.err"
    shell: idx_cmd

rule STAR_idx_contaminants:
    shadow:"shallow"
    input:
        "res/contaminants/contaminants.fasta"
    output:
        f=config["genomes"]["contaminants"]["indexes"]["STAR"]+ "/genomeParameters.txt"
    threads: 8
    params:
        outdir=config["genomes"]["contaminants"]["indexes"]["STAR"]
    resources:
        mem = 16000
    conda: "envs/star.yaml"
    log:
        stdout="log/STAR_idx/contaminants.out",
        stderr="log/STAR_idx/contaminants.err"
    shell: idx_cmd

rule cutadapt_pe:
    input:
        r1=lambda wc: "{}/{}".format(config["data_dirs"]["rna"],config["samples"]["rna"][wc.sid][wc.tp][wc.rep]["1"]["r1"]),
        r2=lambda wc: "{}/{}".format(config["data_dirs"]["rna"],config["samples"]["rna"][wc.sid][wc.tp][wc.rep]["1"]["r2"])
    output:
        r1=temp("results/cutadapt_pe/{sid}-{tp}.{rep}_r1.fastq.gz"),
        r2=temp("results/cutadapt_pe/{sid}-{tp}.{rep}_r2.fastq.gz")
    log:
        stdout="log/cutadapt_pe/{sid}-{tp}.{rep}.out",
        stderr="log/cutadapt_pe/{sid}-{tp}.{rep}.err"
    threads: 3
    conda: "envs/cutadapt.yaml"
    resources:
        mem = 4000
    shell:
        "cutadapt {config[params][rna][cutadapt_pe]} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.stdout} 2> {log.stderr}"

rule cutadapt:
    input:
        lambda wc: ["{d}/{f}".format(d=config["data_dirs"]["srna"],f=f) for f in list(config["samples"]["srna"][wc.sid][wc.tp][wc.rep].values())]
    output:
        "results/cutadapt/{sid}-{tp}.{rep}.fastq.gz"
    log:
        stdout="log/cutadapt/{sid}-{tp}.{rep}.out",
        stderr="log/cutadapt/{sid}-{tp}.{rep}.err"
    threads: 5
    conda: "envs/cutadapt.yaml"
    resources:
        mem = 4000
    shell:
        "cat {input} | gunzip -c | cutadapt {config[params][srna][cutadapt]} -j {threads} -o {output} - > {log.stdout} 2> {log.stderr}"

rule STAR_rna:
    input:
        r1 = "results/cutadapt_pe/{sid}-{tp}.{rep}_r1.fastq.gz",
        r2 = "results/cutadapt_pe/{sid}-{tp}.{rep}_r2.fastq.gz",
        idx = config["genomes"][config["ref_genome"]]["indexes"]["STAR"] + "/genomeParameters.txt"
    output:
        bam="results/STAR_rna/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        log="results/STAR_rna/{sid}-{tp}.{rep}/Log.final.out"
    threads: 6
    log:
        stdout="log/STAR_rna/{sid}-{tp}.{rep}.out",
        stderr="log/STAR_rna/{sid}-{tp}.{rep}.err"
    params:
        outdir=lambda wc: f"results/STAR_rna/{wc.sid}-{wc.tp}.{wc.rep}",
        genome_idx = config["genomes"][config["ref_genome"]]["indexes"]["STAR"]
    conda: "envs/star.yaml"
    resources:
        mem = 64000
    shell: """
        mkdir -p {params.outdir}
        STAR --genomeDir {params.genome_idx} --outFileNamePrefix {params.outdir}/ --readFilesIn {input.r1} {input.r2} --runThreadN {threads} {config[params][rna][STAR]} > {log.stdout} 2> {log.stderr}
    """

rule STAR_idx__strain_wise_genome:
    shadow:"shallow"
    input:
        "results/genomes/ref_genome.fasta"
    output:
        f=config["genomes"][config["ref_genome"]]["indexes"]["STAR"] + "/genomeParameters.txt"
    threads: 8
    params:
        outdir=config["genomes"][config["ref_genome"]]["indexes"]["STAR"]
    resources:
        mem = 64000
    conda: "envs/star.yaml"
    log:
        stdout="log/STAR_idx/genome.out",
        stderr="log/STAR_idx/genome.err"
    shell: idx_cmd

rule STAR_rna_strain_wise:
    input:
        r1 = "results/cutadapt_pe/{sid}-{tp}.{rep}_r1.fastq.gz",
        r2 = "results/cutadapt_pe/{sid}-{tp}.{rep}_r2.fastq.gz",
        idx = config["genomes"][config["ref_genome"]]["indexes"]["STAR"] + "/genomeParameters.txt"
    output:
        bam="results/STAR_rna/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        log="results/STAR_rna/{sid}-{tp}.{rep}/Log.final.out"
    threads: 6
    log:
        stdout="log/STAR_rna/{sid}-{tp}.{rep}.out",
        stderr="log/STAR_rna/{sid}-{tp}.{rep}.err"
    params:
        outdir=lambda wc: f"results/STAR_rna/{wc.sid}-{wc.tp}.{wc.rep}",
        genome_idx = config["genomes"][config["ref_genome"]]["indexes"]["STAR"]
    conda: "envs/star.yaml"
    resources:
        mem = 64000
    shell: """
        mkdir -p {params.outdir}
        STAR --genomeDir {params.genome_idx} --outFileNamePrefix {params.outdir}/ --readFilesIn {input.r1} {input.r2} --runThreadN {threads} {config[params][rna][STAR]} > {log.stdout} 2> {log.stderr}
    """

rule downsample:
    input:
        ["results/cutadapt_pe/{{sid}}-{{tp}}.{rep}_{{r}}.fastq.gz".format(rep=rep) for rep in reps]
    output:
        r_tmp = temp("results/downsample/tmp_{sid}-{tp}_{r}.fastq.gz"),
        r = temp("results/downsample/{sid}-{tp}_{r}.fastq.gz")
    threads: 2
    log:
        stdout="log/downsample/{sid}-{tp}_{r}.out",
        stderr="log/downsample/{sid}-{tp}_{r}.err"
    resources:
        mem = 64000
    conda: "envs/seqtk.yaml" 
    shell:"""
        cat {input} > {output.r_tmp}
        seqtk sample -s123 {output.r_tmp} 50000000 | gzip > {output.r}
    """

rule trinity:
    input:
        r1 = "results/downsample/{sid}-{tp}_r1.fastq.gz",
        r2 = "results/downsample/{sid}-{tp}_r2.fastq.gz"
    output:
        f = "results/trinity/trinity-{sid}-{tp}.Trinity.fasta"
    threads: 12
    log:
        stdout="log/trinity/{sid}-{tp}.out",
        stderr="log/trinity/{sid}-{tp}.err"
    resources:
        mem = 128000
    params:
        d = "results/trinity/trinity-{sid}-{tp}"
    conda: "envs/trinity.yaml" 
    shell:"""
        echo "********************************" | tee {log.stdout} {log.stderr}
        echo "NEW RUN: $(date)" | tee -a {log.stdout} {log.stderr}
        echo "********************************" | tee -a {log.stdout} {log.stderr}
        Trinity --normalize_reads --min_kmer_cov 2 --seqType fq --left {input.r1} --right {input.r2} --output {params.d} --SS_lib_type RF --CPU {threads} --max_memory 90G >> {log.stdout} 2>> {log.stderr}
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

rule filter_trinity:
    input:
        f = "results/trinity/trinity-{sid}-{tp}.Trinity.fasta"
    output:
        f = "results/trinity/trinity-{sid}-{tp}.Trinity.{thres}.fasta"
    threads: 1
    log:
        stderr="log/trinity/{sid}-{tp}.{thres}.err"
    resources:
        mem = 4000
    conda: "envs/seqkit.yaml"
    shell:"""
        seqkit seq -m 500 {input.f} > {output.f} 2> {log.stderr}
    """

rule STAR_idx_trinity:
    shadow:"shallow"
    input:
        "results/trinity/trinity-{sid}-{tp}.Trinity.{thres}.fasta"
    output:
        d=directory("results/trinity_idx/{sid}-{tp}.{thres}")
    threads: 8
    params:
        outdir="results/trinity_idx/{sid}-{tp}.{thres}"
    resources:
        mem = 96000
    conda: "envs/star.yaml"
    log:
        stdout="log/trinity_idx/{sid}-{tp}.{thres}.out",
        stderr="log/trinity_idx/{sid}-{tp}.{thres}.err"
    shell: idx_cmd

rule STAR_trinity:
    shadow:"shallow"
    input:
        sample = "results/cutadapt/{sid}-{tp}.{rep}.fastq.gz",
        genome_idx="results/trinity_idx/{sid}-{tp}.{thres}"
    output:
        bam="results/STAR_trinity/{sid}-{tp}.{rep}.{thres}/Aligned.sortedByCoord.out.bam",
        log="results/STAR_trinity/{sid}-{tp}.{rep}.{thres}/Log.final.out"
    threads: 6
    log:
        stdout="log/STAR_trinity/{sid}-{tp}.{rep}.{thres}.out",
        stderr="log/STAR_trinity/{sid}-{tp}.{rep}.{thres}.err"
    resources:
        mem = 96000
    conda: "envs/star.yaml"
    params:
        outdir=lambda wc: f"results/STAR_trinity/{wc.sid}-{wc.tp}.{wc.rep}.{wc.thres}"
    shell: """
        mkdir -p {params.outdir}
        STAR --genomeChrBinNbits 9 --genomeDir {input.genome_idx} --outFileNamePrefix {params.outdir}/ --readFilesIn {input.sample} --runThreadN {threads} {config[params][srna][STAR]} > {log.stdout} 2> {log.stderr}
    """

rule trinity_idx_bam:
    input:
        bam="results/STAR_trinity/{sid}-{tp}.{rep}.{thres}/Aligned.sortedByCoord.out.bam",
    output:
        bam="results/STAR_trinity/{sid}-{tp}.{rep}.{thres}/Aligned.sortedByCoord.out.bam.bai",
    threads: 6
    log:
        stdout="log/trinity_idx_bam/{sid}-{tp}.{rep}.{thres}.out",
        stderr="log/trinity_idx_bam/{sid}-{tp}.{rep}.{thres}.err"
    resources:
        mem = 4000
    conda: "envs/samtools.yaml"
    shell: """
        samtools index {input.bam} > {log.stdout} 2> {log.stderr}
    """

rule filter_precursors:
    input:
        bam="results/STAR_trinity/{sid}-{tp}.{rep}.{thres}/Aligned.sortedByCoord.out.bam",
        bai="results/STAR_trinity/{sid}-{tp}.{rep}.{thres}/Aligned.sortedByCoord.out.bam.bai",
        log="results/STAR_trinity/{sid}-{tp}.{rep}.{thres}/Log.final.out"
    output:
        "results/filter_precursors/{sid}-{tp}.{rep}.{thres}.csv",
    threads: 1
    log:
        stdout="log/filter_precursors/{sid}-{tp}.{rep}.{thres}.out",
        stderr="log/filter_precursors/{sid}-{tp}.{rep}.{thres}.err"
    resources:
        mem = 4000
    conda: "envs/filter_precursors.yaml"
    script: "scripts/filter_precursors.py"

rule merge_bams:
    input:
        bams = [f"results/STAR_rna/{{sid}}-{{tp}}.{rep}/Aligned.sortedByCoord.out.bam" for rep in reps]
    output:
        bam="results/merge_bams/{sid}-{tp}.bam",
        header=temp("results/merge_bams/{sid}-{tp}.header")
    threads: 3
    log:
        stdout="log/merge_bams/{sid}-{tp}.out",
        stderr="log/merge_bams/{sid}-{tp}.err"
    conda: "envs/samtools.yaml"
    resources:
        mem = 16000
    shell: """
        samtools view -H {input.bams[0]} > {output.header} 2> {log.stderr} ; samtools cat -h {output.header} {input.bams} 2>> {log.stderr} | samtools view -h - 2>> {log.stderr} | awk -v strType=2 -f workflow/scripts/tagXSstrandedData.awk 2> {log.stderr} | samtools sort -O BAM -o {output.bam} 2> {log.stderr}
    """

rule stringtie:
    input:
        bam="results/merge_bams/{sid}-{tp}.bam",
        ann="results/genomes/ref_genome.gtf"
    output:
        gtf="results/stringtie_tp/{sid}-{tp}/{sid}-{tp}.gtf",
        abundance="results/stringtie_tp/{sid}-{tp}/{sid}-{tp}.abundance.tab"
    threads: 8
    log:
        stdout="log/stringtie_tp/{sid}-{tp}.out",
        stderr="log/stringtie_tp/{sid}-{tp}.err"
    conda: "envs/stringtie.yaml"
    resources:
        mem = 8000
    shell: """
        stringtie {input.bam} -o {output.gtf} -A {output.abundance} -p {threads} -G {input.ann} --rf -l {wildcards.sid}-{wildcards.tp} -B > {log.stdout} 2> {log.stderr}
    """

rule STAR_srna:
    input:
        sample = "results/cutadapt/{sid}-{tp}.{rep}.fastq.gz",
        idx = config["genomes"][config["ref_genome"]]["indexes"]["STAR"]+ "/genomeParameters.txt",
    output:
        bam="results/STAR_srna/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        log="results/STAR_srna/{sid}-{tp}.{rep}/Log.final.out"
    threads: 6
    log:
        stdout="log/STAR_srna/{sid}-{tp}.{rep}.out",
        stderr="log/STAR_srna/{sid}-{tp}.{rep}.err"
    resources:
        mem = 96000
    conda: "envs/star.yaml"
    params:
        outdir=lambda wc: f"results/STAR_srna/{wc.sid}-{wc.tp}.{wc.rep}",
        genome_idx = config["genomes"][config["ref_genome"]]["indexes"]["STAR"]
    shell: """
        mkdir -p {params.outdir}
        STAR --genomeDir {params.genome_idx} --outFileNamePrefix {params.outdir}/ --readFilesIn {input.sample} --runThreadN {threads} {config[params][srna][STAR]} > {log.stdout} 2> {log.stderr}
    """

rule STAR_srna_idx:
    input:
        bam="results/STAR_srna/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
    output:
        bai="results/STAR_srna/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/STAR_srna_idx/{sid}-{tp}.{rep}.log"
    threads: 1
    conda: "envs/samtools.yaml"
    resources:
        mem = 4000
    shell:
        "samtools index {input} > {log} 2>&1"

rule collapse_raw:
    input:
        "results/cutadapt/{sid}-{tp}.{rep}.fastq.gz",
    output:
        "results/collapse/{sid}-{tp}.{rep}.raw.fasta.gz"
    log:
        stdout="log/collapse/{sid}-{tp}.{rep}.raw.out",
        stderr="log/collapse/{sid}-{tp}.{rep}.raw.err"
    conda: "envs/tstk.yaml"
    threads: 1
    resources:
        mem = 64000 
    shell: """ 
        collapse {input} {output} > {log.stdout} 2> {log.stderr}
    """

rule collapse_aligned:
    input:
        bam="results/STAR_srna/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"
    output:
        "results/collapse/{sid}-{tp}.{rep}.aligned.fasta.gz"
    log:
        stdout="log/collapse/{sid}-{tp}.{rep}.aligned.out",
        stderr="log/collapse/{sid}-{tp}.{rep}.aligned.err"
    conda: "envs/tstk.yaml"
    threads: 1
    resources:
        mem = 96000 
    shell: """ 
        collapse {input.bam} {output} > {log.stdout} 2> {log.stderr}
    """

rule peterplots:
    input:
        "results/collapse/{sid}-{tp}.{rep}.{rawali}.fasta.gz"
    #params:
    #    libsize = lambda wc: config["libsizes"][wc.sample]
    output:
        normal="results/peterplots/{sid}-{tp}.{rep}.{rawali}.svg",
        collapsed="results/peterplots/{sid}-{tp}.{rep}.collapsed.{rawali}.svg"
    log:
        stdout="log/peterplots/{sid}-{tp}.{rep}.{rawali}.out",
        stderr="log/peterplots/{sid}-{tp}.{rep}.{rawali}.err"
    conda: "envs/tstk.yaml"
    threads: 1
    resources:
        mem = 16000
    shell: """ 
        peterplot {input} {output.normal} --uncollapse --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 > {log.stdout} 2> {log.stderr}
        peterplot {input} {output.collapsed} --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 >> {log.stdout} 2>> {log.stderr}
    """
