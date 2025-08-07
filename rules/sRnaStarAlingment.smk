"""

alinging sRNA seq data to pacbio assemblies from thomas keen lab 

"""
rule cutadapt_input:
    input:
        lambda wc: [
            "{}/{}".format(config["data_dirs"]["srna"], config["samples"]["srna"][wc.sid][wc.tp][wc.rep][key])
            for key in config["samples"]["srna"][wc.sid][wc.tp][wc.rep]
        ]
    output:
        temp(config['resultDir']  + "/cutadapt_copy/{sid}-{tp}.{rep}.fastq"),
    log:
        stderr=config['logDir']  + "/cutadapt_copy/{sid}-{tp}.{rep}_cutadapt.err"
    priority: 1
    resources:
        mem=200000,
        queue='normal',
        tmpdir=config['tmpDir']
    shell:
        "zcat {input} > {output} 2> {log.stderr}"


rule cutadapt_input_wrong:
    input:
        lambda wc: [
            f"{config['data_dirs']['srna']}/{f}"
            for f in list(config["samples"]["srna"][wc.sid][wc.tp][wc.rep].values())
        ]
    output:
        temp(f"{config['resultDir']}/cutadapt_se/{{sid}}-{{tp}}.{{rep}}.fastq_qe")
    log:
        stderr=f"{config['logDir']}/cutadapt_se/{{sid}}-{{tp}}.{{rep}}_cutadapt.err"
    priority: 1
    resources:
        mem=200000,
        queue='normal',
        tmpdir=config['tmpDir']
    shell:
        "zcat {input} > {output}  2> {log.stderr}"


rule cutadapt_srna:
    input:
        config['resultDir'] + "/cutadapt_copy/{sid}-{tp}.{rep}.fastq"
    output:
        fastq = temp(config['resultDir'] + "/cutadapt_se/{sid}-{tp}.{rep}.fastq"),
        qc = config['resultDir'] + "/cutadapt_se/{sid}-{tp}.{rep}.txt"
    log:
        config['logDir'] + "/cutadapt_se/{sid}-{tp}.{rep}.err"
    threads: 1
    params:
        adapters = "",
        extra = config['params']['srna']['cutadapt'] 
    priority: 2
    resources:
        mem = 200000,
        queue = 'normal',
        tmpdir = config['tmpDir']
    wrapper:
        "v1.10.0/bio/cutadapt/se"



rule STAR_alingment_srna_strain:
    input:
        fq1 =  config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq",
        idx = config['resultDir'] + "/indexs/{sid}",
        #fq1 = config['resultDir']  + "/STAR_srna_strain_wise/tmp_fasta/{sid}-{tp}.{rep}.fasta"
    output:
        aln=temp(config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"),
        log_final=config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    priority:3
    log:
        config['logDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}.log"
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        idx = config['resultDir'] + "/indexs/{sid}",
        extra = config['params']['srna']['STAR'] 
    wrapper:
        "v5.8.2/bio/star/align"


rule STAR_alingment_srna_strain_old_genomes:
    input:
        fq1 =  config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq",
        idx = config['resources'] + '/old_genomes/index/{sid}' ,
        #fq1 = config['resultDir']  + "/STAR_srna_strain_wise/tmp_fasta/{sid}-{tp}.{rep}.fasta"
    output:
        aln=temp(config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"),
        log_final=config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    priority:3
    log:
        config['logDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}.log",
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        idx = config['resultDir'] + "/indexs/{sid}",
        extra = config['params']['srna']['STAR']
    wrapper:
        "v5.8.2/bio/star/align"


"""

alinging sRNA seq data to refrance genome 
"""

rule STAR_alingment_srna_GRCm39_106:
    input:
        fq1 = config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq",
        idx = config['resultDir'] + "/indexs/GRCm39.106",
        #fq1 = config['resultDir']  + "/STAR_srna_strain_wise/tmp_fasta/{sid}-{tp}.{rep}.fasta"
    output:
        bam=temp(config['resultDir']  + "/STAR_srna_GRCm39_106/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"),
        log_final=config['resultDir']  + "/STAR_srna_GRCm39_106/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    priority:3
    log:
        config['logDir']  + "/STAR_srna_GRCm39_106/{sid}-{tp}.{rep}.log",
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        idx = config['resultDir'] + "/indexs/GRCm39.106",
        extra = config['params']['srna']['STAR']
    wrapper:
        "v5.8.2/bio/star/align"

rule STAR_alingment_srna_GRCm38_68:
    input:
        fq1 = config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq",
        idx = config['resultDir']  + "/indexs/GRCm38.68",
        #fq1 = config['resultDir']  + "/STAR_srna_strain_wise/tmp_fasta/{sid}-{tp}.{rep}.fasta"
    output:
        bam=temp(config['resultDir']  + "/STAR_srna_GRCm38_68/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"),
        log_final=config['resultDir']  + "/STAR_srna_GRCm38_68/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    priority:3
    log:
        config['logDir']  + "/STAR_srna_GRCm38_68/{sid}-{tp}.{rep}.log",
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        idx = config['resultDir'] + "/indexs/GRCm38_68",
        extra = config['params']['srna']['STAR']
    wrapper:
        "v5.8.2/bio/star/align"
