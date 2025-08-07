

#####

#RNA seq Analysis (STAR alingment)

####

rule rna_cutadapt_pe:
    localrule: True
    input: 
        [lambda wc: "{}/{}".format(config["data_dirs"]["rna"],config["samples"]["rna"][wc.sid][wc.tp][wc.rep]["1"]["r1"]),
        lambda wc: "{}/{}".format(config["data_dirs"]["rna"],config["samples"]["rna"][wc.sid][wc.tp][wc.rep]["1"]["r2"])]
    output:
        fastq1= temp(config['resultDir']  + "/cutadapt_pe_rna/{sid}-{tp}.{rep}_r1.fastq"),
        fastq2= temp(config['resultDir']  + "/cutadapt_pe_rna/{sid}-{tp}.{rep}_r2.fastq"),
        qc = config['resultDir']  + "/cutadapt_pe_rna/{sid}-{tp}.{rep}.qc.txt"
    params:
        adapters="",
        extra=config['params']['rna']['cutadapt_pe']
    log:
        config['logDir']  + "/cutadapt_pe_rna/{sid}/{sid}-{tp}.{rep}.log"
    threads: 1
    priority:1
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    wrapper:
        "v5.8.2/bio/cutadapt/pe"




rule STAR_alingment_rna_strain_wise:
    localrule: True
    input:
        fq1 = [config['resultDir']  + "/cutadapt_pe_rna/{sid}-{tp}.{rep}_r1.fastq"],
        fq2 = [config['resultDir']  + "/cutadapt_pe_rna/{sid}-{tp}.{rep}_r2.fastq"],
        idx = config['resultDir'] + "/indexs/{sid}",
    output:
        aln=config['resultDir']  + "/STAR_rna_strain_wise/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        log_final=config['resultDir']  + "/STAR_rna_strain_wise/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    resources:
        mem = 400000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    priority:5
    log:
        config['logDir']  + "/STAR_rna_strain_wise/{sid}/{sid}-{tp}.{rep}.log"
    params:
        idx = config['resultDir'] + "/indexs/{sid}",
        extra = config['params']['rna']['STAR']
    wrapper:
        "v5.8.2/bio/star/align"


