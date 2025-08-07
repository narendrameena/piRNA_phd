
rule filter_trinity:
    localrule: True
    input:
        f = config['resultDir'] + "/trinity/trinity-{sid}-{tp}.Trinity.fasta"
    output:
        f = temp(config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta")
    threads: 1
    log:
        stderr=config['logDir']  +"/trinity/{sid}-{tp}.500.err"
    resources:
        mem = 4000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    priority:6
    conda: config['envsDir'] + "/seqkit.yaml"
    shell:
        """
        seqkit seq -m 500 {input.f} > {output.f} 2> {log.stderr}
        """
rule STAR_index_trinity:
    localrule: True
    input:
        fasta=config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta"
    output:
        directory(config['resultDir'] + "/trinity/index/{sid}-{tp}.500")
    threads: 4
    params:
        extra="--limitGenomeGenerateRAM 96000000000"
    message:
        "**** Creating STAR index for {wildcards.sid}-{wildcards.tp} ****"
    resources:
        mem_mb=200000,
        queue='normal'
    log:
        config['logDir'] + "/trinity/index/{sid}-{tp}.log"
    wrapper:
        "v5.8.2/bio/star/index"

rule STAR_trinity:
    localrule: True
    input:
        fq1 =  config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq",
        idx = config['resultDir'] + "/trinity/index/{sid}-{tp}.500",
        #fq1 = config['resultDir']  + "/STAR_srna_strain_wise/tmp_fasta/{sid}-{tp}.{rep}.fasta"
    output:
        aln= config['resultDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        log_final =config['resultDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    priority:8
    log:
        config['logDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}.log",
    resources:
        mem = 400000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        idx = config['resultDir'] + "/trinity/index/{sid}-{tp}.500",
        extra = config['params']['srna']['STAR']
    wrapper:
        "v5.8.2/bio/star/align"
    



rule trinity_idx_bam:
    localrule: True
    input:
        config['resultDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
    output:
        config['resultDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
    log:
        config['logDir']  + "/trinity/star_index/{sid}-{tp}.{rep}.log"
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        extra="",
    threads: 4
    priority:9
    wrapper:
        "v5.8.2/bio/samtools/index"
    

rule filter_precursors:
    localrule: True
    input:
        bam=config['resultDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        bai=config['resultDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
        log=config['resultDir']  + "/trinity/star/{sid}/{sid}-{tp}.{rep}/Log.final.out",  
    output:
        protected(config['resultDir']  + "/trinity/filter_precursors/{sid}/{sid}-{tp}.{rep}.500.csv"),
    threads: 1
    log:
        stdout=config['logDir']  + "/trinity/filter_precursors/{sid}-{tp}.{rep}.500.out",
        stderr=config['logDir']  + "/trinity/filter_precursors/{sid}-{tp}.{rep}.500.err"
    priority:10
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    conda: 
        config['envsDir'] +  "/filter_precursors.yaml"
    script: 
        config['scriptDir'] + "/python/filter_precursors.py" 
    

rule filter_trinity_precursors_bed:
    localrule: True
    input:
        csv = config['resultDir']  + "/trinity/filter_precursors/{sid}/{sid}-{tp}.{rep}.500.csv",
        bed = config['resultDir'] + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed",
        script = config['scriptDir'] + "/sh/piRNA_filter_precursors.sh",  
    output:
        bed = protected(config['resultDir']  + "/filter_precursors_bed/{sid}/{sid}-{tp}.{rep}.bed"),
    threads: 1
    log:
        config['logDir']  + "/trinity/filter_precursors/{sid}-{tp}.{rep}.500.err"
    priority:10
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    conda: 
        config['envsDir'] +  "/filter_precursors.yaml"
    script: 
        """
        bash {input.script} {input.csv} {input.bed} {output.bed}
        """ 
    
