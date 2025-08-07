
"""

seq count from Collapsed sRNA raw sequncing data 

"""

rule collapse_raw_seq_count_srna:
    input:
        fasta = config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.raw.fasta.gz"
    output:
        temp(config['resultDir']  + "/seqCount/{sid}-{tp}.{rep}.raw.count.gz")
    log:
        stderr=config['logDir']  + "/seqCount/{sid}-{tp}.{rep}.raw.err"
    threads: 1
    priority:4
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    conda: config['envsDir'] + "/seqkit.yaml"
    shell:
        r"""
        seqkit fx2tab {input.fasta} | awk '{{split($1,a,"-"); print $2 "\t" a[1] "\t" a[2]}}' | sed  '1i seq \t seq_id \t {wildcards.sid}-{wildcards.tp}.{wildcards.rep}' |  gzip --stdout > {output} 2> {log.stderr}
        """



"""

collapsed alinged BAM files using Tom's Scripts tstk package

"""

rule strains_collapse_aligned_bam_count_srna:
    localrule: True
    input:
        bam = config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.aligned.fasta.gz"
    output:
        config['resultDir']  + "/seqCount/{sid}-{tp}.{rep}.aligned.count.gz"
    log:
        stderr=config['logDir']  + "/SeqCount/{sid}-{tp}.{rep}.aligned.err"
    threads: 1
    priority:4
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    conda: config['envsDir'] + "/seqkit.yaml"
    shell:
        r"""
        seqkit fx2tab {input.bam} | awk '{{split($1,a,"-"); print $2 "\t" a[1] "\t" a[2]}}' | sed  '1i seq \t seq_id \t {wildcards.sid}-{wildcards.tp}.{wildcards.rep}' |  gzip --stdout > {output} 2> {log.stderr}
        """

