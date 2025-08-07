"""

lift annotation from mm11 ot pac_bio genomes 

"""


rule liftoff_strains:
    input:
        ref = config['resultDir']  + "/ref_genome/GRCm39.106.fasta",
        tgt = config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta",
        ann = config['resultDir']  + "/ref_genome/GRCm39.106.gtf", 
    output:
        main = config['resultDir']  + "/ref_genome/{sid}/{sid}.gff3",
        unmapped = config['resultDir']  + "/ref_genome/{sid}/{sid}.unmapped.txt",    
    resources:
        queue  = 'normal',
        mem = 200000
    threads: 4
    priority:1
    params:
        extra="-mm2_options='-ax asm20 --end-bonus 5 --eqx -N 50 -p 0.5' -exclude_partial -p 4 -dir tmp/{sid}"
    log:
        config['logDir']  + "/liftoff/{sid}.log",
    wrapper:
        "v1.7.0/bio/liftoff"

rule liftoff_trinity_precursors:
    input:
        ref = config['resultDir']  + "/ref_genome/GRCm39.106.fasta",
        tgt = config["tom_results"] + "/trinity/trinity-{sid}-{tp}.Trinity.500.fasta",
        ann = config['resultDir']  + "/ref_genome/GRCm39.106.gtf", 
    output:
        main = config['resultDir']  + "/trinity_genome/{sid}/trinity-{sid}-{tp}.Trinity.500.gff3",
        unmapped = config['resultDir']  + "/trinity_genome/{sid}/trinity-{sid}-{tp}.Trinity.500.unmapped.txt",
    threads: 4
    priority:1
    resources:
        queue  = 'normal',
        mem = 200000
    params:
        extra="-mm2_options='-ax asm20 --end-bonus 5 --eqx -N 50 -p 0.5'  -p 4 -dir tmp/{sid}-{tp}",
    log:
        config['logDir']  + "/trinity_genome/{sid}-{tp}.log",
    wrapper:
        "v5.8.2/bio/liftoff"







###        
#    shell:
#        "liftoff {params.extra} -g {input.ann} -o {output.main} -u {output.unmapped} {input.tgt} {input.ref} > {log}"
#conda: config['envsDir'] + "/liftoff.yaml"
