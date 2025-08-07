

"""
betools jaccard score of zamore pirna and and filtred trinity precursors for GRCm38

"""

checkpoint bedtools_jaccard_zamore_ref_GRCm38_68:
    input:
        left = config['resultDir']  + "/zamore/GRCm38.68/piRNA_gene_annotation_GRCm38.68_without_chr_sort_bed6.bed",
        right= config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.sort_bed6.bed",
        genome = config['resultDir']  + "/ref_genome/GRCm38.68.fasta.txt"
    output:
        temp(config['resultDir']  + "/bedJaccard/GRCm38.68/{sid}_{tp}_jaccard.tab")
    params:
        ## Add optional parameters
        extra=" -s" 
    priority:4
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 4000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bedtool_jaccard/GRCm38.68/{sid}_{tp}.log"
    shell:
        "bedtools jaccard {params.extra} -a {input.left} -b {input.right} -g {input.genome}> {output} 2> {log}"


"""
bzip bed intersect files files (refrnace genome )
"""

rule bgzip_bedtools_jaccard_zamore_ref_GRCm38_68:
    input:
        config['resultDir']  + "/bedJaccard/GRCm38.68/{sid}_{tp}_jaccard.tab"
    output:
        config['resultDir']  + "/bedJaccard/GRCm38.68/{sid}_{tp}_jaccard.tab.gz",
    params:
        extra="", # optional
    threads: 1
    priority:4
    resources:
        queue  = 'normal',
        mem = 200000,
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bgzip_bedtools_jaccard_zamore_ref_genome/GRCm38.68/{sid}_{tp}_GRCm38.68_jaccard.bed.log",
    wrapper:
        "v5.8.2/bio/bgzip"


"""
betools jaccard score of zamore pirna and and filtred trinity precursors for GRCm39.106

"""

checkpoint bedtools_jaccard_zamore_ref_GRCm39_106:
    input:
        left = config['resultDir']  + "/zamore/GRCm39.106/piRNA_gene_annotation_GRCm39.106_without_chr_sort_bed6.bed",
        right= config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.sort_bed6.bed",
        genome = config['resultDir']  + "/ref_genome/GRCm39.106.fasta.txt"
    output:
        temp(config['resultDir']  + "/bedJaccard/GRCm39.106/{sid}_{tp}_jaccard.tab")
    params:
        ## Add optional parameters
        extra=" -s" 
    priority:4
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 4000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bedtool_jaccard/GRCm39.106/{sid}_{tp}.log"
    shell:
        "bedtools jaccard {params.extra} -a {input.left} -b {input.right} -g {input.genome} > {output} 2> {log}"


"""
bzip bed intersect files files (refrnace genome )
"""

rule bgzip_bedtools_jaccard_zamore_ref_GRCm39_106:
    input:
        config['resultDir']  + "/bedJaccard/GRCm39.106/{sid}_{tp}_jaccard.tab"
    output:
        config['resultDir']  + "/bedJaccard/GRCm39.106/{sid}_{tp}_jaccard.tab.gz",
    params:
        extra="", # optional
    threads: 1
    priority:4
    resources:
        queue  = 'normal',
        mem = 200000,
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bgzip_bedtools_jaccard_zamore_ref_genome/GRCm39.106/{sid}_{tp}_GRCm39.106_jaccard.bed.log",
    wrapper:
        "v5.8.2/bio/bgzip"



checkpoint bedtools_jaccard_zamore_strains:
    localrule: True
    input:
        left = config['resultDir']  + "/zamore/{sid}/piRNA_gene_annotation_{sid}_without_chr_sort_bed6.bed",
        right = config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.sort_bed6.bed",
        genome = config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta.txt"
    output:
        temp(config['resultDir']  + "/bedJaccard/strains/{sid}/{sid}_{tp}_Jaccard.tab")
    params:
        ## Add optional parameters
        extra=" -s" 
    priority:4
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 4000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bedtool_Jaccard/strains/{sid}_{tp}_{sid}.log"
    shell:
        "bedtools jaccard {params.extra} -a {input.left} -b {input.right} -g {input.genome}> {output} 2> {log}"


"""
bzip bed intersect files (strains genome )
"""

rule bgzip_bedtools_jaccard_zamore_strains:
    localrule: True
    input:
        config['resultDir']  + "/bedJaccard/strains/{sid}/{sid}_{tp}_Jaccard.tab"
    output:
        config['resultDir']  + "/bedJaccard/strains/{sid}/{sid}_{tp}_Jaccard.tab.gz",
    params:
        extra="", # optional
    threads: 1
    priority:4
    resources:
        queue  = 'normal',
        mem = 200000,
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bgzip_bedtools_Jaccard_zamore_strains/{sid}/{sid}_{tp}_Jaccard.bed.log",
    wrapper:
        "v5.8.2/bio/bgzip"

