

"""
betools intersect of zamore pirna and and filtred trinity precursors for GRCm38

"""

checkpoint bedtools_intersection_zamore_ref_GRCm38_68:
    input:
        left = config['resultDir']  + "/zamore/GRCm38.68/piRNA_gene_annotation_GRCm38.68_without_chr.bed",
        right= config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bed12.bed.gz",
    output:
        temp(config['resultDir']  + "/bedIntersect/GRCm38.68/{sid}_{tp}_intersected.bed")
    params:
        ## Add optional parameters
        extra="-wo -s" 
    resources:
        mem = 4000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    priority:6
    log:
        config['logDir']  + "/bedtool_intersect/{sid}_{tp}_GRCm38.68.log"
    wrapper:
        "v1.7.0/bio/bedtools/intersect"


"""
bzip bed intersect files files (refrnace genome )
"""

rule bgzip_bedtools_intersection_zamore_ref_GRCm38_68:
    input:
        config['resultDir']  + "/bedIntersect/GRCm38.68/{sid}_{tp}_intersected.bed"
    output:
        config['resultDir']  + "/bedIntersect/GRCm38.68/{sid}_{tp}_intersected.bed.gz",
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 200000,
        tmpdir= config['tmpDir']
    priority:8
    log:
        config['logDir']  + "/bgzip_bedtools_intersection_zamore_ref_genome/GRCm38.68/{sid}_{tp}_intersected.bed.log",
    wrapper:
        "v1.7.0/bio/bgzip"







"""
betools intersect of zamore pirna and and filtred trinity precursors for GRCm39

"""

checkpoint bedtools_intersection_zamore_ref_GRCm39_106:
    input:
        left = config['resultDir']  + "/zamore/GRCm39.106/piRNA_gene_annotation_GRCm39.106_without_chr.bed",
        right= config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.bed12.bed.gz",
    output:
        temp(config['resultDir']  + "/bedIntersect/GRCm39.106/{sid}_{tp}_intersected.bed")
    params:
        ## Add optional parameters
        extra="-wo -s" 
    resources:
        mem = 4000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    priority:7
    log:
        config['logDir']  + "/bedtool_intersect/{sid}_{tp}_GRCm39.106.log"
    wrapper:
        "v1.7.0/bio/bedtools/intersect"


"""
bzip bed intersect files files (refrnace genome )
"""

rule bgzip_bedtools_intersection_zamore_ref_GRCm39_106:
    input:
        config['resultDir']  + "/bedIntersect/GRCm39.106/{sid}_{tp}_intersected.bed"
    output:
        config['resultDir']  + "/bedIntersect/GRCm39.106/{sid}_{tp}_intersected.bed.gz",
    params:
        extra="", # optional
    threads: 1
    priority:9
    resources:
        queue  = 'normal',
        mem = 200000,
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bgzip_bedtools_intersection_zamore_ref_genome/GRCm39.106/{sid}_{tp}_intersected.bed.log",
    wrapper:
        "v1.7.0/bio/bgzip"









rule bedtools_intersection_zamore_strains:
    localrule: True
    input:
        left = config['resultDir']  + "/zamore/with_05/{sid}/piRNA_gene_annotation_{sid}_without_chr.bed",
        right= config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed.gz",
    output:
        config['resultDir']  + "/zamoreIntersection/strains/{sid}/{sid}_{tp}_intersected.bed"
    params:
        ## Add optional parameters
        extra="-wo -s" 
    resources:
        mem = 4000,
        queue  = 'normal',
    priority:8
    log:
        config['logDir']  + "/bedtool_intersect/{sid}_{tp}.log"
    wrapper:
        "v5.8.2/bio/bedtools/intersect"


"""
bzip bed intersect files (strains genome )
"""

rule bgzip_bedtools_intersection_zamore_strains:
    input:
        config['resultDir']  + "/bedIntersect/strains/{sid}/{sid}_{tp}_intersected.bed"
    output:
        config['resultDir']  + "/bedIntersect/strains/{sid}/{sid}_{tp}_intersected.bed.gz",
    params:
        extra="", # optional
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bgzip_bedtools_intersection_zamore_strains/{sid}/{sid}_{tp}_intersected.bed.log",
    wrapper:
        "v5.8.2/bio/bgzip"




#####

#intersection of strain based mapped bed files into one file (from diffrent timepoints)

#####



checkpoint bedtools_intersection_minimap_diffrent_timepoints:
    input:
        a = config['resultDir']  + "/minimap2/strains/{sid}/{sid}-12.5dpp.500.bed12.bed.gz",
        b = config['resultDir']  + "/minimap2/strains/{sid}/{sid}-16.5dpc.500.bed12.bed.gz",
        c = config['resultDir']  + "/minimap2/strains/{sid}/{sid}-20.5dpp.500.bed12.bed.gz",
        #refrance
        d = config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed.gz",
    output:
        temp(config['resultDir']  + "/strainBedIntersect/{sid}/{sid}_{tp}_intersected.tab")
    params:
        ## Add optional parameters
        extra="-wao -s -names {sid}-12.5dpp {sid}-16.5dpc {sid}-20.5dpp" 
    conda: config['envsDir'] + "/bedtools.yaml"
    priority:8
    resources:
        mem = 4000,
        queue  = 'normal',
    log:
        config['logDir']  + "/strain_bedtool_intersect/{sid}_{tp}.log"
    shell:
        "bedtools intersect -a {input.d} -b {input.a} {input.b} {input.c} {params.extra} > {output} 2> {log}"


