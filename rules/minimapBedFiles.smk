
####

#minimap Alingment pac bio assamblies 

####



"""

getting BED file from minimap BAM files 

"""

rule minimap_bam_to_bed_precursors_pacbio:
    localrule: True
    input:
        bam_file = config['resultDir'] + "/minimap2/strains/{sid}/{sid}-{tp}.500.bam"
    output:
        bed_file = config['resultDir'] + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed"
    threads: 1
    params:
        bedtools_opts = "-bed12 -i"
    conda:
        config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 200000,
        queue = 'normal'
    log:
        log_file = config['logDir'] + "/minimap2_bam_to_bed_strains/{sid}-{tp}.bed.log"
    shell:
        """
        bedtools bamtobed {params.bedtools_opts} {input.bam_file} | awk 'NF>=6' > {output.bed_file} 2> {log.log_file}
        """

rule minimap_bam_to_bed_precursors_6NJ:
    localrule: True
    input:
        bam_file = config['resultDir'] + "/minimap2/C57BL_6NJ/{sid}/{sid}-{tp}.500.bam"
    output:
        bed_file = config['resultDir'] + "/minimap2/C57BL_6NJ/{sid}/{sid}-{tp}.500.bed12.bed"
    threads: 1
    params:
        bedtools_opts = "-bed12 -i"
    conda:
        config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 200000,
        queue = 'normal'
    log:
        log_file = config['logDir'] + "/minimap2_bam_to_bed_C57BL_6NJ/{sid}-{tp}.bed.log"
    shell:
        """
        bedtools bamtobed {params.bedtools_opts} {input.bam_file} | awk 'NF>=6' > {output.bed_file} 2> {log.log_file}
        """

rule minimap_bam_to_bed_precursors_SPRET_EiJ:
    localrule: True
    input:
        bam_file = config['resultDir'] + "/minimap2/SPRET_EiJ/{sid}/{sid}-{tp}.500.bam"
    output:
        bed_file = config['resultDir'] + "/minimap2/SPRET_EiJ/{sid}/{sid}-{tp}.500.bed12.bed"
    threads: 1
    params:
        bedtools_opts = "-bed12 -i"
    conda:
        config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 200000,
        queue = 'normal'
    log:
        log_file = config['logDir'] + "/minimap2_bam_to_bed_SPRET_EiJ/{sid}-{tp}.bed.log"
    shell:
        """
        bedtools bamtobed {params.bedtools_opts} {input.bam_file} | awk 'NF>=6' > {output.bed_file} 2> {log.log_file}
        """


"""

bgzip minimap BED files 

"""

rule bed_bgzip_strains:
    localrule: True
    input:
        config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        protected(config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed.gz"),
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_bgzip_strains/{sid}/{sid}-{tp}.log",
    wrapper:
        "v5.8.2/bio/bgzip"



checkpoint minimap_bam_to_bed_precursors_old_genome:
    input:
        config['resultDir']  + "/minimap2/oldGenome/{sid}/{sid}-{tp}.500.bam"
    output:
        temp(config['resultDir']  + "/minimap2/oldGenome/{sid}/{sid}-{tp}.500.bed12.bed")
    threads: 1
    params:
        "-bed12 -i"
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
       config['logDir']  + "/minimap2_bam_to_bed_oldGenome/{sid}-{tp}.bed.log"
    shell: """
        bedtools bamtobed {params} {input} > {output}  2> {log}
    """

rule bed_bgzip_old_Genome:
    input:
        config['resultDir']  + "/minimap2/oldGenome/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        protected(config['resultDir']  + "/minimap2/oldGenome/{sid}/{sid}-{tp}.500.bed12.bed.gz"),
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_bgzip_strains/{sid}/{sid}-{tp}.log",
    wrapper:
        "v5.8.2/bio/bgzip"


rule bed12_bed6_strains:
    localrule: True
    input:
        config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed.gz"
    output:
        temp(config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.sort_bed6.bed"),
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 20000
    conda: config['envsDir'] + "/bedtools.yaml"
    log:
        config['logDir']  + "/bed12_bed6_strains/{sid}/{sid}-{tp}.log",
    shell:
        "bedtools bed12tobed6 -i {input} | sort -k 1,1 -k2,2n  > {output} 2> {log}"


"""

Indexing minimap BED files 

"""

rule bed_index_tabix_strains:
    localrule: True
    input:
        config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed.gz"
    output:
        config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed.gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a bed)
        "-p bed"
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_index_tabix_strains/{sid}-{tp}.log"
    wrapper:
        "0.51.3/bio/tabix"


####

#minimap2 refrance genome GRCm38.68

####

rule minimap_bam_to_bed_precursors_ref_GRCm38_68:
    localrule: True
    input:
        config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bam"
    output:
        config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bed12.bed"
    threads: 1
    params:
        "-bed12 -i"
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        queue  = 'normal',
        mem = 200000
    log:
       config['logDir']  + "/minimap2_precursor/GRCm38.68/{sid}-{tp}.bed.log"
    shell: """
        bedtools bamtobed {params} {input} > {output}  2> {log}
    """

rule bgzip_minimap_bam_to_bed_precursors_ref_GRCm38_68:
    input:
        config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        protected(config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bed12.bed.gz"),
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.log",
    wrapper:
        "v5.8.2/bio/bgzip"


rule bed12_bed6_GRCm38_68:
    input:
        config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bed12.bed.gz"
    output:
        temp(config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.sort_bed6.bed"),
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 20000
    conda: config['envsDir'] + "/bedtools.yaml"
    log:
        config['logDir']  + "/bed12_bed6_GRCm38_68/{sid}/{sid}-{tp}.log",
    shell:
        "bedtools bed12tobed6 -i {input} | sort -k 1,1 -k2,2n  > {output} 2> {log}"



####

#minimap2 refrance genome GRCm39.106

####

checkpoint minimap_bam_to_bed_precursors_ref_GRCm39_106:
    input:
        config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.bam"
    output:
        temp(config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.bed12.bed")
    threads: 1
    params:
        "-bed12 -i"
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        queue  = 'normal',
        mem = 200000
    log:
       config['logDir']  + "/minimap2_precursor/GRCm39.106/{sid}-{tp}.bed.log"
    shell: """
        bedtools bamtobed {params} {input} > {output}  2> {log}
    """

rule bgzip_minimap_bam_to_bed_precursors_ref_GRCm39_106:
    input:
        config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        protected(config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.bed12.bed.gz"),
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.log",
    wrapper:
        "v5.8.2/bio/bgzip"



rule bed12_bed6_GRCm39_106:
    input:
        config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.bed12.bed.gz"
    output:
        temp(config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.sort_bed6.bed"),
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 20000
    conda: config['envsDir'] + "/bedtools.yaml"
    log:
        config['logDir']  + "/bed12_bed6_GRCm39.106/{sid}/{sid}-{tp}.log",
    shell:
        "bedtools bed12tobed6 -i {input} | sort -k 1,1 -k2,2n  > {output} 2> {log}"





"""

bgzip zamore BED file 

"""

rule bed_bgzip_zamore:
    localrule: True
    input:
        config['resultDir']  + "/zamore/{ref}/piRNA_gene_annotation_mm10_without_chr.bed"
    output:
        protected(config['resultDir']  + "/zamore/{ref}/piRNA_gene_annotation_mm10_without_chr.bed.gz")
    params:
        extra="", # optional
    threads: 1
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_zamore_bgzip/piRNA_gene_annotation_{ref}_without_chr.bed.log",
    wrapper:
        """v5.8.2/bio/bgzip"""














####

#zamore lab piRNA 

####




"""

Indexing zamore BED file 

"""

rule bed_index_zamore:
    input:
        config['resultDir']  + "/zamore/{ref}/piRNA_gene_annotation_{ref}_without_chr.bed.gz"
    output:
        config['resultDir']  + "/zamore/{ref}/piRNA_gene_annotation_{ref}_without_chr.bed.gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p bed"
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_index_zamore_tabix/piRNA_gene_annotation_{ref}_without_chr.bed.log"
    wrapper:
        "v5.8.2/bio/tabix/index"

