####

#bed files from pacbio assemblies 

####


"""

genrating BED files from sRNA alinged BAM files pacbio assemblies 

"""

rule bed_file_from_srna_bam_strains:
    localrule: True
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        scriptDir = config['scriptDir']
    output:
        bed = temp(config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12.bed"),
        sam = temp(config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.sam"),
        bed_seq = temp(config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12_seq.bed"),
    threads: 1 
    resources:
        queue  = 'normal',
        mem = 200000
    priority:7
    log:
        stderr_bed = config['logDir']  + "/bedFileFromSrnaBam/{sid}-{tp}.{rep}.bed.err"
    conda: config['envsDir'] + "/bedtools.yaml"
    shell: """
        bedtools bamtobed -bed12 -i {input.bam} > {output.bed}  2> {log.stderr_bed}
        samtools view {input.bam} > {output.sam}
        paste {output.bed} {output.sam} | awk -f {input.scriptDir}/sh/bed_seq.awk > {output.bed_seq}
    """


"""

bgzip sRNA STAR alinged sRNA BED files 

"""

rule bed_bgzip_srna_strains:
    localrule: True
    input:
        config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12.bed"
    output:
        config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12.bed.gz"
    params:
        extra="", # optional
    threads: 1
    priority:8
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_sRNA_bgzip/{sid}/{sid}-{tp}.{rep}.log",
    wrapper:
        """v5.8.2/bio/bgzip"""



rule bed_seq_bgzip_srna_strains:
    localrule: True
    input:
        config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12_seq.bed"
    output:
        config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12_seq.bed.gz"
    params:
        extra="", # optional
    threads: 1
    priority:7
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_sRNA_bgzip/{sid}/{sid}-{tp}.{rep}_seq.log",
    wrapper:
        """v5.8.2/bio/bgzip"""







"""




Indexing sRNA STAR alinged BED files 

"""

checkpoint bed_index_srna_strains_tabix:
    localrule: True
    input:
        config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12.bed.gz"
    output:
        config['resultDir']  + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12.bed.gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p bed"
    resources:
        queue  = 'normal',
        mem = 200000
    priority:7
    log:
        config['logDir']  + "/bed_index_srna_tabix/{sid}-{tp}.{rep}.log"
    wrapper:
        """
        "v5.8.2/bio/tabix/index"

        """

"""

trinity assemblies 
"""



rule bed_file_from_srna_bam_trinity:
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise/trinity/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"
    output:
        bed = temp(config['resultDir']  + "/bedFileFromSrnaBam/trinity/{sid}/{sid}-{tp}.{rep}.bed12.bed")
    threads: 1 
    resources:
        queue  = 'normal',
        mem = 200000
    priority:9
    log:
        stderr_bed = config['logDir']  + "/bedFileFromSrnaBam/trinity/{sid}-{tp}.{rep}.bed.err"
    conda: config['envsDir'] + "/bedtools.yaml"
    shell: """
        bedtools bamtobed -bed12 -i {input.bam} > {output.bed}  2> {log.stderr_bed}
    """


"""

bgzip sRNA STAR alinged sRNA BED files 

"""

rule bed_bgzip_srna_trinity:
    input:
        config['resultDir']  + "/bedFileFromSrnaBam/trinity/{sid}/{sid}-{tp}.{rep}.bed12.bed"
    output:
        config['resultDir']  + "/bedFileFromSrnaBam/trinity/{sid}/{sid}-{tp}.{rep}.bed12.bed.gz"
    params:
        extra="", # optional
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_sRNA_bgzip/trinity_assemblies/{sid}/{sid}-{tp}.{rep}.log",
    wrapper:
        """v5.8.2/bio/bgzip"""

"""

Indexing sRNA STAR alinged BED files 

"""

checkpoint bed_index_srna_trinity_tabix:
    input:
        config['resultDir']  + "/bedFileFromSrnaBam/trinity/{sid}/{sid}-{tp}.{rep}.bed12.bed.gz"
    output:
        config['resultDir']  + "/bedFileFromSrnaBam/trinity/{sid}/{sid}-{tp}.{rep}.bed12.bed.gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p bed"
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        config['logDir']  + "/bed_index_srna_tabix/trinity/{sid}-{tp}.{rep}.log"
    priority:11
    wrapper:
        """
        "0.51.3/bio/tabix"

        """





###