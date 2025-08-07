
"""

merge bam files from STAR alingmetn of RNA seq data

"""

rule merge_bams_rna:
    input:
        [config['resultDir']  + f"/STAR_rna_strain_wise/{{sid}}-{{tp}}.{rep}/Aligned.sortedByCoord.out.bam" for rep in reps]
    output:
        config['resultDir']  + "/merge_bams/{sid}-{tp}.bam"
    params:
        extra=""
    threads: 4
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/merge_bams/{sid}-{tp}.log"
    wrapper:
        "v5.8.2/bio/samtools/merge"



"""

Index RNA merge bam files 

"""

rule STAR_strains_rna_merge_bam_index:
    input:
        config['resultDir']  + "/merge_bams/{sid}-{tp}.bam",
    output:
        config['resultDir']  + "/merge_bams/{sid}-{tp}.bam.bai",
    log:
        config['logDir']  + "/merge_bams/{sid}-{tp}.index.log"
    resources:
        mem = 200000,
        queue  = 'normal',
    params:
        extra="",
    threads: 4
    wrapper:
        "v5.8.2/bio/samtools/index"