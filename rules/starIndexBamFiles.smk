
"""

Indexing alinged (with refrance to Pac Bio assemblies) BAM files

"""

rule STAR_strains_srna_bam_index:
    input:
        config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
    output:
        config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
    log:
        config['logDir']  + "/STAR_srna_strain_wise_idx/{sid}-{tp}.{rep}.log"
    resources:
        mem = 200000,
        queue  = 'normal',
    params:
        extra="",
    threads: 4
    priority:4
    wrapper:
        "v5.8.2/bio/samtools/index"



rule STAR_old_genome_srna_bam_index:
    input:
        config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
    output:
        config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
    log:
        config['logDir']  + "/STAR_srna_old_genome_wise_idx/{sid}-{tp}.{rep}.log"
    resources:
        mem = 200000,
        queue  = 'normal',
    params:
        extra="",
    threads: 4
    priority:4
    wrapper:
        "v5.8.2/bio/samtools/index"




