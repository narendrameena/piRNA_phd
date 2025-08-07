

"""

coverage of trinity precursors from sRNA STAR alingment BAM files  

"""

rule bed_file_coverage_from_bam_srna:
    localrule: True
    input:
        b = config['resultDir']  + "/STAR_srna_strain_wise/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        a = config['resultDir']  + "/minimap2/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        temp(config['resultDir']  + "/bedFilesRNABamCoverage/{sid}/{sid}-{tp}.{rep}.cov")
    threads: 2
    priority:7
    resources:
        queue  = 'normal',
        mem = 20000
    log:
        config['logDir']  + "/bedFilesRNABamCoverage/{sid}-{tp}.{rep}.cov.log"
    params:
        extra="-s " 
    wrapper:
        "v1.7.0/bio/bedtools/coveragebed"



"""

bgzip minimap BED coverage File files 

"""

rule bed_cov_bgzip_strains:
    localrule: True
    input:
        config['resultDir']  + "/bedFilesRNABamCoverage/{sid}/{sid}-{tp}.{rep}.cov"
    output:
        config['resultDir']  + "/bedFilesRNABamCoverage/{sid}/{sid}-{tp}.{rep}.cov.gz",
    params:
        extra="", # optional
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 20000,
        tmpdir= config['tmpDir']
    log:
        config['logDir']  + "/bedFilesRNABamCoverage/{sid}/{sid}-{tp}.{rep}_bzip.log",
    wrapper:
        "v1.7.0/bio/bgzip"


