

rule STAR_alingment_rna_strain_wise_gtf_testing:
    input:
        fq1 = [config['resultDir']  + "/cutadapt_pe_rna/C57BL_6NJ-{tp}.{rep}_r1.fastq"],
        fq2 = [config['resultDir']  + "/cutadapt_pe_rna/C57BL_6NJ-{tp}.{rep}_r2.fastq"],
        idx = config['resultDir'] + "/indexs/{sid}",
    output:
        bam=config['resultDir']  + "/STAR_rna_strain_wise_gtf_testing/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        log_final=config['resultDir']  + "/STAR_rna_strain_wise_gtf_testing/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    priority:5
    log:
        config['logDir']  + "/STAR_rna_strain_wise_gtf_testing/{sid}/{sid}-{tp}.{rep}.log"
    params:
        idx = config['resultDir'] + "/indexs/{sid}",
        extra = config['params']['rna']['STAR']
    wrapper:
        "v1.7.0/bio/star/align"


rule feature_counts_rna_strains_gtf_testing:
    input:
        sam =config['resultDir']  + "/STAR_rna_strain_wise_gtf_testing/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/featureCount_rna_gtf_testing/strains/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['rna']['featureCounts'] 
    log:
        config['logDir']  +"/featureCount_rna_gtf_testing/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"