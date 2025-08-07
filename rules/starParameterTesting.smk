


rule STAR_alingment_srna_strain_testing:
    input:
        fq1 =  config['resultDir']  + "/cutadapt_se/C57BL_6NJ-16.5dpc.1.fastq",
        idx = config['resultDir'] + "/indexs/C57BL_6NJ",
        #fq1 = config['resultDir']  + "/STAR_srna_strain_wise/tmp_fasta/C57BL_6NJ-16.5dpc.1.fasta"
    output:
        bam=config['resultDir']  + "/STAR_testing_pacbio/C57BL_6NJ/C57BL_6NJ-16.5dpc.1_{mis}_{multi}_{mAnc}/Aligned.sortedByCoord.out.bam",
        log_final=config['resultDir']  + "/STAR_testing_pacbio/C57BL_6NJ/C57BL_6NJ-16.5dpc.1_{mis}_{multi}_{mAnc}/Log.final.out"
    threads: 4
    priority:3
    log:
        config['logDir']  + "/STAR_testing_pacbio/C57BL_6NJ/C57BL_6NJ-16.5dpc.1_{mis}_{multi}_{mAnc}.log",
    resources:
        mem = 200000,
        queue  = 'normal',
    params:
        idx = config['resultDir'] + "/indexs/C57BL_6NJ",
        extra = config['params']['srna']['STAR_test'] + "  --outFilterMismatchNmax  {mis} --outFilterMultimapNmax {multi} --winAnchorMultimapNmax {mAnc}" 
    wrapper:
        "v5.8.2/bio/star/align"



rule feturecount_srna_TE_strains_testing:
    input:
        sam = config['resultDir']  + "/STAR_testing_pacbio/C57BL_6NJ/C57BL_6NJ-16.5dpc.1_{mis}_{multi}_{mAnc}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/TEannotation/pacBio/C57BL_6NJ_v3.gff3",
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/TEtranscriptCount_testing/pacBio/C57BL_6NJ/C57BL_6NJ-16.5dpc.1_{mis}_{multi}_{mAnc}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:1
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['TE_rpkm']['featureCounts'] 
    log:
        config['logDir']  +"/TEtranscriptCount_testing/C57BL_6NJ/C57BL_6NJ-16.5dpc.1_{mis}_{multi}_{mAnc}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


