rule STAR_alingment_uniqe_marker_strain:
    localrule: True
    input:
        fq1 =  config['resources']  + "/uniqMakers/{sid}-{mr}.fasta",
        idx = config['resultDir'] + "/indexs/{sid}",
    output:
        bam=config['resultDir']  + "/STAR_uniqe_marker_strain_wise/uniq_alinged/{sid}/{sid}-{mr}/Aligned.sortedByCoord.out.bam",
        log_final=config['resultDir']  + "/STAR_uniqe_marker_strain_wise/uniq_alinged/{sid}/{sid}-{mr}/Log.final.out"
    threads: 4
    priority:3
    log:
        config['logDir']  + "/STAR_uniqe_marker_strain_wise/uniq_alinged/{sid}/{sid}-{mr}.log"
    resources:
        mem = 400000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        idx = config['resultDir'] + "/indexs/{sid}",
        extra = config['params']['srna']['STAR'] 
    wrapper:
        "v5.8.2/bio/star/align"


rule extract_uniq_marker_cordinates_reads_index:
    localrule: True
    input:
        config['resultDir']   + "/STAR_uniqe_marker_strain_wise/uniq_alinged/{sid}/{sid}-{mr}/Aligned.sortedByCoord.out.bam"
    output:
        config['resultDir']  + "/STAR_uniqe_marker_strain_wise/uniq_alinged/{sid}/{sid}-{mr}/Aligned.sortedByCoord.out.bam.bai",
    log:
        config['logDir']  + "/STAR_uniqe_marker_strain_wise/indexbam/{sid}-{mr}.log"
    resources:
        mem = 200000,
        queue  = 'normal',
    params:
        extra="",
    threads: 1
    priority:4
    wrapper:
        "v5.8.2/bio/samtools/index"

rule extract_uniq_marker_cordinates_reads:
    localrule: True
    input:
        aligned_bam = config['resultDir'] + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        bam = config['resultDir'] + "/STAR_uniqe_marker_strain_wise/uniq_alinged/{sid}/{sid}-{mr}/Aligned.sortedByCoord.out.bam",
        bai = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/uniq_alinged/{sid}/{sid}-{mr}/Aligned.sortedByCoord.out.bam.bai",
        script = config['scriptDir'] + "/sh/extractBam.sh"
    output:
        bam = config['resultDir'] + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam"
    log:
        config['logDir'] + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}.log"
    resources:
        mem = 200000,
        queue = 'normal',
        tmpdir = config['tmpDir']
    priority: 15
    conda: config['envsDir'] + "/samtools.yaml"
    shell:
        """
        bash {input.script} {input.bam} {input.aligned_bam} {output.bam} > {log} 2>&1
        """


rule STAR_strains_uniqe_marker_bam_index:
    localrule: True
    input:
        config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam"
    output:
        config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam.bai",
    log:
        config['logDir']  + "/STAR_uniqe_marker_strain_wise_idx/indexbam/{sid}-{tp}.{rep}-{mr}.log"
    resources:
        mem = 200000,
        queue  = 'normal',
    params:
        extra="",
    threads: 4
    priority:4
    wrapper:
        "v5.8.2/bio/samtools/index"


rule strains_uniq_marker_collapse_aligned_bam_srna:
    localrule: True
    input:
        bam=config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam",
        scriptDir = config['scriptDir']
    output:
        temp(config['resultDir']  + "/STAR_uniqe_marker_strain_wise/collapse/{sid}-{tp}.{rep}-{mr}.aligned.fasta.gz")
    log:
        stdout=config['logDir']  + "/STAR_uniqe_marker_strain_wise/collapse/{sid}-{tp}.{rep}-{mr}.aligned.out",
        stderr=config['logDir']  + "/STAR_uniqe_marker_strain_wise/collapse/{sid}-{tp}.{rep}-{mr}.aligned.err"
    threads: 1
    priority:2
    resources:
        mem = 200000,
        queue  = 'normal',
    conda: config['envsDir'] + "/tstk.yaml"
    shell:
        "collapse  {input.bam} {output} > {log.stdout} 2> {log.stderr}"
   

rule strains_collapse_uniqe_markers_bam_count_srna:
    localrule: True
    input:
        bam = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/collapse/{sid}-{tp}.{rep}-{mr}.aligned.fasta.gz",
        bai = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam.bai"
    output:
        config['resultDir']  + "/STAR_uniqe_marker_strain_wise/seqCount/{sid}/{sid}-{tp}.{rep}-{mr}.uniqe_marker.count.gz"
    log:
        stderr=config['logDir']  + "/STAR_uniqe_marker_strain_wise/seqCount/{sid}/{sid}-{tp}.{rep}-{mr}.uniqe_marker.err"
    threads: 1
    priority:15
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    conda: config['envsDir'] + "/seqkit.yaml"
    shell:
        r"""
        seqkit fx2tab {input.bam} | awk '{{split($1,a,"-"); print $2 "\t" a[1] "\t" a[2]}}' | sed  '1i seq \t seq_id \t {wildcards.sid}-{wildcards.tp}.{wildcards.rep}-{wildcards.mr}' |  gzip --stdout > {output} 2> {log.stderr}
        """




rule bam_to_bed_uniqe_marker:
    localrule: True
    input:
        bam=config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam",
        bai = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam.bai"
    output:
        bed=config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedbed/{sid}/{sid}-{tp}.{rep}-{mr}/{sid}-{tp}.{rep}-{mr}.bed"
    log:
        config['logDir']  + "/STAR_uniqe_marker_strain_wise/extractedbed/{sid}/{sid}-{tp}.{rep}-{mr}.log"
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    conda: config['envsDir'] + "/bedtools.yaml"
    priority:15
    conda: config['envsDir'] + "/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input.bam} > {output.bed} 2> {log}"



rule convert_bam_to_bigwig_uniqe_marker:
    localrule: True
    input:
        bam = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        bai = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam.bai",
    output:
        forward  = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedbigwig/{sid}/{sid}-{tp}.{rep}-{mr}_plusStrand.bw",
        reverse1 =  config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedbigwig/{sid}/{sid}-{tp}.{rep}-{mr}_minusStrand.bw"
    conda: config['envsDir'] + "/bigwig.yaml"
    params:
        normlize = '--normalizeUsing CPM',
        script = config['scriptDir'] + "/R/" +"bamTobigWig.R"
    threads: 1
    priority:15
    resources:
        queue  = 'normal',
        mem = 200000
    log:
        forward = config['logDir']  +"/STAR_uniqe_marker_strain_wise/extractedbigwig/{sid}/{sid}-{tp}.{rep}-{mr}_forward.log",
        reverse1 = config['logDir']  +"/STAR_uniqe_marker_strain_wise/extractedbigwig/{sid}/{sid}-{tp}.{rep}-{mr}_reverse.log"
    shell:
        """
        bamCoverage  {params.normlize} -b {input.bam} -o {output.forward}  --filterRNAstrand forward > {log.forward} &
        bamCoverage  {params.normlize} -b {input.bam} -o {output.reverse1} --filterRNAstrand reverse > {log.reverse1}
        """

rule feturecount_uniqe_marker_TE_strains:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/TEannotation/pacBio/{sid}_v3.gff3",
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/STAR_uniqe_marker_strain_wise/fetureCountTE/{sid}/{sid}-{tp}.{rep}-{mr}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:1
    priority:15
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['TE_rpkm']['featureCounts'] 
    log:
        config['logDir']  +"/STAR_uniqe_marker_strain_wise/fetureCountTE/{sid}/{sid}-{tp}.{rep}-{mr}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"



rule intron_exon_uniqe_marker_counts_pacBio_list1:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/intergenic/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/STAR_uniqe_marker_strain_wise/fetureCountlist1/{sid}/{sid}-{tp}.{rep}-{mr}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:15
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list1'] 
    log:
        config['logDir']  +"/STAR_uniqe_marker_strain_wise/fetureCountlist1/{sid}/{sid}-{tp}.{rep}-{mr}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


rule intron_exon_uniqe_marker_counts_pacBio_list2:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/intron/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/STAR_uniqe_marker_strain_wise/fetureCountlist2/{sid}/{sid}-{tp}.{rep}-{mr}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:15
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list2'] 
    log:
        config['logDir']  +"/STAR_uniqe_marker_strain_wise/fetureCountlist2/{sid}/{sid}-{tp}.{rep}-{mr}.log"
    wrapper:
        "0.72.0/bio/subread/featurecounts"

rule intron_exon_uniqe_marker_counts_pacBio_list3:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/gene_intron/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/STAR_uniqe_marker_strain_wise/fetureCountlist3/{sid}/{sid}-{tp}.{rep}-{mr}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:15
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list3'] 
    log:
        config['logDir']  +"/STAR_uniqe_marker_strain_wise/fetureCountlist3/{sid}/{sid}-{tp}.{rep}-{mr}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


rule intron_exon_uniqe_marker_counts_pacBio_list4:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_uniqe_marker_strain_wise/extractedBam/{sid}/{sid}-{tp}.{rep}-{mr}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/lncRNA_intron/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/STAR_uniqe_marker_strain_wise/fetureCountlist4/{sid}/{sid}-{tp}.{rep}-{mr}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:15
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list4'] 
    log:
        config['logDir']  +"/STAR_uniqe_marker_strain_wise/fetureCountlist4/{sid}/{sid}-{tp}.{rep}-{mr}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"