####
## feature count intron and exon level 
####

rule intergenic_list1_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/intergenic_list1_gff.py"
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/intergenic/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/genicRegionGff3/intergenic/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"


rule intron_exon_srna_counts_pacBio_list1:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/intergenic/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list1/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list1'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list1/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


rule intergenic_list2_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/intron_list2_gff.py"
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/intron/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/genicRegionGff3/intron/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"

rule intron_exon_srna_counts_pacBio_list2:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/intron/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list2/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list2'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list2/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"



#list2 exon 

rule list2_exon_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/list2_exon_gff.py"
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/list2_exon/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/genicRegionGff3/list2_exon/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"

rule exon_srna_counts_pacBio_list2:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/list2_exon/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list2_exon/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:      2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list2_exon'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list2_exon/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"



#list2 CDS 

rule list2_CDS_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/list2_CDS_gff.py",
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/list2_CDS/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/genicRegionGff3/list2_CDS/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"

rule cds_srna_counts_pacBio_list2:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/list2_CDS/{sid}_v3.2.gff3",
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list2_CDS/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:      2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list2_cds'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list2_CDS/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"

#list2 3utr

rule list2_3utr_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/list2_3utr_gff.py"
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/list2_3utr/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/genicRegionGff3/list2_3utr/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"

rule utr3_srna_counts_pacBio_list2:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/list2_3utr/{sid}_v3.2.gff3",
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list2_3utr/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:      2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list2_3utr'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list2_3utr/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


#list2 5utr

rule list2_5utr_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/list2_5utr_gff.py"
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/list2_5utr/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/genicRegionGff3/list2_3utr/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"

rule utr5_srna_counts_pacBio_list2:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/list2_5utr/{sid}_v3.2.gff3",
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list2_5utr/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:      2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['featureCounts_list2_5utr'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list2_5utr/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"



#list 3 

rule intergenic_list3_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/gene_intron_list3_gff.py"
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/gene_intron/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/genicRegionGff3/intron/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"


rule intron_exon_srna_counts_pacBio_list3:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/gene_intron/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list3/{sid}/{sid}-{tp}.{rep}",
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
        extra = config['params']['intron_rpkm']['featureCounts_list3'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list3/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


rule intergenic_list4_pacBio_gff3:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/python/lncRNA_intron_list4_gff.py"
    output:
        annotation = config['resultDir']  + "/genicRegionGff3/lncRNA_intron/{sid}_v3.2.gff3",
    threads:
        1
    priority:8
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['srna']['transcriptCounts'] 
    log:
        config['logDir']  +"/intronExongff/lncRNA_intron/{sid}.log"
    shell:"python {input.script} {input.annotation} {output.annotation} > {log}"

rule intron_exon_srna_counts_pacBio_list4:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/genicRegionGff3/lncRNA_intron/{sid}_v3.2.gff3",
        #awk '{ sub(/Parent=transcript/, "regex Parent=transcript:"$3) }1'
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/IntronExonfeatureCount/pacBio/list4/{sid}/{sid}-{tp}.{rep}",
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
        extra = config['params']['intron_rpkm']['featureCounts_list4'] 
    log:
        config['logDir']  +"/IntronExonfeatureCount/pacBio/list4/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"