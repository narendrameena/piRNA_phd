rule bed2_file_gff3_piRNA:
    localrule: True
    input:
        script = config['scriptDir'] + "/python/bedtogff3Convert.py",
        bed = config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        annotation = config['resultDir']  + "/bedToGff3PiRNA/{sid}_{tp}_piRNA_v3.2.gff3",
    threads:
        1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/bed_ggf3_piRNA/intergenic/{sid}_{tp}.log"
    shell:"python {input.script} {input.bed} {output.annotation} > {log}"


rule piRNA_proporation_featurecount:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/bedToGff3PiRNA/{sid}_{tp}_piRNA_v3.2.gff3",
    output:
        protected(multiext(config['resultDir']  + "/piRNA_proporation/{sid}_{tp}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:17
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['piRNA_count'] 
    log:
        config['logDir']  +"/piRNA_proporation/{sid}_{tp}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


rule zamore_09_bed2_file_gff3_piRNA:
    localrule: True
    input:
        script = config['scriptDir'] + "/python/bedtogff3Convert.py",
        bed = config['resultDir']  + "/zamore/with_09/{sid}/piRNA_gene_annotation_{sid}_without_chr.bed"
    output:
        annotation = config['resultDir']  + "/bedToGff3PiRNA/zamore_09/{sid}_{tp}_piRNA_v3.2.gff3",
    threads:
        1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/bed_ggf3_piRNA/zamore_09/{sid}_{tp}.log"
    shell:"python {input.script} {input.bed} {output.annotation} > {log}"


rule zamore_09_piRNA_proporation_featurecount:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/bedToGff3PiRNA/zamore_09/{sid}_{tp}_piRNA_v3.2.gff3",
    output:
        protected(multiext(config['resultDir']  + "/piRNA_proporation/zamore_09/{sid}_{tp}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:17
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['piRNA_count'] 
    log:
        config['logDir']  +"/piRNA_proporation/zamore_09/{sid}_{tp}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


rule zamore_05_bed2_file_gff3_piRNA:
    localrule: True
    input:
        script = config['scriptDir'] + "/python/bedtogff3Convert.py",
        bed = config['resultDir']  + "/zamore/with_05/{sid}/piRNA_gene_annotation_{sid}_without_chr.bed"
    output:
        annotation = config['resultDir']  + "/bedToGff3PiRNA/zamore_05/{sid}_{tp}_piRNA_v3.2.gff3",
    threads:
        1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/bed_ggf3_piRNA/zamore_05/{sid}_{tp}.log"
    shell:"python {input.script} {input.bed} {output.annotation} > {log}"


rule zamore_05_piRNA_proporation_featurecount:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/bedToGff3PiRNA/zamore_05/{sid}_{tp}_piRNA_v3.2.gff3",
    output:
        protected(multiext(config['resultDir']  + "/piRNA_proporation/zamore_05/{sid}_{tp}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:17
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['piRNA_count'] 
    log:
        config['logDir']  +"/piRNA_proporation/zamore_05/{sid}_{tp}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


rule zamore_02_bed2_file_gff3_piRNA:
    localrule: True
    input:
        script = config['scriptDir'] + "/python/bedtogff3Convert.py",
        bed = config['resultDir']  + "/zamore/with_02/{sid}/piRNA_gene_annotation_{sid}_without_chr.bed"
    output:
        annotation = config['resultDir']  + "/bedToGff3PiRNA/zamore_02/{sid}_{tp}_piRNA_v3.2.gff3",
    threads:
        1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/bed_ggf3_piRNA/zamore_02/{sid}_{tp}.log"
    shell:"python {input.script} {input.bed} {output.annotation} > {log}"


rule zamore_02_piRNA_proporation_featurecount:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/bedToGff3PiRNA/zamore_02/{sid}_{tp}_piRNA_v3.2.gff3",
    output:
        protected(multiext(config['resultDir']  + "/piRNA_proporation/zamore_02/{sid}_{tp}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:
        2
    priority:17
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['intron_rpkm']['piRNA_count'] 
    log:
        config['logDir']  +"/piRNA_proporation/zamore_02/{sid}_{tp}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"

rule final_filtred_precursors_gtf:
    localrule: True
    input:
        bed = config['resultDir'] + "/filtred_final_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}.bed",
        script = config['scriptDir'] + "/sh/bedToGtf.sh"
    output:
        config['resultDir'] + "/filtred_final_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}.gtf"
    threads: 1
    conda: config['envsDir'] + "/refrancebasedprecursors.yaml"
    resources:
        mem = 200000,
        queue = 'normal',
        tmpdir = config['tmpDir']
    priority: 10
    log:
        config['logDir'] + "/filtred_precursors/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        bash {input.script} {input.bed} {output} > {log}
        """

rule final_piRNA_proporation_featurecount:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir'] + "/filtred_final_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}.gtf"
    output:
        protected(multiext(config['resultDir']  + "/final_precursors_piRNA_proporation/{sid}/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:  2
    priority:17
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        extra = config['params']['intron_rpkm']['piRNA_count'] 
    log:
        config['logDir']  +"/final_precursors_piRNA_proporation/{sid}_{tp}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"
