
rule process_bam_file_for_sliding_window_strains:
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        script = config['scriptDir'] + "/python/" +"slidingWindowOnBam.py",
        bai = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai"
    output:
        temp(config['resultDir']  + "/slidingWindow/pacBio/{sid}/{sid}-{tp}.{rep}.tab")
    params:
        max_quality = 60
    threads: 1
    priority:9
    resources:
        queue  = 'normal',
        mem = 200000
    conda: 
        config['envsDir'] + "/slidingWindow.yaml"
    log:
        config['logDir']  +"/slidingWindow/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        "python {input.script} --input_file {input.bam} --output_file {output} > {log}"

        

"""

bgzip sliding window file 

"""

rule bgzip_bed_bam_file_for_sliding_window_strains:
    input:
        config['resultDir']  + "/slidingWindow/pacBio/{sid}/{sid}-{tp}.{rep}.tab"
    output:
        protected(config['resultDir']  + "/slidingWindow/pacBio/{sid}/{sid}-{tp}.{rep}.tab.gz"),
    params:
        extra="", # optional
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 20000
    log:
        config['logDir']  + "/sliding_window/{sid}/{sid}-{tp}.{rep}_bzip.log",
    wrapper:
        "v1.7.0/bio/bgzip"



### old genomes 

rule process_bam_file_for_sliding_window_old_genomes:
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        script = config['scriptDir'] + "/python/" +"slidingWindowOnBam.py",
        bai = config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai"
    output:
        temp(config['resultDir']  + "/slidingWindow/old_genome/{sid}/{sid}-{tp}.{rep}.tab")
    params:
        max_quality = 60
    threads: 1
    priority:9
    resources:
        queue  = 'normal',
        mem = 200000
    conda: 
        config['envsDir'] + "/slidingWindow.yaml"
    log:
        config['logDir']  +"/slidingWindow/old_genome/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        "python {input.script} --input_file {input.bam} --output_file {output} > {log}"

        

"""

bgzip sliding window file 

"""

rule bgzip_bed_bam_file_for_sliding_window_old_genomes:
    input:
        config['resultDir']  + "/slidingWindow/old_genome/{sid}/{sid}-{tp}.{rep}.tab"
    output:
        protected(config['resultDir']  + "/slidingWindow/old_genome/{sid}/{sid}-{tp}.{rep}.tab.gz"),
    params:
        extra="", # optional
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 20000
    log:
        config['logDir']  + "/sliding_window/old_genome/{sid}/{sid}-{tp}.{rep}_bzip.log",
    wrapper:
        "v1.7.0/bio/bgzip"



#shortstack for pirna cluster identifications 


rule shortstack_pirna_clusters_pacBio:
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        genomefile = config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta"
    output:
       dir = directory(config['resultDir']  + "/shortstack/pacBio/{sid}/{sid}-{tp}.{rep}"),
       result = config['resultDir']  + "/shortstack/pacBio/{sid}/{sid}-{tp}.{rep}/Results.txt"
    params:
        max_quality = 60
    threads: 4
    priority:10
    resources:
        queue  = 'normal',
        mem = 200000
    conda: 
        config['envsDir'] + "/shortstack.yaml"
    log:
        config['logDir']  +"/shortstack/pacBio/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        """
        rm -r {output.dir}
        ShortStack --genomefile  {input.genomefile} --bamfile {input.bam} --outdir {output.dir}  > {log}
        """



rule shortstack_pirna_clusters_old_genome:
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        genomefile = config['resources'] + '/old_genomes/REL-1509-Assembly/' + "{sid}.chromosomes.unplaced.gt2k.fa"
    output:
       dir = directory(config['resultDir']  + "/shortstack/oldGenome/{sid}/{sid}-{tp}.{rep}"),
       result = config['resultDir']  + "/shortstack/oldGenome/{sid}/{sid}-{tp}.{rep}/Results.txt"
    params:
        max_quality = 60
    threads: 4
    priority:10
    resources:
        queue  = 'normal',
        mem = 200000
    conda: 
        config['envsDir'] + "/shortstack.yaml"
    log:
        config['logDir']  +"/shortstack/oldGenome/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        """
        rm -r {output.dir}
        ShortStack --genomefile  {input.genomefile} --bamfile {input.bam} --outdir {output.dir}  > {log}
        """



rule shortstack_pirna_clusters_mm39:
    input:
        bam = config['resultDir']  + "/STAR_srna_GRCm39_106/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        genomefile = config['resultDir']  + "/ref_genome/GRCm39.106.fasta"
    output:
       dir = directory(config['resultDir']  + "/shortstack/mm39/{sid}/{sid}-{tp}.{rep}"),
       result = config['resultDir']  + "/shortstack/mm39/{sid}/{sid}-{tp}.{rep}/Results.txt"
    params:
        max_quality = 60
    threads: 4
    priority:10
    resources:
        queue  = 'normal',
        mem = 200000
    conda: 
        config['envsDir'] + "/shortstack.yaml"
    log:
        config['logDir']  +"/shortstack/mm39/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        """
        rm -r {output.dir}
        ShortStack --genomefile  {input.genomefile} --bamfile {input.bam} --outdir {output.dir}  > {log}
        """
