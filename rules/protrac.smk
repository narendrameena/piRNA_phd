rule protrac_pirna_clusters_pacBio:
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        genomefile = config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta",
        gtf =  config['resultDir']  + "/annotationGTF/{sid}_v3.2.gtf", 
        repeatmasker = config['resources']  + "/repeatMasker/{sid}_chromosomes_MT_unplaced.fasta.out",
        script = config['scriptDir'] + "/perl/proTrac/proTRAC_2.4.3.pl", 
    output:
       sam = temp(config['resultDir']  + "/protrac/pacBio/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.sam"),
       dir = directory(config['resultDir']  + "/protrac/pacBio/{sid}/{sid}-{tp}.{rep}"),
       pti = protected(config['resultDir']  + "/protrac/pacBio/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.sam.sorted.ELAND3.pTi"),
    params:
        max_quality = 60
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 400000
    conda: 
        config['envsDir'] + "/protrac.yaml"
    log:
        config['logDir']  + "/protrac/pacBio/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        """
        samtools view -h -o {output.sam} {input.bam}
        cd {output.dir}
        perl {input.script} -map {output.sam} -genome {input.genomefile} -format SAM -repeatmasker {input.repeatmasker}  -geneset {input.gtf}  -pti > {log}    
        """



rule protrac_pirna_clusters_old_genome:
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        genomefile = config['resources'] + '/old_genomes/REL-1509-Assembly/' + "{sid}.chromosomes.unplaced.gt2k.fa",
        gtf = config['resultDir']  + "/annotationGTF/old_genome/{sid}_v3.2.gtf",
        repeatmasker = config['resources']  + "/old_genomes/repeatmasker/{sid}.chromosomes.unplaced.gt2k.fa.out",
        script = config['scriptDir'] + "/perl/proTrac/proTRAC_2.4.3.pl",
    output:
       sam = temp(config['resultDir']  + "/protrac/oldGenome/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.sam"),
       dir = directory(config['resultDir']  + "/protrac/oldGenome/{sid}/{sid}-{tp}.{rep}"),
       pti = protected(config['resultDir']  + "/protrac/oldGenome/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.sam.sorted.ELAND3.pTi"),
    params:
        max_quality = 60
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 200000
    conda: 
        config['envsDir'] + "/shortstack.yaml"
    log:
       config['logDir']  +"/protrac/oldGenome/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        """
        samtools view -h -o {output.sam} {input.bam}
        cd {output.dir}
        perl {input.script} -map {output.sam} -genome {input.genomefile} -format SAM -repeatmasker {input.repeatmasker}  -geneset {input.gtf}  -pti > {log}
        """



rule protrac_pirna_clusters_mm39:
    localrule: True
    input:
        bam = config['resultDir']  + "/STAR_srna_GRCm39_106/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        genomefile = config['resultDir']   + "/ref_genome/GRCm39.106.fasta",
        gtf =  config['resultDir']  + "/ref_genome/GRCm39.106.gtf", 
        repeatmasker = config['resultDir']  + "/ref_genome/mm39.fa.out",
        script = config['scriptDir'] + "/perl/proTrac/proTRAC_2.4.3.pl", 
    output:
       sam = temp(config['resultDir']  + "/protrac/mm39/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.sam"),
       dir = directory(config['resultDir']  + "/protrac/mm39/{sid}/{sid}-{tp}.{rep}"),
       pti = protected(config['resultDir']  + "/protrac/mm39/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.sam.sorted.ELAND3.pTi"),
    params:
        max_quality = 60
    threads: 1
    priority:10
    resources:
        queue  = 'normal',
        mem = 200000
    conda: 
        config['envsDir'] + "/protrac.yaml"
    log:
        config['logDir']  + "/protrac/mm39/{sid}/{sid}-{tp}.{rep}.log"    
    shell:
        """
        samtools view -h -o {output.sam} {input.bam}
        cd {output.dir}
        perl {input.script} -map {output.sam} -genome {input.genomefile} -format SAM -repeatmasker {input.repeatmasker}  -geneset {input.gtf}  -pti > {log}    
        """
