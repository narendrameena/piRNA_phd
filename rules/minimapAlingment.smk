
"""

Mapping filtred trinity precursors to PacBio assemblies using minimap2 

"""

checkpoint minimap_precursor_alingment_strains:
    localrule: True
    input:
        target = config['strains_genomes_path'] + "/{sid}_chromosomes_MT.fasta",
        query = [config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta"]
    output:
        bam = temp(config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bam"),
    threads: 4
    params:
        extra="-ax splice",          
        sorting="coordinate",          
        sort_extra="" 
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/minimap2/strains/{sid}-{tp}.bam.log"
    wrapper:
        "v5.8.2/bio/minimap2/aligner"



####

#minimap2 with respect to refrenace genomes 

####


rule minimap_precursor_alingment_ref_GRCm38_68:
    localrule: True
    input:
        target = config['resultDir']  + "/ref_genome/GRCm38.68.fasta",
        query = [config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta"]
    output:
        bam = temp(config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bam"),
    threads: 4
    params:
        extra="-ax splice",          
        sorting="coordinate",          
        sort_extra="",
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/minimap2/GRCm38.68/{sid}-{tp}.bam.log"
    wrapper:
        "v5.8.2/bio/minimap2/aligner"

rule minimap_precursor_alingment_ref_GRCm39_106:
    localrule: True
    input:
        target = config['resultDir']  + "/ref_genome/GRCm39.106.fasta",
        query = [config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta"]
    output:
        bam = temp(config['resultDir']  + "/minimap2/GRCm39.106/{sid}/{sid}-{tp}.500.bam"),
    threads: 4
    params:
        extra="-ax splice",          
        sorting="coordinate",          
        sort_extra="",
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/minimap2/GRCm39.106/{sid}-{tp}.bam.log"
    wrapper:
        "v5.8.2/bio/minimap2/aligner"



rule minimap_precursor_alingment_6NJ:
    localrule: True
    input:
        target = config['strains_genomes_path'] + "/C57BL_6NJ_chromosomes_MT.fasta",
        query = [config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta"]
    output:
        bam = temp(config['resultDir']  + "/minimap2/C57BL_6NJ/{sid}/{sid}-{tp}.500.bam"),
    threads: 4
    params:
        extra="-ax splice",          
        sorting="coordinate",          
        sort_extra="" 
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/minimap2/C57BL_6NJ/{sid}-{tp}.bam.log"
    wrapper:
        "v5.8.2/bio/minimap2/aligner"

rule minimap_precursor_alingment_SPRET_EiJ:
    localrule: True
    input:
        target = config['strains_genomes_path'] + "/SPRET_EiJ_chromosomes_MT.fasta",
        query = [config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta"]
    output:
        bam = temp(config['resultDir']  + "/minimap2/SPRET_EiJ/{sid}/{sid}-{tp}.500.bam"),
    threads: 4
    params:
        extra="-ax splice",          
        sorting="coordinate",          
        sort_extra="" 
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/minimap2/SPRET_EiJ/{sid}-{tp}.bam.log"
    wrapper:
        "v5.8.2/bio/minimap2/aligner"

checkpoint minimap_precursor_alingment_old_genome:
    input:
        target = config['resources'] + '/old_genomes/REL-1509-Assembly/' + "{sid}.chromosomes.unplaced.gt2k.fa",
        query = [config['resultDir']  + "/trinity/trinity_thres/trinity-{sid}-{tp}.Trinity.500.fasta"]
    output:
        bam = temp(config['resultDir']  + "/minimap2/oldGenome/{sid}/{sid}-{tp}.500.bam"),
    threads: 4
    params:
        extra="-x asm20",          
        sorting="coordinate",          
        sort_extra="" 
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/minimap2/oldGenome/{sid}-{tp}.bam.log"
    wrapper:
        "v5.8.2/bio/minimap2/aligner"