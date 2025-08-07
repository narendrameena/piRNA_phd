"""

samtools index 

"""


rule samtools_index_strains:
    input:
        config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta",
    output:
        temp(config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta.fai"),
    log:
        config['logDir']  + "/genomeFile/{sid}/{sid}.log"
    resources:
        queue  = 'normal',
        mem = 20000
    params:
        extra="",  # optional params string
    wrapper:
        "v1.7.0/bio/samtools/faidx"

rule genome_file_strains:
    input:
        config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta.fai"
    output:
        config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta.txt",
    log:
        config['logDir']  + "/genomeFile/{sid}/{sid}_awk.log"
    resources:
        queue  = 'normal',
        mem = 20000
    params:
        extra="",  # optional params string
    shell:
        "awk '{{print \$1,"",\$2}}'   {input} > {output} 2> {log}"


"""
GRCm38.68
"""

rule samtools_index_GRCm38_68:
    input:
        config['resultDir']  + "/ref_genome/GRCm38.68.fasta",
    output:
        temp(config['resultDir']  + "/ref_genome/GRCm38.68.fasta.fai"),
    resources:
        queue  = 'normal',
        mem = 20000
    log:
        config['logDir']  + "/genomeFile/GRCm38.68.log"
    params:
        extra="",  # optional params string
    wrapper:
        "v1.7.0/bio/samtools/faidx"

rule genome_file_GRCm38_68:
    input:
        config['resultDir']  + "/ref_genome/GRCm38.68.fasta.fai",
    output:
        config['resultDir']  + "/ref_genome/GRCm38.68.fasta.txt",
    log:
        config['logDir']  + "/genomeFile/GRCm38.68_awk.log"
    resources:
        queue  = 'normal',
        mem = 20000
    params:
        extra="",  # optional params string
    shell:
        "awk '{{print \$1,"",\$2}}'   {input} > {output} 2> {log}"

  
""" 
  GRCm39.106
"""  
  
      
rule samtools_index_GRCm39_106:
    input:
        config['resultDir']  + "/ref_genome/GRCm39.106.fasta",
    output:
        temp(config['resultDir']  + "/ref_genome/GRCm39.106.fasta.fai"),
    resources:
        queue  = 'normal',
        mem = 20000
    log:
        config['logDir']  + "/genomeFile/GRCm39.106.log"
    params:
        extra="",  # optional params string
    wrapper:
        "v5.8.2/bio/samtools/faidx"

rule genome_file_GRCm39_106:
    input:
        config['resultDir']  + "/ref_genome/GRCm39.106.fasta.fai",
    output:
        config['resultDir']  + "/ref_genome/GRCm39.106.fasta.txt",
    log:
        config['logDir']  + "/genomeFile/GRCm39.106_awk.log"
    resources:
        queue  = 'normal',
        mem = 20000
    params:
        extra="",  # optional params string
    shell:
        "awk '{{print \$1,"",\$2}}'   {input} > {output} 2> {log}"