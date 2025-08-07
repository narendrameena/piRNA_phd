"""
Indexing of reference genomes using STAR
"""

rule STAR_index_reference_genome:
    input:
        fasta=config['resultDir'] + "/ref_genome/{ref_genome}.fasta"
    output:
        directory(config['resultDir'] + "/indexes/{ref_genome}")
    threads: 4
    params:
        extra="--limitGenomeGenerateRAM 96000000000"
    message:
        "**** STAR index for {wildcards.ref_genome} ****"
    resources:
        mem_mb=200000,
        queue='normal'
    log:
        config['logDir'] + "/indexes/{ref_genome}.log"
    wrapper:
        "v5.8.2/bio/star/index"

"""
Indexing old assemblies from Thomas Keen's lab
"""

rule STAR_index_old_genomes:
    input:
        fasta=config['resources'] + "/old_genomes/REL-1509-Assembly/{sid}.chromosomes.unplaced.gt2k.fa"
    output:
        directory(config['resources'] + "/old_genomes/index/{sid}")
    threads: 4
    params:
        extra="--limitGenomeGenerateRAM 96000000000"
    message:
        "**** STAR index for old genome {wildcards.sid} ****"
    resources:
        mem_mb=200000,
        queue='normal'
    log:
        config['logDir'] + "/STAR_index_strain_wise_old_genomes/{sid}.log"
    wrapper:
        "v5.8.2/bio/star/index"

"""
Indexing PacBio assemblies using STAR
"""

rule remove_MT_from_genome:
    input:
        fasta=config['strains_genomes_path'] + "/{sid}_chromosomes_MT.fasta"
    output:
        fasta=temp(config['strains_genomes_path'] + "/{sid}_chromosomes.fasta")
    log:
        config['logDir'] + "/remove_MT_from_genome/{sid}.log"
    resources:
        mem_mb=2000,
        queue='normal',
        tmpdir=config['tmpDir']
    shell:
        """
        set -euo pipefail
        grep -v 'MT' "{input.fasta}" > "{output.fasta}" 2> "{log}"
        """


rule STAR_index_pacbio_genome:
    input:
        fasta=config['resources'] + "/REL-2205-Assembly/{sid}_chromosomes_MT.fasta",
        gtf=lambda wildcards: config['resources'] + f"/annotation/{wildcards.sid}_v3.3.gff3"
    output:
        directory(config['resultDir']  + "/indexs/{sid}")
    threads: 4
    params:
        extra = "--limitGenomeGenerateRAM  96000000000",
    message:
        "****STAR index****"
    priority:2
    resources:
        mem = 200000,
        queue  = '2004*',
    log:
        config['logDir']  + "/PICB/index/{sid}.log",
    wrapper:
        "v5.8.2/bio/star/index"