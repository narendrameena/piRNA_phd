"""
downloding refrenace genome from ensamble;
0-based cordinates  
"""
rule get_strains_genomes:
    localrule: True
    output:
        config['resultDir'] + "/strain_genome/{sid}.fasta"
    cache: True
    resources:
        mem = 4000,
        queue = 'normal',
    threads: 1        
    params:
        url = lambda w: "ftp://mousegenomes_ftp:8ShghQQe@ftp-private.ebi.ac.uk/upload/REL-2205-Assembly/" + w.sid + "_chromosomes_MT.fasta.gz"
    shell:
        """
        wget -O {output}.gz {params.url}
        gunzip {output}.gz
        """


"""
getting annotation from ensamble 
0-based cordinates 
"""
rule get_strains_genome_annotation:
    localrule: True
    output:
        config['resultDir'] + "/strains_genome_annotation/{sid}.gff3"
    cache: True
    resources:
        mem = 4000,
        queue = 'normal',
    threads: 1        
    params:
        url = lambda w: "ftp://mousegenomes_ftp:8ShghQQe@ftp-private.ebi.ac.uk/upload/REL-2205-Assembly-Annotation/v3.5/" + w.sid + "_v3.5.gff3"
    shell:
        """
        wget -O {output} {params.url}
        """


rule get_strains_genome_repeatmasker:
    localrule: True
    output:
        config['resultDir'] + "/ref_genome_masked/{sid}_chromosomes_MT_unplaced.fasta.out"
    cache: True
    resources:
        mem = 4000,
        queue = 'normal',
    threads: 1        
    params:
        url = lambda w: "ftp://mousegenomes_ftp:8ShghQQe@ftp-private.ebi.ac.uk/upload/REL-2205-Assembly/masked/" + w.sid + "_chromosomes_MT_unplaced.fasta.out"
    shell:
        """
        wget -O {output} {params.url}
        """


rule get_strains_genome_liftOver_chain:
    localrule: True
    output:
        config['resultDir'] + "/ref_genome_chain/GRCm39_{sid}_chromosomes_MT_unplaced.chain"
    cache: True
    resources:
        mem = 4000,
        queue = 'normal',
    threads: 1        
    params:
        url = lambda w: "ftp://mousegenomes_ftp:8ShghQQe@ftp-private.ebi.ac.uk/upload/REL-2205-Assembly/GRCm39_chains/GRCm39_" + w.sid + "_chromosomes_MT_unplaced.chain"
    shell:
        """
        wget -O {output} {params.url}
        """
