

rule gff3_to_gtf_pacBio:
    input:
        gff3 = config['resultDir'] + "/strains_genome_annotation/{sid}.gff3",
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        gtf = temp(config['resultDir']  + "/annotationGTF/{sid}.gtf")
    threads:1
    priority:10
    conda: config['envsDir'] + "/tpmCalculator.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/gff3_to_gtf__rna/{sid}.log"
    shell:
        "agat_convert_sp_gff2gtf.pl --gff {input.gff3} -o {output.gtf} > {log} "




#RNA seq data count gene TPM values 
rule TPMCount_counts_rna_strains:
    input:
        sam = config['resultDir']  + "/STAR_rna_strain_wise/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/annotationGTF/{sid}.gtf", 
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        dir = config['resultDir']  + "/TPMCount_rna/strains/{sid}-{tp}.{rep}/",
        tpm = protected(config['resultDir']  + "/TPMCount_rna/strains/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out_genes.uni")
    threads:1
    priority:10
    conda: config['envsDir'] + "/tpmCalculator.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/TPMCount_rna/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        cd {output.dir}
        TPMCalculator -g {input.annotation} -k gene_id -e -p -a  -b {input.sam}  
        """



rule gff3_to_gtf_old_Genome:
    input:
        gff3 = config['resources']  + "/old_genomes/annotation/Mus_musculus_{sid}_v1.109.gff3", 
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        gtf = temp(config['resultDir']  + "/annotationGTF/old_genome/{sid}_v3.2.gtf")
    threads:1
    priority:10
    conda: config['envsDir'] + "/tpmCalculator.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/gff3_to_gtf_rna_oldGenome/{sid}.log"
    shell:
        "agat_convert_sp_gff2gtf.pl --gff {input.gff3} -o {output.gtf} > {log} "



rule TPMCount_counts_rna_old_genome:
    input:
        sam = config['resultDir']  + "/STAR_rna_old_genome_strain_wise/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/annotationGTF/old_genome/{sid}_v3.2.gtf",
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        dir = config['resultDir']  + "/TPMCount_rna/oldGenome/{sid}-{tp}.{rep}/",
        tpm = protected(config['resultDir']  + "/TPMCount_rna/oldGenome/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out_genes.out")
    threads:1
    priority:10
    conda: config['envsDir'] + "/tpmCalculator.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/TPMCount_rna_oldGenome/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        cd {output.dir}
        TPMCalculator -g {input.annotation} -k gene_id -e -p -b {input.sam}  
        """
