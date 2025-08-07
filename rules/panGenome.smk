"""
copy fasta and gff file 
"""

rule copy_genome_for_panGenome:
    input:
        fasta = config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta",
        gtf = config['resources']  + "/annotation/{sid}_v3.2.gff3"
    output:
        genome_gff_dir  = dir(config['resultDir']  + "/pangenome/fasta/"),
        fasta = temp(config['resultDir']  + "/pangenome/fasta/{sid}.fa")
        gtf = temp(config['resultDir']  + "/pangenome/fasta/{sid}.gff")
    threads: 1
    params:
        extra = "--limitGenomeGenerateRAM  96000000000",
    priority:1
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/copy_pagenome/{sid}.log",
    shell:
        """
        cp {input.fasta} {output.fasta}
        cp {input.gtf} {input.gtf}
        """


rule panGenome_and_SV: 
    input:
        genome_gff_dir  = dir(config['resultDir']  + "pangenome/fasta/"),
        genome_list  = config['resources']  + "/strains_list.lst",
        script_dir = config['scriptDir'] + "/psvcp_v1.01/Construct_pan_and_Call_sv.py",  
    output:
        outdir = config['resultDir']  + "/pangenome"
        pan = protected(config['resultDir']  + "/pangenome/pan.fa")
    threads: 12
    params:
        prefix = 'mouse_strains',
        extra = "-t 12",
    priority:10
    conda: config['envsDir'] + "/pangenome.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/copy_pagenome/{sid}.log",
    shell:
        """
        cd {output.outdir}
        python3 {input.script_dir}/Construct_pan_and_Call_sv.py {params.extra} script_dir {input.script_dir} genome_dir {input.genome_gff_dir} genome_list {input.genome_list} -o {params.prefix}
        """

