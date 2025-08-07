
rule minimap2_and_paf2chain:
    input:
        query=config['strains_genomes_path'] + "/{sid}_chromosomes_MT.fasta",
        ref="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/black6/genome/Mus_musculus.GRCm39.dna.toplevel.fa"
    output:
        paf=temp(config['resultDir']  + "minimap2_chain_alignments/{sid}_to_GRCm39.paf"),
        chain=config['resultDir']  + "minimap2_chain_alignments/{sid}_to_GRCm39.chain"
    params:
        threads=16,
        preset="asm5"
    log:
        log_file = config['logDir'] + "/chain_files/{sid}_minimap_chain_file.log"
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    shell:
        """
        # Step 1: Generate PAF file using minimap2
        minimap2 -t {params.threads} -cx {params.preset} --cs {input.ref} {input.query} > {output.paf}

        # Step 2: Convert PAF to chain file using paf2chain
        paf2chain -i {output.paf} > {output.chain}
        """