"""
LiftOver genome coordinates and prepare BED files
"""

# Remove 'chr' prefix from chromosome names in BED files
rule remove_chr_prefix_genome:
    input:
        bed=config['resultDir'] + "/zamore/{genome}/piRNA_gene_annotation_{genome}_with_chr.bed"
    output:
        bed=config['resultDir'] + "/zamore/{genome}/piRNA_gene_annotation_{genome}_without_chr.bed"
    log:
        config['logDir'] + "/remove_chr_prefix/{genome}.log"
    shell:
        """
        sed 's/^chr//' {input.bed} > {output.bed} 2> {log}
        """

# Convert BED12 to BED6 and sort
rule bed12_to_bed6_genome:
    input:
        bed=config['resultDir'] + "/zamore/{genome}/piRNA_gene_annotation_{genome}_without_chr.bed"
    output:
        bed=config['resultDir'] + "/zamore/{genome}/piRNA_gene_annotation_{genome}_without_chr_sort_bed6.bed"
    log:
        config['logDir'] + "/bed12_to_bed6/{genome}.log"
    conda:
        config['envsDir'] + "/bedtools.yaml"
    shell:
        """
        bedtools bed12tobed6 -i {input.bed} | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """

# LiftOver from GRCm38 to GRCm39
rule liftover_GRCm38_to_GRCm39:
    input:
        bed=config['resultDir'] + "/zamore/GRCm38.68/piRNA_gene_annotation_GRCm38.68_with_chr.bed",
        chain=config['resultDir'] + "/liftOverChainFiles/GRCm38.68_to_GRCm39.106.over.chain"
    output:
        bed=config['resultDir'] + "/zamore/GRCm39.106/piRNA_gene_annotation_GRCm39.106_with_chr.bed",
        unmapped=config['resultDir'] + "/zamore/GRCm39.106/piRNA_gene_annotation_GRCm39.106_with_chr_unmapped.bed"
    log:
        config['logDir'] + "/liftover/GRCm38_to_GRCm39.log"
    benchmark:
        config['resultDir'] + "/benchmarks/liftover/GRCm38_to_GRCm39.txt"
    params:
        extra="-fudgeThick"
    resources:
        mem_mb=4000,
        queue='normal'
    conda:
        config['envsDir'] + "/liftover.yaml"
    shell:
        """
        liftOver {params.extra} {input.bed} {input.chain} {output.bed} {output.unmapped} &> {log}
        """

# General rule for LiftOver from GRCm39 to strains with varying minMatch
rule liftover_GRCm39_to_strains:
    input:
        bed=config['resultDir'] + "/zamore/GRCm39.106/piRNA_gene_annotation_GRCm39.106_without_chr.bed",
        chain=config['resultDir'] + "/liftOverChainFiles/GRCm39_chains/GRCm39_{sid}_chromosomes_MT_unplaced.chain"
    output:
        bed=config['resultDir'] + "/zamore/with_{minMatch}/{sid}/piRNA_gene_annotation_{sid}_without_chr.bed",
        unmapped=config['resultDir'] + "/zamore/with_{minMatch}/{sid}/piRNA_gene_annotation_{sid}_without_chr_unmapped.bed"
    params:
        extra=lambda wildcards: f"-fudgeThick -minMatch=0.{wildcards.minMatch}"
    log:
        config['logDir'] + "/liftover/with_{minMatch}/{sid}.log"
    benchmark:
        config['resultDir'] + "/benchmarks/with_{minMatch}/liftover/GRCm39_to_{sid}.txt"
    resources:
        mem_mb=4000,
        queue='normal'
    conda:
        config['envsDir'] + "/liftover.yaml"
    shell:
        """
        liftOver {params.extra} {input.bed} {input.chain} {output.bed} {output.unmapped} &> {log}
        """

# Convert BED12 to BED6 and sort for strains after LiftOver
rule bed12_to_bed6_liftover_strains:
    input:
        bed=config['resultDir'] + "/zamore/with_{minMatch}/{sid}/piRNA_gene_annotation_{sid}_without_chr.bed"
    output:
        bed=config['resultDir'] + "/zamore/with_{minMatch}/{sid}/piRNA_gene_annotation_{sid}_without_chr_sort_bed6.bed"
    log:
        config['logDir'] + "/bed12_to_bed6_liftover_strains/with_{minMatch}/{sid}.log"
    conda:
        config['envsDir'] + "/bedtools.yaml"
    shell:
        """
        bedtools bed12tobed6 -i {input.bed} | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """
