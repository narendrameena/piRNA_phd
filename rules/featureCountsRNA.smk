"""
RNA Seq Feature Count with PacBio Assemblies Using Snakemake
"""

# RNA Strains Feature Count Rule
rule feature_counts_rna_strains:
    localrule: True
    input:
        samples=lambda wildcards: config['resultDir'] + f"/STAR_rna_strain_wise/{wildcards.sid}-{wildcards.tp}.{wildcards.rep}/Aligned.sortedByCoord.out.bam",
        annotation=lambda wildcards: config['resultDir'] + f"/strains_genome_annotation/{wildcards.sid}.gff3",
    output:
        protected(multiext(config['resultDir'] + "/featureCount_rna/strains/{sid}-{tp}.{rep}",
                           ".featureCounts",
                           ".featureCounts.summary",
                           ".featureCounts.jcounts"))
    threads: 2
    priority: 10
    resources:
        mem=20000,
        queue='normal',
        tmpdir=config['tmpDir']
    params:
        extra=config['params']['rna']['featureCounts']
    log:
        config['logDir'] + "/featureCount_rna/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"

# Transposable Elements Annotation Conversion Rule
rule TE_repeatmasker_to_gff3_strains:
    input:
        repeatmasker=lambda wildcards: config['resultDir'] + f"/ref_genome_masked/{wildcards.sid}_chromosomes_MT_unplaced.fasta.out",
        script=config['scriptDir'] + "/sh/rm2gff3.sh",
    output:
        annotation=config['resultDir'] + "/TEannotation/pacBio/{sid}.gff3",
    threads: 1
    priority: 9
    resources:
        mem=20000,
        queue='normal',
        tmpdir=config['tmpDir']
    shell:
        "sh {input.script} {input.repeatmasker} > {output.annotation}"

# sRNA TE Strains Feature Count Rule
rule featurecount_srna_TE_strains:
    input:
        samples=lambda wildcards: config['resultDir'] + f"/STAR_srna_strain_wise/{wildcards.sid}/{wildcards.sid}-{wildcards.tp}.{wildcards.rep}/Aligned.sortedByCoord.out.bam",
        annotation=lambda wildcards: config['resultDir'] + f"/TEannotation/pacBio/{wildcards.sid}.gff3",
    output:
        protected(multiext(config['resultDir'] + "/TEtranscriptCount/pacBio/{sid}/{sid}-{tp}.{rep}",
                           ".featureCounts",
                           ".featureCounts.summary",
                           ".featureCounts.jcounts"))
    threads: 1
    priority: 10
    resources:
        mem=20000,
        queue='normal',
        tmpdir=config['tmpDir']
    params:
        extra=config['params']['TE_rpkm']['featureCounts']
    log:
        config['logDir'] + "/TEtranscriptCount/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"


