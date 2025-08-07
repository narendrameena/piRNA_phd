"""
Generate bed6 files from RepeatMasker output files
"""

rule bed6_from_repeatmasker_Files:
    input:
        repeatmasker_out=config['resources'] + "/repeatMasker/{sid}_chromosomes_MT_unplaced.fasta.out",
        script=config['scriptDir'] + "/sh/repeatbedFile.awk"
    output:
        config['resultDir'] + "/repeatMaskerBedFiles/{sid}.bed6"
    log:
        config['logDir'] + "/repeatMaskerBedFiles/{sid}.log"
    resources:
        mem=1000,
        queue='normal',
        tmpdir=config['tmpDir']
    shell:
        """
        awk -f "{input.script}" "{input.repeatmasker_out}" > "{output}" 2> "{log}"
        """

"""
Intersect piRNA and TE elements from Trinity minimap alignments
"""

rule bedtools_intersection_piRNA_repeatMasker_trinity_minimap2:
    input:
        left=config['resultDir'] + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed.gz",
        right=config['resultDir'] + "/repeatMaskerBedFiles/{sid}.bed6"
    output:
        config['resultDir'] + "/repeatMaskerBedFiles/piRNA_intersection_minimap/{sid}-{tp}.{rep}_intersected.bed.gz"
    params:
        extra="-wo -s"
    log:
        config['logDir'] + "/repeatMaskerBedFiles/piRNA_intersection_minimap/{sid}-{tp}.{rep}_trinity.log"
    resources:
        mem=4000,
        queue='normal',
        tmpdir=config['tmpDir']
    wrapper:
        "v5.8.2/bio/bedtools/intersect"

rule bedtools_intersection_piRNA_repeatMasker_strains:
    input:
        left=config['resultDir'] + "/bedFileFromSrnaBam/{sid}/{sid}-{tp}.{rep}.bed12_seq.bed.gz",
        right=config['resultDir'] + "/repeatMaskerBedFiles/{sid}.bed6"
    output:
        config['resultDir'] + "/repeatMaskerBedFiles/piRNA_intersection_strains/{sid}-{tp}.{rep}_with_seq_intersected.bed.gz"
    params:
        extra="-wo -s"
    log:
        config['logDir'] + "/repeatMaskerBedFiles/piRNA_intersection/{sid}-{tp}.{rep}_with_seq_strains.log"
    resources:
        mem=4000,
        queue='normal',
        tmpdir=config['tmpDir']
    wrapper:
        "v5.8.2/bio/bedtools/intersect"
