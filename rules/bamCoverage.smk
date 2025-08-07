rule convert_bam_to_bigwig_strain:
    input:
        bam = config['resultDir'] + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        bai = config['resultDir'] + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai"
    output:
        forward = config['resultDir'] + "/bamCoverageSrna/pacBio/{sid}/{sid}-{tp}.{rep}_plusStrand.bw",
        reverse1 = config['resultDir'] + "/bamCoverageSrna/pacBio/{sid}/{sid}-{tp}.{rep}_minusStrand.bw"
    conda: config['envsDir'] + "/bigwig.yaml"
    params:
        normlize = '--normalizeUsing CPM'
    threads: 2
    priority: 10
    resources:
        mem = 20000,  # Adjust based on your system's capabilities
        queue = 'normal',
        tmpdir = config['tmpDir']
    log:
        forward = config['logDir'] + "/bamCoverageSrnaBW/pacBio/{sid}/{sid}-{tp}.{rep}_forward.log",
        reverse1 = config['logDir'] + "/bamCoverageSrnaBW/pacBio/{sid}/{sid}-{tp}.{rep}_reverse.log" 
    shell:
        """
        bamCoverage {params.normlize} -b {input.bam} -o {output.forward} --filterRNAstrand forward > {log.forward} &&
        bamCoverage {params.normlize} -b {input.bam} -o {output.reverse1} --filterRNAstrand reverse > {log.reverse1}
        wait
        """
