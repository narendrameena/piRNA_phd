rule bed_precursors_annotation_ref_GRCm38_68:
    input:
        config['resultDir']  + "/minimap2/GRCm38.68/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        config['resultDir']  + "/minimap2/GRCm38.68/annotation/{sid}/{sid}-{tp}.500.bed12.bed"
    threads: 1
    params:
        "--extended --ambiguities best_all"
    conda: config['envsDir'] + "/bedannotation.yaml"
    resources:
        queue  = 'normal',
        mem = 200000,
        tmpdir= config['tmpDir']
    log:
       config['logDir']  + "/minimap2_precursor/GRCm38.68/annotation/{sid}-{tp}.bed.log"
    shell: 
        "bed_annotation  {input} -g mm10 -o {output}  {params} > {log}"
