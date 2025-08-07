rule ping_pong_z_score_srna_strains:
    localrule: True
    input:
        bam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        bai = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
        script = config['scriptDir'] + "/sRNAplot/signature_modified.py",
    output:
        hscore = config['resultDir']  + "/pingpong/{sid}/{sid}-{tp}.{rep}/hscore.csv",
        zscore = protected(config['resultDir']  + "/pingpong/{sid}/{sid}-{tp}.{rep}/zscore.csv")
    threads: 4
    priority:15
    conda:
        config["envsDir"] +  "/pingpong.yaml"
    resources:
        mem = 400000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",  
    log:
        config['logDir']  +"/pingpong/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        conda run -n snakemake python {input.script} --input {input.bam} --minquery 25 --maxquery 33 --mintarget 25 --maxtarget 33 --minscope 1 --maxscope 21 --output_h {output.hscore}  --output_z {output.zscore}  > {log}
        """

