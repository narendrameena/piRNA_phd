
rule downsample:
    localrule: True
    input:
        r1 = [config['resultDir']  + "/cutadapt_pe_rna/{{sid}}-{{tp}}.{{rep}}_r1.fastq".format(rep=rep) for rep in reps],
        r2 = [config['resultDir']  + "/cutadapt_pe_rna/{{sid}}-{{tp}}.{{rep}}_r2.fastq".format(rep=rep) for rep in reps]
    output:
        r1 = "results/downsample/{{sid}}-{{tp}}_r1.fastq.gz",
        r2 = "results/downsample/{{sid}}-{{tp}}_r2.fastq.gz"
    threads: 2
    log:
        stdout = "log/downsample/{{sid}}-{{tp}}.out",
        stderr = "log/downsample/{{sid}}-{{tp}}.err"
    resources:
        mem = 64000
    conda: "envs/seqtk.yaml"
    shell:
        """
        # Concatenate all replicates for r1 and downsample
        cat {input.r1} > results/downsample/tmp_{{sid}}-{{tp}}_r1.fastq.gz
        seqtk sample -s123 results/downsample/tmp_{{sid}}-{{tp}}_r1.fastq.gz 50000000 | gzip > {output.r1}
        
        # Concatenate all replicates for r2 and downsample
        cat {input.r2} > results/downsample/tmp_{{sid}}-{{tp}}_r2.fastq.gz
        seqtk sample -s123 results/downsample/tmp_{{sid}}-{{tp}}_r2.fastq.gz 50000000 | gzip > {output.r2}
        """

rule trinity:
    localrule: True
    input:
        r1 = config['resultDir'] + "/downsample/{sid}-{tp}_r1.fastq.gz",
        r2 = config['resultDir'] + "/downsample/{sid}-{tp}_r2.fastq.gz"
    output:
        f = config['resultDir'] + "/trinity/trinity-{sid}-{tp}.Trinity.fasta"
    threads: 12
    log:
        stdout = "log/trinity/{sid}-{tp}.out",
        stderr = "log/trinity/{sid}-{tp}.err"
    resources:
        mem = 128000
    params:
        d = "results/trinity/trinity-{sid}-{tp}"
    conda: "envs/trinity.yaml"
    shell:
        """
        echo "********************************" | tee {log.stdout} {log.stderr}
        echo "NEW RUN: $(date)" | tee -a {log.stdout} {log.stderr}
        echo "********************************" | tee -a {log.stdout} {log.stderr}
        Trinity --normalize_reads --min_kmer_cov 2 --seqType fq --left {input.r1} --right {input.r2} --output {params.d} --SS_lib_type RF --CPU {threads} --max_memory 90G >> {log.stdout} 2>> {log.stderr}
        echo "*******MOVE FASTA FILE**********" | tee -a {log.stdout} {log.stderr}
        mv {params.d}/Trinity.fasta {output.f}
        echo "*******CREATE EMPTY DIR*********" | tee -a {log.stdout} {log.stderr}
        mkdir {params.d}_empty >> {log.stdout} 2>> {log.stderr}
        echo "*******DELETE TRINITY FILES*********" | tee -a {log.stdout} {log.stderr}
        rsync -a --delete {params.d}_empty/ {params.d}/ >> {log.stdout} 2>> {log.stderr}
        echo "*******DELETE EMPTY DIRS*********" | tee -a {log.stdout} {log.stderr}
        rmdir {params.d} >> {log.stdout} 2>> {log.stderr}
        rmdir {params.d}_empty >> {log.stdout} 2>> {log.stderr}
        """