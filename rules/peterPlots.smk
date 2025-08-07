



"""

Collapsed sRNA raw sequncing data usig Tom's tool tstk

"""

rule collapse_raw_seq_srna:
    localrule: True
    input:
        fasta = config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq",
        scriptDir = config['scriptDir']
    output:
        temp(config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.raw.fasta.gz")
    log:
        stdout=config['logDir']  + "/collapse/{sid}-{tp}.{rep}.raw.out",
        stderr=config['logDir']  + "/collapse/{sid}-{tp}.{rep}.raw.err"
    threads: 1
    priority:2
    resources:
        mem = 200000,
        queue  = 'normal',
    conda: config['envsDir'] + "/tstk.yaml"
    shell:
        "collapse {input.fasta} {output} > {log.stdout} 2> {log.stderr}"



"""

collapsed alinged BAM files using Tom's Scripts tstk package

"""

rule strains_collapse_aligned_bam_srna:
    localrule: True
    input:
        bam=config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        scriptDir = config['scriptDir']
    output:
        temp(config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.aligned.fasta.gz")
    log:
        stdout=config['logDir']  + "/collapse/{sid}-{tp}.{rep}.aligned.out",
        stderr=config['logDir']  + "/collapse/{sid}-{tp}.{rep}.aligned.err"
    threads: 1
    priority:2
    resources:
        mem = 200000,
        queue  = 'normal',
    conda: config['envsDir'] + "/tstk.yaml"
    shell:
        "collapse  {input.bam} {output} > {log.stdout} 2> {log.stderr}"
    


"""

plot Pterplots and raw output of Peterplot using Tom's script; tstk package 

"""

rule strains_peterplots_srna:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.{rawali}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/peterplots/pacBio/{sid}-{tp}.{rep}.{rawali}.svg"),
        collapsed = protected(config['resultDir']  + "/peterplots/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.svg"),
        normal_csv = protected(config['resultDir']  + "/peterplots/pacBio/{sid}-{tp}.{rep}.{rawali}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/peterplots/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.tab")
    log:
        stdout=config['logDir']  + "/peterplots/pacBio/{sid}-{tp}.{rep}.{rawali}.out",
        stderr=config['logDir']  + "/peterplots/pacBio/{sid}-{tp}.{rep}.{rawali}.err"
    threads: 1
    priority:3
    conda: config['envsDir'] + "/tstk.yaml"
    resources:
        mem = 16000,
        queue  = 'normal',
    shell: r""" 
        conda run -n snakemake python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.normal} --uncollapse --noN --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.normal_csv}> {log.stdout} 2> {log.stderr}
        conda run -n snakemake python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.collapsed} --normnreads --noN --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.collapsed_csv} > {log.stdout} 2>> {log.stderr}
    """


rule strains_peterplots_10th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.{rawali}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.{rawali}.svg"),
        collapsed = protected(config['resultDir']  + "/seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.svg"),
        normal_csv = protected(config['resultDir']  + "/seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.{rawali}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.tab")
    log:
        stdout=config['logDir']  + "/seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.{rawali}.out",
        stderr=config['logDir']  + "/seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.{rawali}.err"
    threads: 1
    priority:3
    conda: config['envsDir'] + "/tstk.yaml"
    resources:
        mem = 16000,
        queue  = 'normal',
    shell: r""" 
        python {input.scriptDir}/python/peterplot_10th_nt.py {input.fasta} {output.normal} --uncollapse --noN --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 5 --maxlength 36 --outraw {output.normal_csv}> {log.stdout} 2> {log.stderr}
        python {input.scriptDir}/python/peterplot_10th_nt.py  {input.fasta} {output.collapsed} --normnreads --noN --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 5 --maxlength 36 --outraw {output.collapsed_csv} > {log.stdout} 2>> {log.stderr}
    """



rule strains_peterplots_9th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.{rawali}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.{rawali}.svg"),
        collapsed = protected(config['resultDir']  + "/seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.svg"),
        normal_csv = protected(config['resultDir']  + "/seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.{rawali}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.tab")
    log:
        stdout=config['logDir']  + "/seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.{rawali}.out",
        stderr=config['logDir']  + "/seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.{rawali}.err"
    threads: 1
    priority:3
    conda: config['envsDir'] + "/tstk.yaml"
    resources:
        mem = 16000,
        queue  = 'normal',
    shell: r""" 
        python {input.scriptDir}/python/peterplot_9th_nt.py {input.fasta} {output.normal} --uncollapse --noN --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 5 --maxlength 36 --outraw {output.normal_csv}> {log.stdout} 2> {log.stderr}
        python {input.scriptDir}/python/peterplot_9th_nt.py  {input.fasta} {output.collapsed} --normnreads --noN --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 5 --maxlength 36 --outraw {output.collapsed_csv} > {log.stdout} 2>> {log.stderr}
    """

rule strains_peterplots_11th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/{sid}-{tp}.{rep}.{rawali}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.{rawali}.svg"),
        collapsed = protected(config['resultDir']  + "/seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.svg"),
        normal_csv = protected(config['resultDir']  + "/seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.{rawali}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.collapsed.{rawali}.tab")
    log:
        stdout=config['logDir']  + "/seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.{rawali}.out",
        stderr=config['logDir']  + "/seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.{rawali}.err"
    threads: 1
    priority:3
    conda: config['envsDir'] + "/tstk.yaml"
    resources:
        mem = 16000,
        queue  = 'normal',
    shell: r""" 
        python {input.scriptDir}/python/peterplot_11th_nt.py {input.fasta} {output.normal} --uncollapse --noN --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 5 --maxlength 36 --outraw {output.normal_csv}> {log.stdout} 2> {log.stderr}
        python {input.scriptDir}/python/peterplot_11th_nt.py  {input.fasta} {output.collapsed} --normnreads --noN --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 5 --maxlength 36 --outraw {output.collapsed_csv} > {log.stdout} 2>> {log.stderr}
    """

rule old_Genome_collapse_aligned_bam_srna:
    input:
        bam=config['resultDir']  + "/STAR_srna_strain_wise_old_genomes/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        scriptDir = config['scriptDir']
    output:
        temp(config['resultDir']  + "/collapse/old_genome/{sid}-{tp}.{rep}.aligned.fasta.gz")
    log:
        stdout=config['logDir']  + "/collapse/old_genome/{sid}-{tp}.{rep}.aligned.out",
        stderr=config['logDir']  + "/collapse/old_genome/{sid}-{tp}.{rep}.aligned.err"
    threads: 1
    priority:2
    resources:
        mem = 200000,
        queue  = 'normal',
    conda: config['envsDir'] + "/tstk.yaml"
    shell:
        "collapse  {input.bam} {output} > {log.stdout} 2> {log.stderr}"
    


"""

plot Pterplots and raw output of Peterplot using Tom's script; tstk package 

"""

rule old_genome_peterplots_srna:
    input:
        fasta = config['resultDir']  + "/collapse/old_genome/{sid}-{tp}.{rep}.aligned.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/peterplots/old_genome/{sid}-{tp}.{rep}.{rawali}.svg"),
        collapsed = protected(config['resultDir']  + "/peterplots/old_genome/{sid}-{tp}.{rep}.collapsed.{rawali}.svg"),
        normal_csv = protected(config['resultDir']  + "/peterplots/old_genome/{sid}-{tp}.{rep}.{rawali}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/peterplots/old_genome/{sid}-{tp}.{rep}.collapsed.{rawali}.tab")
    log:
        stdout=config['logDir']  + "/peterplots/old_genome/{sid}-{tp}.{rep}.{rawali}.out",
        stderr=config['logDir']  + "/peterplots/old_genome/{sid}-{tp}.{rep}.{rawali}.err"
    threads: 1
    priority:3
    conda: config['envsDir'] + "/tstk.yaml"
    resources:
        mem = 16000,
        queue  = 'normal',
    shell: r""" 
        python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.normal} --uncollapse  --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.normal_csv}> {log.stdout} 2> {log.stderr}
        python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.collapsed} --normnreads  --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.collapsed_csv} > {log.stdout} 2>> {log.stderr}
    """