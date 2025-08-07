rule sense_antisense_srna_counts_strains:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation =  config['resultDir']  + "/TEannotation/pacBio/{sid}_v3.gff3",
        script = config['scriptDir'] + "/sh/sense_antisense.sh",
        ref = config['strains_genomes_path']  + "/{sid}_chromosomes_MT.fasta",
    output:
        file = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/results.tsv",
        sense = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/sense.fasta",
        antisense = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/antisense.fasta",
        sense_bam = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/all_sense.bam",
        antisense_bam = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/all_antisense.bam",
    threads:
        1
    priority:10
    conda: 
        config['envsDir'] + "/sense_antisense.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        out = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/",
    log:
        config['logDir']  +"/sesne_antisense/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        "sh {input.script} {input.sam} {input.annotation} {input.ref} {params.out} > {log}"


rule feturecount_sense_TE_strains:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/TEannotation/pacBio/{sid}_v3.gff3",
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/sense_fetureCount/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:1
    priority:15
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['TE_rpkm']['featureCounts']  + " -s 1 "
    log:
        config['logDir']  +"/sense_antisense/fetureCount/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "0.72.0/bio/subread/featurecounts"


rule feturecount_antiSense_TE_strains:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",# list of sam or bam files
        annotation = config['resultDir']  + "/TEannotation/pacBio/{sid}_v3.gff3",
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        protected(multiext(config['resultDir']  + "/antisensefetureCount/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:1
    priority:15
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['TE_rpkm']['featureCounts']  + " -s 2 "
    log:
        config['logDir']  +"/antisense/fetureCount/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "0.72.0/bio/subread/featurecounts"



rule sense_ping_pong_z_score_srna_strains:
    localrule: True
    input:
        bam = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/all_sense.bam", # list of sam or bam files
        script = config['scriptDir'] + "/sRNAplot/signature_modified.py",
    output:
        hscore = config['resultDir']  + "/sense_pingpong/{sid}/{sid}-{tp}.{rep}/hscore.csv" ,
        zscore = protected(config['resultDir']  + "/sense_pingpong/{sid}/{sid}-{tp}.{rep}/zscore.csv")
    threads: 4
    priority:10
    conda:
        config["envsDir"] +  "/pingpong.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",  
    log:
        config['logDir']  +"/sense_pingpong/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        python {input.script} --input {output.bam} --minquery 25 --maxquery 33 --mintarget 25 --maxtarget 33 --minscope 1 --maxscope 21 --output_h {output.hscore}  --output_z {output.zscore}  > {log}
        """



rule antisense_ping_pong_z_score_srna_strains:
    localrule: True
    input:
        bam = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/allall_antisense.bam" ,# list of sam or bam files
        script = config['scriptDir'] + "/sRNAplot/signature_modified.py",
    output:
        hscore = config['resultDir']  + "/antisense_pingpong/{sid}/{sid}-{tp}.{rep}/hscore.csv" ,
        zscore = protected(config['resultDir']  + "/antisense_pingpong/{sid}/{sid}-{tp}.{rep}/zscore.csv")
    threads: 4
    priority:10
    conda:
        config["envsDir"] +  "/pingpong.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",  
    log:
        config['logDir']  +"/antisense_pingpong/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        python {input.script} --input {output.bam} --minquery 25 --maxquery 33 --mintarget 25 --maxtarget 33 --minscope 1 --maxscope 21 --output_h {output.hscore}  --output_z {output.zscore}  > {log}
        """



rule collapse_sense_seq_srna:
    localrule: True
    input:
        fasta = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/sense.fasta",
        scriptDir = config['scriptDir']
    output:
        config['resultDir']  + "/collapse/sense/{sid}-{tp}.{rep}.fasta.gz"
    log:
        stdout=config['logDir']  + "/collapse/sense/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/collapse/sense/{sid}-{tp}.{rep}.err"
    threads: 1
    priority:2
    resources:
        mem = 200000,
        queue  = 'normal',
    conda: config['envsDir'] + "/tstk.yaml"
    shell:
        "collapse {input.fasta} {output} > {log.stdout} 2> {log.stderr}"




rule strains_peterplots_sense_9th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/sense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/sens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/sens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/sens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/sens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/sens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/sens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.err"
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

rule strains_peterplots_sense_10th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/sense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/sens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/sens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/sens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/sens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/sens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/sens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.err"
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


rule strains_peterplots_sense_11th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/sense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/sens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/sens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/sens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/sens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/sens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/sens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.err"
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


rule strains_peterplots_sense:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/sense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/sens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/sens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/sens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/sens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/sens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/sens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.err"
    threads: 1
    priority:3
    conda: config['envsDir'] + "/tstk.yaml"
    resources:
        mem = 16000,
        queue  = 'normal',
    shell: r""" 
        python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.normal} --uncollapse --noN --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.normal_csv}> {log.stdout} 2> {log.stderr}
        python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.collapsed} --normnreads --noN --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.collapsed_csv} > {log.stdout} 2>> {log.stderr}
    """




rule collapse_antisense_seq_srna:
    localrule: True
    input:
        fasta = config['resultDir']  + "/sense_antisense/{sid}/{sid}-{tp}.{rep}/antisense.fasta",
        scriptDir = config['scriptDir']
    output:
        config['resultDir']  + "/collapse/antisense/{sid}-{tp}.{rep}.fasta.gz"
    log:
        stdout=config['logDir']  + "/collapse/antisense/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/collapse/antisense/{sid}-{tp}.{rep}.err"
    threads: 1
    priority:2
    resources:
        mem = 200000,
        queue  = 'normal',
    conda: config['envsDir'] + "/tstk.yaml"
    shell:
        "collapse {input.fasta} {output} > {log.stdout} 2> {log.stderr}"


rule strains_peterplots_antisense:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/antisense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/antisens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/antisens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/antisens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/antisens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/antisens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/antisens_seqlogo/peterplots/pacBio/{sid}-{tp}.{rep}.err"
    threads: 1
    priority:3
    conda: config['envsDir'] + "/tstk.yaml"
    resources:
        mem = 16000,
        queue  = 'normal',
    shell: r""" 
        python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.normal} --uncollapse --noN --normnreads --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.normal_csv}> {log.stdout} 2> {log.stderr}
        python {input.scriptDir}/tstk/tstk/peterplot.py {input.fasta} {output.collapsed} --normnreads --noN --title {wildcards.sid}-{wildcards.tp}.{wildcards.rep} --minlength 19 --maxlength 36 --outraw {output.collapsed_csv} > {log.stdout} 2>> {log.stderr}
    """



rule strains_peterplots_antisense_9th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/antisense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/antisens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/antisens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/antisens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/antisens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/antisens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/antisens_seqlogo/A9th/pacBio/{sid}-{tp}.{rep}.err"
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

rule strains_peterplots_antisense_10th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/antisense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/antisens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/antisens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/antisens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/antisens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/antisens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/antisens_seqlogo/A10th/pacBio/{sid}-{tp}.{rep}.err"
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


rule strains_peterplots_antisense_11th_logo:
    localrule: True
    input:
        fasta = config['resultDir']  + "/collapse/antisense/{sid}-{tp}.{rep}.fasta.gz",
        scriptDir = config['scriptDir']
    output:
        normal = protected(config['resultDir']  + "/antisens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.svg"),
        collapsed = protected(config['resultDir']  + "/antisens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.collapsed.svg"),
        normal_csv = protected(config['resultDir']  + "/antisens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.tab"),
        collapsed_csv = protected(config['resultDir']  + "/antisens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.collapsed.tab")
    log:
        stdout=config['logDir']  + "/antisens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.out",
        stderr=config['logDir']  + "/antisens_seqlogo/A11th/pacBio/{sid}-{tp}.{rep}.err"
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