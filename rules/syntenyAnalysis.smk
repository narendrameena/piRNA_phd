"""
genomic structures predicted by SyRI

"""

rule pac_bio_syntany_analysis_GRCm39_106:
    input:
        ref = config['resultDir']  + "/ref_genome/GRCm39.106.fasta",
        query = config['strains_genomes_path']  + "/C57BL_6NJ_chromosomes_MT.fasta",
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/GRCm39.106/GRCm39.106_C57BL_6NJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/GRCm39.106/GRCm39.106_C57BL_6NJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/GRCm39.106/")
    threads: 4
    params:
        extra = "",
        prefix = "GRCm39.106_C57BL_6NJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/GRCm39.106.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_C57BL_6NJ:
    input:
        ref = config['strains_genomes_path']  + "/C57BL_6NJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/NZO_HlLtJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/C57BL_6NJ/C57BL_6NJ_NZO_HlLtJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/C57BL_6NJ/C57BL_6NJ_NZO_HlLtJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/C57BL_6NJ/")
    threads: 4
    params:
        extra = "",
        prefix = "C57BL_6NJ_NZO_HlLtJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/C57BL_6NJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_NZO_HlLtJ:
    input:
        ref = config['strains_genomes_path']  + "/NZO_HlLtJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/NOD_ShiLtJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/NZO_HlLtJ/NZO_HlLtJ_NOD_ShiLtJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/NZO_HlLtJ/NZO_HlLtJ_NOD_ShiLtJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/NZO_HlLtJ/")
    threads: 4
    params:
        extra = "",
        prefix = "NZO_HlLtJ_NOD_ShiLtJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/NZO_HlLtJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_NOD_ShiLtJ:
    input:
        ref = config['strains_genomes_path']  + "/NOD_ShiLtJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/AKR_J_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/NOD_ShiLtJ/NOD_ShiLtJ_AKR_J.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/NOD_ShiLtJ/NOD_ShiLtJ_AKR_J_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/NZO_HlLtJ/")
    threads: 4
    params:
        extra = "",
        prefix = "NOD_ShiLtJ_AKR_J_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/NOD_ShiLtJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """



rule pac_bio_syntany_analysis_AKR_J:
    input:
        ref = config['strains_genomes_path']  + "/AKR_J_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/A_J_chromosomes_MT.fasta",
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/AKR_J/AKR_J_A_J.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/AKR_J/AKR_J_A_J_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/AKR_J/")
    threads: 4
    params:
        extra = "",
        prefix = "AKR_J_A_J_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/AKR_J.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_A_J:
    input:
        ref = config['strains_genomes_path']  + "/A_J_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/BALB_cJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/A_J/A_J_BALB_cJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/A_J/A_J_BALB_cJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/A_J/")
    threads: 4
    params:
        extra = "",
        prefix = "A_J_BALB_cJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/A_J.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_BALB_cJ:
    input:
        ref = config['strains_genomes_path']  + "/BALB_cJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/C3H_HeJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/BALB_cJ/BALB_cJ_C3H_HeJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/BALB_cJ/BALB_cJ_C3H_HeJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/BALB_cJ/")
    threads: 4
    params:
        extra = "",
        prefix = "BALB_cJ_C3H_HeJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/BALB_cJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_C3H_HeJ:
    input:
        ref = config['strains_genomes_path']  + "/C3H_HeJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/CBA_J_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/C3H_HeJ/C3H_HeJ_CBA_J.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/C3H_HeJ/C3H_HeJ_CBA_J_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/C3H_HeJ/")
    threads: 4
    params:
        extra = "",
        prefix = "C3H_HeJ_CBA_J_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/C3H_HeJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_CBA_J:
    input:
        ref = config['strains_genomes_path']  + "/CBA_J_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/DBA_2J_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/CBA_J/CBA_J_DBA_2J.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/CBA_J/CBA_J_DBA_2J_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/CBA_J/")
    threads: 4
    params:
        extra = "",
        prefix = "CBA_J_DBA_2J_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/CBA_J.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_DBA_2J:
    input:
        ref = config['strains_genomes_path']  + "/DBA_2J_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/129S1_SvImJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/DBA_2J/DBA_2J_129S1_SvImJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/DBA_2J/DBA_2J_129S1_SvImJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/DBA_2J/")
    threads: 4
    params:
        extra = "",
        prefix = "DBA_2J_129S1_SvImJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/DBA_2J.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """



rule pac_bio_syntany_analysis_129S1_SvImJ:
    input:
        ref = config['strains_genomes_path']  + "/129S1_SvImJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/LP_J_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/129S1_SvImJ/129S1_SvImJ_LP_J.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/129S1_SvImJ/129S1_SvImJ_LP_J_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/129S1_SvImJ/")
    threads: 4
    params:
        extra = "",
        prefix = "129S1_SvImJ_LP_J_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/129S1_SvImJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_LP_J:
    input:
        ref = config['strains_genomes_path']  + "/LP_J_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/FVB_NJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/LP_J/LP_J_FVB_NJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/LP_J/LP_J_FVB_NJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/LP_J/")
    threads: 4
    params:
        extra = "",
        prefix = "LP_J_FVB_NJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/LP_J.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_FVB_NJ:
    input:
        ref = config['strains_genomes_path']  + "/FVB_NJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/WSB_EiJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/LP_J/LP_J_WSB_EiJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/LP_J/LP_J_WSB_EiJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/FVB_NJ/")
    threads: 4
    params:
        extra = "",
        prefix = "LP_J_WSB_EiJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/FVB_NJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_WSB_EiJ:
    input:
        ref = config['strains_genomes_path']  + "/WSB_EiJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/CAST_EiJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/WSB_EiJ/WSB_EiJ_CAST_EiJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/WSB_EiJ/WSB_EiJ_CAST_EiJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/WSB_EiJ/")
    threads: 4
    params:
        extra = "",
        prefix = "WSB_EiJ_CAST_EiJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/WSB_EiJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_CAST_EiJ:
    input:
        ref = config['strains_genomes_path']  + "/CAST_EiJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/PWK_PhJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/CAST_EiJ/CAST_EiJ_PWK_PhJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/CAST_EiJ/CAST_EiJ_PWK_PhJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/CAST_EiJ/")
    threads: 4
    params:
        extra = "",
        prefix = "CAST_EiJ_PWK_PhJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/CAST_EiJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """


rule pac_bio_syntany_analysis_PWK_PhJ:
    input:
        ref = config['strains_genomes_path']  + "/PWK_PhJ_chromosomes_MT.fasta",
        query = config['strains_genomes_path']  + "/SPRET_EiJ_chromosomes_MT.fasta"
    output:
        sam = temp(config['resultDir']  + "/syntanyAnalysis/PWK_PhJ/PWK_PhJ_SPRET_EiJ.sam"),
        out = config['resultDir']  + "/syntanyAnalysis/PWK_PhJ/PWK_PhJ_SPRET_EiJ_syri.out",
        dir_out = directory(config['resultDir']  + "/syntanyAnalysis/PWK_PhJ/")
    threads: 4
    params:
        extra = "",
        prefix = "PWK_PhJ_SPRET_EiJ_"
    priority:10
    conda: 
        config['envsDir'] +  "/syntanyAnalysis.yaml"
    resources:
        mem = 200000,
        queue  = 'normal',
    log:
        config['logDir']  + "/syntanyAnalysis/PWK_PhJ.log",
    shell:
        """
        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.query} > {output.sam)
        syri -c {output.sam) -r {input.ref} -q {input.query} -k -F S --nc {threads} --prefix {params.prefix} --dir {output.dir_out}
        """



rule piRNA_naming_by_syntany:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/sh/piTNASynatanyOnlyOnoverlap.sh",
        bed = config['resultDir']  + "/minimap2/strains/{sid}/{sid}-{tp}.500.bed12.bed"
    output:
        output = directory(config['resultDir']  + "/piRNABasedOnSyntany/{sid}_{tp}/"),
        pirna = config['resultDir']  + "/piRNABasedOnSyntany/{sid}_{tp}/final_output.bed",
    threads:
        1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  + "/piRNABasedOnSyntany/{sid}_{tp}.log"
    shell:
        "bash {input.script} {input.bed} {input.annotation} {output.output} > {log}"



rule utr_gene_name:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/sh/utr3GeneList.sh",
        fetureFile = config['resultDir']  + "/IntronExonfeatureCount/pacBio/list3/{sid}/{sid}-{tp}.{rep}.featureCounts",
    output:
        utr= config['resultDir']  + "/UTR_3_gene/{sid}/{sid}-{tp}.{rep}_gene_list.txt",
    threads:
        1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  + "/piRNABasedOnSyntany/UTR_3_gene/{sid}/{sid}-{tp}.{rep}_gene_list.log"
    shell:
        "bash {input.script} {input.fetureFile} {input.annotation} {output.utr} > {log}"


rule gene_name_from_list1:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/sh/GenefromList1.sh",
        fetureFile = config['resultDir']  + "/IntronExonfeatureCount/pacBio/list1/{sid}/{sid}-{tp}.{rep}.featureCounts",
    output:
        utr= config['resultDir']  + "/genefromlist1/{sid}/{sid}-{tp}.{rep}_gene_list.txt",
    threads:
        1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/piRNABasedOnSyntany/gene_3_gene/{sid}/{sid}-{tp}.{rep}_gene_list.log"
    shell:
        "bash {input.script} {input.fetureFile} {input.annotation} {output.utr} > {log}"


rule gene_name_based_syntany:
    localrule: True
    input:
        annotation = config['resources']  + "/annotation/{sid}_v3.2.gff3",
        script = config['scriptDir'] + "/sh/piRNASuntanyBasedNaming.sh",
        piRNA =  config['resultDir'] + "/filtred_final_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}.gtf"
    output:
        syntany = config['resultDir']  + "/syntany/{sid}/{sid}-{tp}.{rep}.txt",
    threads:    1
    priority:16
    conda: config['envsDir'] + "/bedtools.yaml"
    resources:
        mem = 20000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
    log:
        config['logDir']  +"/syntany/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        "bash {input.script} {input.annotation} {input.piRNA} {output.syntany} > {log}"