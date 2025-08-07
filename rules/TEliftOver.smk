rule TE_liftover_GRCHm38toGRCHm39:
    input:
        bed = ancient(config['resources']  + "/TE_Annotation/TE_all_del_insr_final.bed"),
        chain = config['resultDir']  + "/liftOverChainFiles/GRCm38.68_to_GRCm39.106.over.chain",
    output:
        bed =  protected(config['resultDir']  + "/TE_Ins_Del/GRCm39.106/TE_all_del_insr_with_chr.bed"),
        unmapped = config['resultDir']  + "/TE_Ins_Del/GRCm39.106/TE_all_del_insr_with_chr_unmapped.bed",
    log:
        config['logDir']  + "/TE_Ins_Del/GRCm39.106/GRCHm38toGRCHm39.log"
    benchmark:
        config['resultDir'] +"/benchmarks/TE_Ins_Del/GRCm39.106/GRCHm38toGRCHm39.txt",
    params:
        extra = "-fudgeThick",
    resources:
        mem = 4000,
        queue  = 'normal',
    conda:
        config["envsDir"] +  "/liftover.yaml"
    shell:
        "liftOver  {params.extra} {input.bed} {input.chain} {output.bed} {output.unmapped} &> {log} "

rule TE_rmove_chr:
    input:
        config['resultDir']  + "/TE_Ins_Del/GRCm39.106/TE_all_del_insr_with_chr.bed",
    output:
        config['resultDir']  + "/TE_Ins_Del/GRCm39.106/TE_all_del_insr_without_chr.bed",
    resources:
        mem = 4000,
        queue  = 'normal',
    shell:"""
        cut -c 4- {input} > {output}
    """

rule TE_liftover_GRCHm39_to_strains_with_09:
    input:
        bed = config['resultDir']  + "/TE_Ins_Del/GRCm39.106/TE_all_del_insr_without_chr.bed",
        chain = config['resultDir']  + "/liftOverChainFiles/GRCm39_chains/GRCm39_{sid}_chromosomes_MT_unplaced.chain",
    output:
        bed =  protected(config['resultDir']  + "/zamore/with_09/{sid}/TE_all_del_insr_{sid}_without_chr.bed"), 
        unmapped = config['resultDir']  + "/TE_Ins_Del/with_09/{sid}/TE_all_del_insr_{sid}_without_chr_unmapped.bed",
    log:
        config['logDir']  + "/TE_Ins_Del/with_09GRCHm39_to_{sid}.log"
    benchmark:
        config['resultDir'] +"/benchmarks/TE_Ins_Del/with_09/liftover/GRCHm39_to_{sid}.txt"
    params:
        extra = "-fudgeThick",
    resources:
        mem = 4000,
        queue  = 'normal',
    conda:
        config["envsDir"] +  "/liftover.yaml"
    shell:
        "liftOver  {params.extra}  {input.bed} {input.chain} {output.bed} {output.unmapped} &> {log} "


rule TE_liftover_GRCHm39_to_strains_with_02:
    input:
        bed = config['resultDir']  + "/TE_Ins_Del/GRCm39.106/TE_all_del_insr_without_chr.bed",
        chain = config['resultDir']  + "/liftOverChainFiles/GRCm39_chains/GRCm39_{sid}_chromosomes_MT_unplaced.chain",
    output:
        bed =  protected(config['resultDir']  + "/TE_Ins_Del/with_02/{sid}/TE_all_del_insr_{sid}_without_chr.bed"), 
        unmapped = config['resultDir']  + "/TE_Ins_Del/with_02/{sid}/TE_all_del_insr_{sid}_without_chr_unmapped.bed",
    log:
        config['logDir']  + "/TE_Ins_Del/with_02/GRCHm39_to_{sid}.log"
    benchmark:
        config['resultDir'] +"/benchmarks/TE_Ins_Del/with_02/liftover/GRCHm39_to_{sid}.txt"
    params:
        extra = "-fudgeThick -minMatch=0.2",
    resources:
        mem = 4000,
        queue  = 'normal',
    conda:
        config["envsDir"] +  "/liftover.yaml"
    shell:
        "liftOver  {params.extra}  {input.bed} {input.chain} {output.bed} {output.unmapped} &> {log} "


rule TE_liftover_GRCHm39_to_strains_with_05:
    input:
        bed = config['resultDir']  + "/TE_Ins_Del/GRCm39.106/TE_all_del_insr_without_chr.bed",
        chain = config['resultDir']  + "/liftOverChainFiles/GRCm39_chains/GRCm39_{sid}_chromosomes_MT_unplaced.chain",
    output:
        bed =  protected(config['resultDir']  + "/TE_Ins_Del/with_05/{sid}/TE_all_del_insr_{sid}_without_chr.bed"), 
        unmapped = config['resultDir']  + "/TE_Ins_Del/with_05/{sid}/TE_all_del_insr_{sid}_without_chr_unmapped.bed",
    log:
        config['logDir']  + "/TE_Ins_Del/with_05/GRCHm39_to_{sid}.log"
    benchmark:
        config['resultDir'] +"/benchmarks/TE_Ins_Del/with_05/liftover/GRCHm39_to_{sid}.txt"
    params:
        extra = "-fudgeThick -minMatch=0.5",
    resources:
        mem = 4000,
        queue  = 'normal',
    conda:
        config["envsDir"] +  "/liftover.yaml"
    shell:
        "liftOver  {params.extra}  {input.bed} {input.chain} {output.bed} {output.unmapped} &> {log} "
