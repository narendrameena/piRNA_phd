

"""

Downloading chain files 

"""


rule get_liftover_chain_files_for_refrence_genomes:
    output:
       config['resultDir']  + "/liftOverChainFiles/GRCm38.68_to_GRCm39.106.over.chain"
    resources:
        mem = 4000,
        queue  = 'normal',
    threads: 1        
    params:
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"
    shell:"""
        wget -O {output}.gz {params.url}
        gunzip {output}.gz
    """
